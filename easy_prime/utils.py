import warnings
warnings.filterwarnings("ignore")
import os 
import argparse
import shutil
import datetime
import getpass
import uuid
import mechanize
import re
import pandas as pd
import yaml
import itertools
import pickle
import subprocess
from skbio.alignment import global_pairwise_align_nucleotide
from skbio import DNA
import RNA
from copy import deepcopy as dp
import numpy as np
import sys
import math
#------------------ FILE IO ---------------------------

def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()
	
def print_parameters(myDict):
	myGroup = {}
	myGroup['Easy-Prime'] = ['genome_fasta','scaffold','n_jobs','debug','PE2_model','PE3_model','extend_length']
	myGroup['PBS searching'] = ['min_PBS_length','max_PBS_length']
	myGroup['RTT searching'] = ['min_RTT_length','max_RTT_length','min_distance_RTT5','max_max_RTT_length']
	myGroup['sgRNA searching'] = ['gRNA_search_space','sgRNA_length','offset','PAM','max_target_to_sgRNA','max_max_target_to_sgRNA']
	myGroup['ngRNA searching'] = ['max_ngRNA_distance']
	for k in myGroup:
		print_group(myDict,myGroup[k],k)


def print_group(myDict,myList,group_title):
	print ("-------- Parameter Group: %s --------"%(group_title))
	for l in myList:
		print ("%s: %s"%(l,myDict[l]))

def get_parameters(config):
	p_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
	# return dict
	parameters = {}
	# default parameters
	pre_defined_list = {}
	#------------ EASY-PRIME related-----------
	pre_defined_list["genome_fasta"] = "/home/yli11/Data/Human/hg19/fasta/hg19.fa"
	pre_defined_list["n_jobs"] = -1
	pre_defined_list["scaffold"] = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
	pre_defined_list["debug"] = 0
	pre_defined_list["extend_length"] = 1000 # extracting +- 1000bp center at target pos from the genome, in 99.9% cases, you don't need to change this. If change to less than 500, will trigger fasta input mode, may cause error.
	pre_defined_list["PE2_model"] = p_dir+"../model/PE2_model_final.py"
	pre_defined_list["PE3_model"] = p_dir+"../model/PE3_model_final.py"
	
	#------------ PBS -----------
	pre_defined_list["min_PBS_length"] = 10
	pre_defined_list["max_PBS_length"] = 15
	
	#------------ RTT -----------
	pre_defined_list["min_RTT_length"] = 10
	pre_defined_list["max_RTT_length"] = 20 # if no candidate is found, this value will be increased by 5, max to max_max_RTT_length
	pre_defined_list["max_max_RTT_length"] = 50 
	pre_defined_list["min_distance_RTT5"] = 5

	#------------ sgRNA -----------
	pre_defined_list["gRNA_search_space"] = 200
	pre_defined_list["sgRNA_length"] = 20 
	pre_defined_list["offset"] = -3
	pre_defined_list["PAM"] = "NGG"
	pre_defined_list["max_target_to_sgRNA"] = 10 # if no candidate is found, this value will be increased by 5, max to max_max_target_to_sgRNA
	pre_defined_list["max_max_target_to_sgRNA"] = 30
	
	#------------ ngRNA ------------
	pre_defined_list["max_ngRNA_distance"] = 100 # if no candidate is found, this value will be increased by 20, max to max_max_ngRNA_distance
	pre_defined_list["max_max_ngRNA_distance"] = 200 
	pre_defined_list["search_iteration"] = 1 # not affect anything
	
	try:
		with open(config, 'r') as f:
			manifest_data = yaml.load(f,Loader=yaml.FullLoader)
	except:
		print ("Config data is not provided or not parsed successfully, Default parameters were used.")
		
	for p in pre_defined_list:
		try:
			parameters[p] = manifest_data[p]	
		except:
			parameters[p] = pre_defined_list[p]
	return parameters

def to_bed3(chr,start,end):
	outfile = str(uuid.uuid4()).split("-")[-1]
	write_file(outfile,"\t".join([chr,str(start),str(end)]))
	return outfile

def write_fasta(file_name,myDict):
	out = open(file_name,"wt")
	for k in myDict:
		out.write(">"+k+"\n")
		out.write(myDict[k]+"\n")
	out.close()



#------------------ sgRNA finder ---------------------------
from Bio import SeqUtils
def run_pam_finder(target_fa,seq,PAM,abs_start_pos,chr):
	
	# SeqUtils.nt_search("AGGCGGGGG", "NGG")
	# SeqUtils.nt_search("CCACCA", "NGG")
	# forward
	rev_seq = revcomp(target_fa)
	fwd_search = SeqUtils.nt_search(target_fa, seq+PAM)
	rev_search = SeqUtils.nt_search(rev_seq, seq+PAM)
	out = []
	if len(fwd_search) > 1:
		for s in fwd_search[1:]:
			# out.append([chr,s+abs_start_pos,s+abs_start_pos+len(seq),target_fa[s:(s+len(seq))],".","+"])
			out.append([chr,s+abs_start_pos,s+abs_start_pos+len(seq),target_fa[s:(s+len(seq))],target_fa[s:(s+len(seq)+len(PAM))],"+"])
	if len(rev_search) > 1:
		for s in rev_search[1:]:
			# out.append([chr,(len(target_fa)-s)+abs_start_pos-len(seq),(len(target_fa)-s)+abs_start_pos,rev_seq[s:(s+len(seq))],".","-"])
			out.append([chr,(len(target_fa)-s)+abs_start_pos-len(seq),(len(target_fa)-s)+abs_start_pos,rev_seq[s:(s+len(seq))],rev_seq[s:(s+len(seq)+len(PAM))],"-"])
	return pd.DataFrame(out)

def force_recommend_dPAM_PE3b(r,max_eff):
	rank = 0
	if "dPAM" in r.index:
		rank += 1
	if "PE3b" in r.index:
		if not max_eff-r.predicted_efficiency>max_eff*0.1:
			rank += 1
	return rank

#------------------ Fasta Operators ---------------------------
from Bio import SeqIO

def find_pos_ref_alt(x,y):
	"""find variant pos, ref and alt given ref(x) and alt(y) sequences

	Method

	trim 5 end, then trim 3 end, what left is the ref and alt, pos the trim 5 position

	"""

	for i in range(len(x)):
		if x[i]!=y[i]:
			new_x = x[i:]
			new_y = y[i:]
			pos = i+1
			break
	# print (new_x,new_y)
	for i in range(min(len(new_x),len(new_y))):
		i=i+1
		# print ("second for",i,new_x[-i],new_y[-i])
		if new_x[-i]!=new_y[-i]:
			ref = new_x[:-i+1]
			alt = new_y[:-i+1]
			break
		else:

			ref = new_x[:-i]
			alt = new_y[:-i]
			
	## check
	check_y = x[:pos-1]+alt+x[pos+len(ref)-1:]
	# print ("debug",pos,i,ref)
	# print ("new y",x[:pos-1],x[pos+len(ref)-1:],alt)
	if y!= check_y:
		print ("something is wrong, please fix it")
		print (x)
		print (y)
		return -1,ref,alt
	return pos,ref,alt

def fasta2vcf(f):
	"""convert fasta to vcf dataframe

	Input
	-----

	Fasta file, _ref is recognized as ref and _alt is used as alt, these are two keywords

	Output
	------

	vcf dataframe: chr, pos, name, ref, alt, reference sequence


	"""
	my_dict = {}
	for r in SeqIO.parse(f, "fasta"):
		my_dict[r.id] = str(r.seq).upper()
	print (my_dict)
	vcf = pd.DataFrame()
	index_list = []
	chr_list = []
	pos_list = []
	ref_list = []
	alt_list = []
	seq_list = []
	for k in my_dict:
		if not "_ref" in k:
			continue
		name = k.replace("_ref","")
		if not name+"_alt" in my_dict:
			print (k,"alt sequence not found. Please use _ref and _alt keywords. Skip...")
			continue
		ref_seq,alt_seq = my_dict[k],my_dict[name+"_alt"]
		if len(ref_seq) < 30:
			print (k,"Please input sequence length at least 30bp. Skip...")
			continue
		if ref_seq == alt_seq:
			print (k,"Ref and Alt sequence is the same. Please check. Skip...")
			continue
		pos,ref,alt = find_pos_ref_alt(ref_seq,alt_seq)
		index_list.append(name)
		chr_list.append(k)
		seq_list.append(ref_seq)
		pos_list.append(pos)
		ref_list.append(ref)
		alt_list.append(alt)
	vcf[0] = chr_list
	vcf[1] = pos_list
	vcf[2] = index_list
	vcf[3] = ref_list
	vcf[4] = alt_list
	vcf[5] = seq_list
	vcf = vcf[vcf[1]!=-1]
	if vcf.shape[0] == 0:
		print ("no valid sequences in:",f)
		print ("Exit...")
		sys.exit(1)

	return vcf	

def vcf2fasta(vcf,extend_length=None,genome_fasta=None,**kwargs):
	"""extracting +- extend_length given vcf target mutation
	
	input
	-------
	
	vcf is a dataframe, chr, pos, id, ref, alt
	
	return
	--------
	
	a list of sequences, same order
	
	"""
	
	out_bed = str(uuid.uuid4()).split("-")[-1]+".bed"
	out_fa = out_bed+".tab"
	df = vcf.copy()
	df['chr'] = df[0]
	df['start'] = df[1]-extend_length
	df['end'] = df[1]+extend_length
	df[['chr','start','end']].to_csv(out_bed,sep="\t",header=False,index=False)

	p1 = subprocess.Popen(['bedtools','getfasta','-fi',genome_fasta,'-bed',out_bed,'-fo',out_fa,'-tab'],bufsize=0)
	p1.communicate()
	df = pd.read_csv(out_fa,sep="\t",header=None,index_col=0)

	os.remove(out_bed)
	os.remove(out_fa)

	return [x.upper() for x in df[1]]

def get_opposite_strand(x):
	if x == "+":
		return "-"
	return "+"


def sub_fasta_single(target_fa,target_pos, abs_start,abs_end):
	"""given the target_fa we extracted from target_pos, we get sub fasta
	
	target_fa: we extended +- N bp of the target set, this length is 2N
	
	user sequence, abs start and end
	
	target_pos is 1 index, the target pos that we used to get target_fa
	
	Assumption and user query, target_fa on the same chr
	
	"""

	N = int(len(target_fa)/2)
	start =  N-(target_pos-abs_start)
	end = abs_end - abs_start + start
	seq = target_fa[start:end]
	if len(target_fa)<1000:
		## user input fasta:
		seq = target_fa[abs_start:abs_end]
	return seq



def get_fasta_simple(target_fa,df, target_pos,strand=False):
	"""save time and memory get fasta
	
	target_fa: we extended +- N bp of the target set, this length is 2N
	
	df is a normal bed file
	
	target_pos is 1 index
	
	df and target_pos, all on the same chr
	
	"""
	temp = df.copy()
	# print ("len target_fa",len(target_fa))
	N = int(len(target_fa)/2)
	temp.columns = list(range(len(df.columns)))
	temp.index = temp[0]+":"+temp[1].astype(str)+"-"+temp[2].astype(str)
	if len(target_fa)<1000:
		## user input fasta:
		seq_list = []
		for r in temp.values.tolist():
			seq = target_fa[r[1]:r[2]]
			if strand:
				if r[5] == "-":
					seq = revcomp(seq)
			seq_list.append(seq)
		temp[3] = seq_list	
		return temp	
	temp[2] = temp[2]-temp[1]
	temp[1] = N-(target_pos-temp[1])
	temp[2] = temp[2]+temp[1]
	seq_list = []
	for r in temp.values.tolist():
		seq = target_fa[r[1]:r[2]]
		if strand:
			if r[5] == "-":
				seq = revcomp(seq)
		seq_list.append(seq)
	temp[3] = seq_list
	return temp


	
tab = str.maketrans("ACTG", "TGAC")
def revcomp(seq):
	return seq.translate(tab)[::-1]
		

	
#------------------ pegRNA Operators ---------------------------

	
def distance_matrix(lines):
	"""given sgRNA bed dataframe (the 5th is the name), get gRNA distance matrix
	
	index: chr_start_end_strand_seq (5th column)
	
	use the last element as the pos to calculate distance
	
	comparing to D Liu distance, this is always 1 larger than them, because we use the 4th nucleotide as the cut position, cas9 cut between 3rd and 4th
	
	return
	-------
	
	2d dict
	
	"""
	
	dist_dict={}
	for x in lines:
		dist_dict[x[4]]={}
		for y in lines:
			dist_dict[x[4]][y[4]] = x[-1]-y[-1]
	return dist_dict



def is_gRNA_valid(cas9_cut_position,target_mutation,strand,user_target_mutation_pos,diff):
	# cas9_cut_position=[chr,pos], pos is 1-based position, same as in vcf file
	# target_mutation=[chr,pos], pos is 1-based position, same as in vcf file
	"""
	PBS sequence can't be on the same side to the target site in terms of the cut site
	
	The position of the target mutation should be:
	On the right of the cas9 cut position if gRNA strand is +
	On the left of the cas9 cut position if gRNA strand is -

	user_target_mutation_pos, mutation correct cause bug when strand = -
	
	Return
	------
	
	is gRNA valid, target mutation distance to cut site

	
	"""
	if cas9_cut_position[0] != target_mutation[0]:
		return -1
	distance = int(target_mutation[1]-cas9_cut_position[1])
	# print (cas9_cut_position,target_mutation,strand,user_target_mutation_pos,diff,distance)
	if strand=="+":
		if distance>=0:
			return distance
	if strand=="-":
		if distance == 0:
			if user_target_mutation_pos != target_mutation[1]:
				distance = -1
		if distance<=0:
			distance = -distance
			if diff > 0:
				distance  = distance - diff + 1
			return distance	
	return -1
	
	


def get_gRNA_cut_site(start,end,strand,offset=-3):
	# return 1-index pos
	# cut site = the first base before cut near PAM direction

	if strand == "+":
		return int(end + offset)
	if strand == "-":
		return int(start - offset +1)
		


#------------------ featurize ---------------------------

	
def GC_content(seq):
	GC=["G","C"]
	count=0
	for i in seq:
		if i in GC:
			count+=1
	return float(count)/len(seq)
	pass
	


ulength = 31
md = RNA.md()
md.max_bp_span = 70
md.window_size = 70

def call_RNAplfold(seq,scaffold_length):
	"""RNA fold python binding
	"""
	# scaffold - RTT - PBS
	# we only care about RTT+PBS (3' extension) pairs with scaffold

	data = []
	fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
	fc.probs_window(ulength, RNA.PROBS_WINDOW_BPP | RNA.PROBS_WINDOW_UP, pf_window_callback2, data)
	
	# parse output
	df = pd.DataFrame(data)
	df2 = df.copy()
	df2[0] = df[1]
	df2[1] = df[0]
	df = pd.concat([df,df2])
	seq_length = 10
	RTT_start = scaffold_length + 1
	RTT_end = scaffold_length + seq_length
	scaffold_start = 1
	scaffold_end =  scaffold_length
	
	## subset df
	df = df[df[0]>=RTT_start]
	df = df[df[0]<=RTT_end]
	df = df[df[1]>=scaffold_start]
	df = df[df[1]<=scaffold_end]
	df = df[[0,2]]
	df = df.groupby(0).max()
	myDict = df[2].to_dict()
	
	value_list = []
	for i in range(seq_length):
		RTT_pos = scaffold_length +i+ 1
		if RTT_pos in myDict:
			value_list.append(myDict[RTT_pos])
		else:
			value_list.append(0)
	return value_list


def pf_window_callback2(v, v_size, i, maxsize, what, data=None):
	if what & RNA.PROBS_WINDOW_UP:
		pass
	else:
		data+=[[i,j,p] for j, p in enumerate(v) if (p is not None) and (p >= 0.01)]

def local_alignments(ref,q):
	query = StripedSmithWaterman(ref)
	alignment = query(q)
	return alignment['optimal_alignment_score']


def alignments_to_cigar(ref,q):
	Match = 0
	Mis = 0
	D = 0
	I = 0
	for i in range(len(ref)):
		a=ref[i]
		b=q[i]
		if a=="-":
			I=I+1
		elif b=="-":
			D=D+1
		elif a==b:
			Match+=1
		else:
			Mis+=1
	return [Match,Mis,D,I]	

def global_alignments(ref,q):

	s1 = DNA(ref)
	s2 = DNA(q)
	alignment, score, start_end_positions = global_pairwise_align_nucleotide(s1,s2,match_score=4,mismatch_score=1)
	return alignments_to_cigar(alignment[0]._string.decode("utf-8"),alignment[1]._string.decode("utf-8"))

def is_dPAM(PAM_seq, RTT, cut_offset=-3):
	# Assuming no N is RTT, which should be true
	# match PAM seq to RTT, should be abs(cut_offset)
	# print (PAM_seq, RTT)
	# will need to do revcomp no matter what, because RTT is always xxxxxxxPAM

	seq = revcomp(RTT)
	fwd_search = SeqUtils.nt_search(seq, PAM_seq)
	flag = 1
	if len(fwd_search) > 1:
		if abs(cut_offset) in fwd_search:
			flag = 0

	return flag	
	
def target_to_RTT5_feature(pegRNA,nick_gRNA,target_loc,RTS_length,alt_length):
	# 1. target mutation distance to cut 1 (pegRNA)
	# 2. target mutation distance to cut 2 (nick-gRNA)
	# target_loc [chr,pos]
	# 3. cut 1 to cut 2
	cut1 = get_gRNA_cut_site(pegRNA[1],pegRNA[2],pegRNA[3])
	## if nick_gRNA not exist return a big number
	if nick_gRNA[0] == "0":
		cut2=0
	else:
		cut2 = get_gRNA_cut_site(nick_gRNA[1],nick_gRNA[2],nick_gRNA[3])
	
	# cut2 to cut1
	a=cut2-cut1-1


	# target to cut1
	b=target_loc[1]-cut1
	
	if pegRNA[3]=="-":
		a+=2 # match to coordinate system
		a=-a
		b=-b
	c=cut2-target_loc[1]	
	if nick_gRNA[0] == "0":
		a=np.nan	
		c=np.nan
	if b <0:
		print ("pegRNA to target is less than 0!")
		exit()
	d = RTS_length - alt_length- b + 1 ## number of nucleotide to the RTT 5' end, from the target mutation (not including)
	if d <0:
		print ("pegRNA to target is less than 0!")
		exit()	
	index_list = ["nick_to_pegRNA","target_to_pegRNA","target_to_ngRNA","target_to_RTT5"]
	out = pd.DataFrame([a,b,c,d])
	out.index = index_list
	return out
	

	pass


#------------------ DeepSpCas9 score  ---------------------------


def list_to_fasta(l):
	out = []
	for i in l:
		if len(i) != 23:
			print ("something is wrong")
		out.append(">%s"%(i[:20]))
		out.append("AAAA%sAAA"%(i))
	return "\n".join(out)
def parse_webpage(c):
	a=re.findall('"([^"]*)"', str(c))
	# get the random id
	for x in a:
		# random id format could be changed during version updates
		if "DeepSpCas9/job/user" in x:
			random_id = x.split("/")[-1]
	url = "http://deepcrispr.info/DeepSpCas9/data/%s/Results.zip"%(random_id)
	df = pd.read_csv(url,sep="\t")
	df[df.ID==df['Guide Sequence (20bp)']]
	df.index = df.ID.tolist()
	return df['DeepSpCas9 Score'].to_dict()
	
def get_DeepSpCas9_score(gRNA_list):
	"""Grab score from web server
	
	input gRNA list should include PAM sequence, PAM is NGG
	"""
	url="http://deepcrispr.info/DeepSpCas9/"
	br = mechanize.Browser()
	br.set_handle_robots(False) # ignore robots
	br.open(url)
	br.select_form(nr=0)
	br["ENTER_FASTA"] = list_to_fasta(gRNA_list)
	res = br.submit()
	output = res.read()
	res = parse_webpage(output)
	flag = False
	for i in gRNA_list:
		if not i[:20] in res:
			print ("gRNA: %s NOT FOUND!"%(i[:20]))
			print ("DeepSpCas9 API error!")
			flag = True
			res[i[:20]] = 0
	if flag:
		print (output)
		print (res)
	return res
# gRNA_list = ['GGAATCCCTTCTGCAGCACCAGG','GGCCCAGACTGAGCACGTGAAGG']
# res = get_DeepSpCas9_score(gRNA_list)
# GGAATCCCTTCTGCAGCACC    55.750
# GGCCCAGACTGAGCACGTGA    52.897
# Name: DeepSpCas9 Score, dtype: float64
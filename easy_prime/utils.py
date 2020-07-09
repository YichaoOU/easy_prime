import warnings
warnings.filterwarnings("ignore")
import os 
import argparse
import shutil
import datetime
import getpass
import uuid
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

#------------------ FILE IO ---------------------------

def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()
	
def print_parameters(myDict):
	myGroup = {}
	myGroup['Easy-Prime'] = ['genome_fasta','scaffold','n_jobs','debug','ML_model','extend_length']
	myGroup['PBS searching'] = ['min_PBS_length','max_PBS_length']
	myGroup['RTT searching'] = ['min_RTT_length','max_RTT_length','min_distance_RTT5']
	myGroup['sgRNA searching'] = ['gRNA_search_space','sgRNA_length','offset','PAM']
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
	pre_defined_list["extend_length"] = 1000
	pre_defined_list["ML_model"] = p_dir+"../model/xgb_model_final.pkl"
	
	#------------ PBS -----------
	pre_defined_list["min_PBS_length"] = 8
	pre_defined_list["max_PBS_length"] = 17
	
	#------------ RTT -----------
	pre_defined_list["min_RTT_length"] = 10
	pre_defined_list["max_RTT_length"] = 25
	pre_defined_list["min_distance_RTT5"] = 3

	#------------ sgRNA -----------
	pre_defined_list["gRNA_search_space"] = 150
	pre_defined_list["sgRNA_length"] = 20
	pre_defined_list["offset"] = -3
	pre_defined_list["PAM"] = "NGG"
	pre_defined_list["max_target_to_sgRNA"] = 10
	
	#------------ ngRNA ------------
	pre_defined_list["max_ngRNA_distance"] = 100
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
			out.append([chr,s+abs_start_pos,s+abs_start_pos+len(seq),target_fa[s:(s+len(seq))],".","+"])
	if len(rev_search) > 1:
		for s in rev_search[1:]:
			out.append([chr,(len(target_fa)-s)+abs_start_pos-len(seq),(len(target_fa)-s)+abs_start_pos,rev_seq[s:(s+len(seq))],".","-"])
	return pd.DataFrame(out)


#------------------ Fasta Operators ---------------------------

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



def is_gRNA_valid(cas9_cut_position,target_mutation,strand):
	# cas9_cut_position=[chr,pos], pos is 1-based position, same as in vcf file
	# target_mutation=[chr,pos], pos is 1-based position, same as in vcf file
	"""
	PBS sequence can't be on the same side to the target site in terms of the cut site
	
	The position of the target mutation should be:
	On the right of the cas9 cut position if gRNA strand is +
	On the left of the cas9 cut position if gRNA strand is -
	
	Return
	------
	
	is gRNA valid, target mutation distance to cut site

	
	"""
	if cas9_cut_position[0] != target_mutation[0]:
		return -1
	distance = int(target_mutation[1]-cas9_cut_position[1])
	if strand=="+":
		if distance>=0:
			return distance
	if strand=="-":
		if distance<=0:
			return -distance	
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
	seq_length = 14
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



def is_dPAM(PAM_seq, target_pos,ref,alt,pegRNA_loc):
	"""wheter target mutation affects PAM sequence
	
	pegRNA_loc is [chr.start,end,strand]
	
	"""
	# currently only work for NG or NGG, and SNV
	# report a bug in Biopython
	# https://github.com/biopython/biopython/issues/3023
	# currently, assuming ref do not contain IUPAC letter
	# get PAM location
	flag = 0
	# print ("pegRNA_loc---",list(pegRNA_loc))
	if len(ref)==len(alt)==1:
		PAM_rel_pos = list(range(len(PAM_seq)))
		PAM_abs_pos = []
		
		if pegRNA_loc[3] == "-":
			for i in range(len(PAM_seq)):
				PAM_abs_pos.append(pegRNA_loc[1]-i)
		else:
			for i in range(len(PAM_seq)):
				PAM_abs_pos.append(pegRNA_loc[2]+1+i)
		if target_pos in PAM_abs_pos:
			for i in PAM_rel_pos:
				PAM_pos = PAM_abs_pos[i]
				if target_pos == PAM_pos:
					PAM_nuc = PAM_seq[i]
					if PAM_nuc != "N":
						flag = 1
						if ref != PAM_nuc:
							print ("PAM something is wrong")
			

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
	a=cut2-cut1


	# target to cut1
	b=target_loc[1]-cut1
	
	if pegRNA[3]=="-":
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


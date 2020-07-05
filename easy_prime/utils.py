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

#------------------ FILE IO ---------------------------

def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()
	
def print_parameters(myDict):
	myGroup = {}
	myGroup['Easy-Prime'] = ['genome_fasta','scaffold','n_jobs','debug','ML_model']
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
	pre_defined_list["ML_model"] = p_dir+"../model/xgb_model_final.pkl"
	
	#------------ PBS -----------
	pre_defined_list["min_PBS_length"] = 7
	pre_defined_list["max_PBS_length"] = 18
	
	#------------ RTT -----------
	pre_defined_list["min_RTT_length"] = 7
	pre_defined_list["max_RTT_length"] = 40
	pre_defined_list["min_distance_RTT5"] = 5

	#------------ sgRNA -----------
	pre_defined_list["gRNA_search_space"] = 150
	pre_defined_list["sgRNA_length"] = 20
	pre_defined_list["offset"] = -3
	pre_defined_list["PAM"] = "NGG"
	
	#------------ ngRNA ------------
	pre_defined_list["max_ngRNA_distance"] = 150
	
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

def delete_files(myList):
	for f in myList:
		# os.system("rm %s"%(f))
		subprocess.call("rm %s"%(f),shell=True)

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


def get_opposite_strand(x):
	if x == "+":
		return "-"
	return "+"

def get_fasta_given_bed(genome_fa,extended_file):
	out = extended_file+".for_cas_input.fa"
	command = "bedtools getfasta -fi %s -bed %s -fo %s"%(genome_fa,extended_file,out)
	# os.system(command)
	subprocess.call(command,shell=True)
	return out

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


def get_fasta_single(chr,start,end,genome_fasta=None):
	
	out_bed = str(uuid.uuid4()).split("-")[-1]
	out_fa = str(uuid.uuid4()).split("-")[-1]
	write_file(out_bed,"%s\t%s\t%s"%(chr,start,end))
	command = "bedtools getfasta -fi %s -bed %s -fo %s -tab"%(genome_fasta,out_bed,out_fa)
	# os.system(command)
	subprocess.call(command,shell=True)
	lines = open(out_fa).readlines()[0]
	seq = lines.split()[-1]
	# os.system("rm %s;rm %s"%(out_bed,out_fa))
	subprocess.call("rm %s;rm %s"%(out_bed,out_fa),shell=True)
	return seq
	

	
tab = str.maketrans("ACTG", "TGAC")
def revcomp(seq):
	return seq.translate(tab)[::-1]
		
def get_seq_from_bed(bed_file,genome_fasta):
	temp_file = addon_string = str(uuid.uuid4()).split("-")[-1]
	command = "bedtools getfasta -fi %s -bed %s -fo %s -tab -s -name"%(genome_fasta,bed_file,temp_file)
	# os.system(command)
	subprocess.call(command,shell=True)
	df = pd.read_csv(temp_file,header=None,sep="\t",index_col=0)
	
	# df.index = [x.split("::")[0] for x in df.index.tolist()]
	# print (df.head())
	# os.system("rm %s"%(temp_file))
	subprocess.call("rm %s"%(temp_file),shell=True)
	return df
	
	
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

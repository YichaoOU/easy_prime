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

import inspect
import itertools
import pickle
import subprocess

#------------------ FILE IO ---------------------------

def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()

def get_parameters(config):
	p_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
	parameters = {}
	pre_defined_list = {}
	pre_defined_list["genome_fasta"] = "/home/yli11/Data/Human/hg19/fasta/hg19.fa"
	pre_defined_list["gRNA_search_space"] = 500
	pre_defined_list["n_jobs"] = -1
	pre_defined_list["min_PBS_length"] = 8
	pre_defined_list["max_PBS_length"] = 18
	pre_defined_list["min_RTS_length"] = 9
	pre_defined_list["max_RTS_length"] = 60
	pre_defined_list["PBS_length_step"] = 2
	pre_defined_list["RTS_length_step"] = 3
	pre_defined_list["N_nick_gRNA_per_sgRNA"] = 5
	pre_defined_list["N_combinations_for_optimization"] = 800
	pre_defined_list["N_top_pegRNAs"] = 3 ## top pegRNAs per sgRNA
	pre_defined_list["scaffold"] = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
	## if you have large deletion in your file, for example
	## 80bp deletion, you should set this parameter > 80, say 100

	pre_defined_list["search_iteration"] = 1
	pre_defined_list["sgRNA_length"] = 20
	pre_defined_list["offset"] = -3
	pre_defined_list["debug"] = 0
	pre_defined_list["PAM"] = "NGG"
	# pre_defined_list["model1"] = "/home/yli11/Tools/easy_prime/model/model1_final_0.4.pkl"
	pre_defined_list["model1"] = p_dir+"../model/model1_final_0.4.pkl"
	pre_defined_list["w1"] = 0.4
	pre_defined_list["model2"] = p_dir+"../model/model2_final_0.6.pkl"
	pre_defined_list["w2"] = 0.6
	pre_defined_list["PCA_model"] = p_dir+"../model/3mer_PCA.pkl"


	try:
		with open(config, 'r') as f:
			manifest_data = yaml.load(f,Loader=yaml.FullLoader)
		# print (manifest_data)
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
	
	
#------------------ CasOFFinder ---------------------------

def run_casOFFinder(genome_fasta,PAM,your_seq_list,nMisMatch=0):
	cas_input = str(uuid.uuid4()).split("-")[-1]
	cas_output = str(uuid.uuid4()).split("-")[-1]
	pattern = "N"*len(your_seq_list[0])+PAM
	config = [genome_fasta,pattern]
	for i in your_seq_list:
		config.append(i+PAM+" %s"%(nMisMatch))
	write_file(cas_input,"\n".join(config))
	command = "cas-offinder %s C %s > /dev/null 2>&1;rm %s"%(cas_input,cas_output,cas_input)
	subprocess.call(command, shell=True)	
	# os.system(command)
	return cas_output


def cas_local_to_df(cas_output,PAM,sgRNA_length):
	"""convert chr1:185056572-185056972	374	TCTATACAAGCTTTTGTGGTAGG	+	0 to bed"""
	df = pd.read_csv(cas_output,sep="\t",header=None)
	
	df = df.dropna()
	df['chr'] = [x.split(":")[0] for x in df[1]]
	df['start'] = df.apply(lambda x:cas_local_row_apply(x,len(PAM)),axis=1)
	df['end'] = df['start']+sgRNA_length
	df['seq'] = [x[:-len(PAM)] for x in df[3].tolist()]
	df[5] = "."
	return df[['chr','start','end','seq',5,4]]


def cas_local_row_apply(x,PAM_length):
	chr,temp = x[1].split(":")
	start,_ = temp.split("-")
	start = int(start)
	if x[4] == "-":
		start = start+PAM_length+x[2]
	else:
		start = start+x[2]
	return start


#------------------ Fasta Operators ---------------------------


def get_opposite_strand(x):
	if x == "+":
		return "-"
	return "+"

def get_fasta_given_bed(genome_fa,extended_file):
	out = extended_file+".fa"
	command = "bedtools getfasta -fi %s -bed %s -fo %s"%(genome_fa,extended_file,out)
	# os.system(command)
	subprocess.call(command,shell=True)
	return out


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

	
def distance_matrix(lines,offset):
	# pos to dict
	# chr_start_end_seq is the key
	dist_dict={}
	for x in lines:
		dist_dict[x[4]]={}
		for y in lines:
			dist_dict[x[4]][y[4]] = abs(x[-1]-y[-1])

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


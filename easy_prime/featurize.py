from .utils import *
import pandas as pd
import glob
from skbio.alignment import StripedSmithWaterman
from skbio.alignment import global_pairwise_align_nucleotide
from skbio import DNA
import uuid
import itertools
import argparse
import os
import pickle
import numpy as np

"""

part of easy_prime


convert raw(x) to X

standard format is 7 columns

8th columns if avail, it's Y

any additional columns will be ignored

no need header

sample_ID,type,seq,chr,start,end,strand,editing_frequency


"""

"""
sample_ID,type,seq,chr,start,end,strand,editing_frequency
EDFIG9B_HEXA_21_REP2,pegRNA,ATCCTTCCAGTCAGGGCCAT,chr15,72638898,72638918,+,8.76e-06
EDFIG9B_HEXAS_2_REP2,pegRNA,ATCCTTCCAGTCAGGGCCAT,chr15,72638898,72638918,+,9.21e-06
EDFIG9B_HEXAS_2_REP3,pegRNA,ATCCTTCCAGTCAGGGCCAT,chr15,72638898,72638918,+,4.59e-05

"""

	
def GC_content(seq):
	GC=["G","C"]
	count=0
	for i in seq:
		if i in GC:
			count+=1
	return float(count)/len(seq)
	pass
	
def parse_RNAplfold_output(inFile,pegRNA_length,scaffold_length):
	keyword = "%start of base pair probability data"
	data = []
	flag = False
	with open(inFile) as f:
		for line in f:
			line = line.strip()
			if flag:
				if "ubox" in line:
					line = line.split()
					# print (line)
					data.append([int(line[0]),int(line[1]),float(line[2])])
					data.append([int(line[1]),int(line[0]),float(line[2])])
			if keyword == line:
				flag = True
	df = pd.DataFrame(data)
	# print (df.head())
	pegRNA = df[df[0]<=pegRNA_length]
	pegRNA = pegRNA[pegRNA[1]>pegRNA_length]
	# pegRNA_PBS_RTS = pegRNA[pegRNA[1]>pegRNA_length+scaffold_length]
	pegRNA_scaffold = pegRNA[pegRNA[1]<=pegRNA_length+scaffold_length]
	# if pegRNA_PBS_RTS.shape[0] == 0:
		# pegRNA_PBS_RTS_mean = 0
	# else:
		# pegRNA_PBS_RTS_mean = pegRNA_PBS_RTS[2].mean()
	pegRNA_PBS_RTS_mean = 0.5
	# max average
	pegRNA_scaffold_max_mean = pegRNA_scaffold.groupby(0).max()[2].mean()
	# percent
	pegRNA_scaffold = pegRNA_scaffold[pegRNA_scaffold[2]>pegRNA_PBS_RTS_mean]
	pegRNA_scaffold = pegRNA_scaffold.drop_duplicates(1)
	pegRNA_scaffold_percent_above_mean = pegRNA_scaffold.shape[0]/scaffold_length
	

	PBS_RTS = df[df[0]>pegRNA_length+scaffold_length]
	PBS_RTS = PBS_RTS[PBS_RTS[1]<=pegRNA_length+scaffold_length]
	PBS_RTS_scaffold = PBS_RTS[PBS_RTS[1]>pegRNA_length]

	# max average
	PBS_RTS_scaffold_max_mean = PBS_RTS_scaffold.groupby(0).max()[2].mean()
	# percent
	PBS_RTS_scaffold = PBS_RTS_scaffold[PBS_RTS_scaffold[2]>pegRNA_PBS_RTS_mean]
	PBS_RTS_scaffold = PBS_RTS_scaffold.drop_duplicates(1)
	PBS_RTS_scaffold_percent_above_mean = PBS_RTS_scaffold.shape[0]/scaffold_length
	
	return [pegRNA_scaffold_max_mean,pegRNA_scaffold_percent_above_mean,PBS_RTS_scaffold_max_mean,PBS_RTS_scaffold_percent_above_mean]
def sequence_kmer(k,seq):
	"""
	
	sequence is not long, use k=3 or k=4
	"""
	bases=['A','T','G','C']
	k_mer = {''.join(p):0 for p in itertools.product(bases, repeat=k)}
	for i in range(0,len(seq)-k+1):
		# print (i,i+k)
		k_mer[seq[i:i+k]]+=1
	return pd.DataFrame.from_dict(k_mer,orient="index")


def call_RNAplfold(pegRNA,PBS,RTS,scaffold):
	
	# to fasta
	total_sequence = pegRNA+scaffold+RTS+PBS
	total_sequence = total_sequence.replace("T","U")
	uid = str(uuid.uuid4()).split("-")[-1]
	write_fasta("%s.fa"%(uid),{uid:total_sequence})
	# run RNAplfold
	command = "RNAplfold < %s.fa"%(uid)
	# os.system(command)
	subprocess.call(command,shell=True)
	
	
	# parse file
	inFile = "%s_dp.ps"%(uid)
	
	
	# pegRNA/PBS-RTS percent of bp with pairing prob above mean
	# pegRNA/PBS-RTS max avearge
	out = pd.DataFrame(parse_RNAplfold_output(inFile,len(pegRNA),len(scaffold)))
	out.index = ['sgRNA_max_mean','sgRNA_percent','PBS_RTS_max_mean','PBS_RTS_percent']
	# os.system("rm %s*"%(uid))
	subprocess.call("rm %s*"%(uid),shell=True)
	
	
	
	return out
	

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


	
def distance_features(pegRNA,nick_gRNA,target_loc):
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
	if nick_gRNA[0] == "0":
		a=np.nan	
	index_list = ["nick_to_sgRNA","target_to_sgRNA"]
	out = pd.DataFrame([a,b])
	out.index = index_list
	return out
	

	pass
from copy import deepcopy as dp
def sequence_features(pegRNA,PBS,RTS,nick_gRNA,ref,alt,target_loc,seq_kmer=None,ref_alt_global_alignment=None):


	out = seq_kmer

	ref_alt_cigar = dp(ref_alt_global_alignment)

	# GC content
	pegRNA_GC = GC_content(pegRNA)
	PBS_GC = GC_content(PBS)
	RTS_GC = GC_content(RTS)
	if nick_gRNA == "0":
		nick_gRNA_GC = np.nan
	else:
		nick_gRNA_GC = GC_content(nick_gRNA)
	
	ref_alt_cigar += [pegRNA_GC,PBS_GC,RTS_GC,nick_gRNA_GC]
	# DNA length

	PBS_length = len(PBS)
	RTS_length = len(RTS)
	RTS_PBS_length_ratio = float(RTS_length)/PBS_length

	ref_alt_cigar += [PBS_length,RTS_length,RTS_PBS_length_ratio]

	columns = ["PC1","PC2","PC3","PC4","PC5","aln_ref_alt_match","aln_ref_alt_mis","aln_ref_alt_del","aln_ref_alt_ins","sgRNA_GC","PBS_GC","RTS_GC","nick_gRNA_GC","PBS_length","RTS_length","RTS_PBS_length_ratio"]

	out = pd.DataFrame(out+ref_alt_cigar)

	out.index = columns
	
	return out

def get_X_feature(s,d,seq_kmer=None,ref_alt = None,scaffold=None,**kwargs):
	# print (s)
	
	# Distance (2) , kmer (5), GC content (4), sequence length (3), RNA fold features(4)

	# for s,d in df.groupby("sample_ID")
	"""
	pegRNA_loc, nick_gRNA_loc, [chr,start,end,strand]
	target_loc [chr,pos]
	others are sequences
	"""
	
	## genomic positions
	loc_cols=['chr','start','end','strand']
	pegRNA_loc = d[d['type']=="sgRNA"][loc_cols].iloc[0].tolist()
	nick_gRNA_loc = d[d['type']=="nick-gRNA"][loc_cols].iloc[0].tolist()
	target_loc = d[['CHROM','POS']].iloc[0].tolist()
	
	## sequences
	pegRNA = d[d['type']=="sgRNA"]['seq'].tolist()[0]
	PBS = d[d['type']=="PBS"]['seq'].tolist()[0]
	RTS = d[d['type']=="RTS"]['seq'].tolist()[0]
	nick_gRNA = d[d['type']=="nick-gRNA"]['seq'].tolist()[0]
	
	## variants
	ref = d['REF'].tolist()[0]
	alt = d['ALT'].tolist()[0]
	### make sure when you specific ref and alt, it should be minimal
	###	e.g., if the variant is a CTT insertion, then ref can be G, and alt can be GCTT
	###	not ref: GGG alt: GGGCTT.
	### otherwise the featurize.py program will output an incorrect distance.
	
	if len(ref) != len(alt):
		target_loc[1] = target_loc[1]+1
	sequence_df = sequence_features(pegRNA,PBS,RTS,nick_gRNA,ref,alt,target_loc,seq_kmer,ref_alt)
	distance_df = distance_features(pegRNA_loc,nick_gRNA_loc,target_loc)
	

	energy_df = call_RNAplfold(pegRNA,PBS,RTS,scaffold)
	
	out = pd.concat([distance_df,sequence_df,energy_df])
	out.columns=[s]
	return out
	




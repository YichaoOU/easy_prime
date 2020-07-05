
# exec(open("vis_pegRNA.py").read())
from dna_features_viewer import GraphicFeature, GraphicRecord
import pandas as pd
import uuid
import os
import subprocess
import matplotlib.pyplot as plt
def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()

def get_fasta_single(chr,start,end,genome_fasta="/home/yli11/Data/Human/hg19/fasta/hg19.fa"):
	
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
	
npg_colors = ["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#464d4f"]
my_colors = {}
my_colors['sgRNA'] = npg_colors[0]
my_colors['PBS'] = npg_colors[1]
my_colors['RTS'] = npg_colors[2]
my_colors['nick-gRNA'] = npg_colors[3]
my_colors['variant'] = "#e6fc3f"

import sys
rawX = "test.top_pegRNAs.csv"

df = pd.read_csv(rawX,index_col=0)
pegRNA_id = df.index.tolist()[0]
variant_id = pegRNA_id.split("_")[0]
df = df.loc[pegRNA_id]

chr = df['CHROM'].tolist()[0]
start = df['start'].min()
start -= start%10-1
end = df['end'].max()
variant_pos = df.POS.min()
ref = df.REF.tolist()[0]
alt = df.ALT.tolist()[0]
predicted_efficiency = df.predicted_efficiency.tolist()[0]
pos = df.POS.min()-start
sequence = get_fasta_single(chr,start,end)
sequence = sequence.upper()

def get_strand(x):
	if x == "+":
		return 1
	else:
		return -1
feature_list = []
for s,r in df.iterrows():
	r_start = r.start-start
	r_end = r_start+(r.end-r.start)
	r_strand = get_strand(r.strand)
	gf = GraphicFeature(start=r_start, end=r_end, strand=r_strand, 
		color=my_colors[r.type],label=r.type)
	feature_list.append(gf)
	
record = GraphicRecord(sequence=sequence, features=feature_list)

ax, _ = record.plot(figure_width=int(len(sequence)/5))
record.plot_sequence(ax)
ax.fill_between((pos-1.5, pos-0.5), +1000, -1000, alpha=0.5,color=my_colors['variant'])
locs, labels = plt.xticks()
new_labels = []
flag = True
for i in locs:
	if flag:
		new_labels.append("%s %s"%(chr,int(start+i+1)))
		flag=False
	else:
		new_labels.append(int(start+i+1))
plt.xticks(locs,new_labels)
plt.title("ID: %s, CHR: %s, POS: %s, REF: %s, ALT: %s \n Predicted efficiency: %.1f"%(variant_id,chr,variant_pos,ref,alt,predicted_efficiency)+"%")
ax.figure.savefig('test.pdf', bbox_inches='tight')












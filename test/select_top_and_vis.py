
# exec(open("vis_pegRNA.py").read())
from dna_features_viewer import GraphicFeature, GraphicRecord
import pandas as pd
import uuid
import os
import subprocess
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import os
os.system("rm *.pdf")
def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()


def get_gRNA_cut_site(start,end,strand,offset=-3):
	# return 1-index pos
	# cut site = the first base before cut near PAM direction

	if strand == "+":
		return int(end + offset)
	if strand == "-":
		return int(start - offset +1)



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
	

def plot(df,pegRNA_id):

	variant_id = pegRNA_id.split("_chr")[0]
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
	ax.figure.savefig('%s.top_pegRNA.pdf'%(pegRNA_id), bbox_inches='tight')
	return df

def get_strand(x):
	if x == "+":
		return 1
	else:
		return -1


npg_colors = ["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#464d4f"]
my_colors = {}
my_colors['sgRNA'] = npg_colors[0]
my_colors['PBS'] = npg_colors[1]
my_colors['RTS'] = npg_colors[2]
my_colors['nick-gRNA'] = npg_colors[3]
my_colors['variant'] = "#e6fc3f"

import sys
rawX = "sankaran_2016_cell.top_pegRNAs.csv"
# rawX = sys.argv[1]

df = pd.read_csv(rawX,index_col=0)
df2 = df.copy()
top_list = []
df2['cut1'] = df2.apply(lambda r:get_gRNA_cut_site(r.start,r.end,r.strand,offset=-3),axis=1)
df2['distance'] = abs(df2['POS']-df2['cut1'])
df3 = df2[df2['type']=='sgRNA']
filter1 = df3[df3['distance']<=20].index.tolist()
df4 = df2[df2['type']=='RTS']
df4['distance'] = df4.apply(lambda r:min(abs(r.start+1-r.POS),abs(r.end-r.POS)),axis=1)
filter2 =  df4[df4['distance']>=5].index.tolist()

top_list = list(set(filter2).intersection(filter1))
df2 = df2.loc[top_list]
df2['group'] = [x.split("_chr")[0] for x in df2.index.tolist()]
# df3['group'] = [x.split("_candidate")[0] for x in df.index.tolist()]
# top_list = []
# for s,d in df3.groupby('group'):
	# d = d.sort_values()
df2 = df2.sort_values("predicted_efficiency",ascending=False)
df2 = df2.drop_duplicates('group')
# exec(open("select_top_and_vis.py").read())
df_list =  Parallel(n_jobs=-1,verbose=10)(delayed(plot)(df,x) for x in df2.index.tolist())

out = pd.concat(df_list)
out.to_csv("%s.top1.pegRNA.csv"%(rawX))

# os.system("")
os.system("mkdir -p vis_files")
os.system("mv *.pdf vis_files")




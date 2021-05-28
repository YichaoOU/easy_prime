# from dna_features_viewer import GraphicFeature, GraphicRecord
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import dash_table
from flask import Flask, send_from_directory
from urllib.parse import quote as urlquote
import dash_bootstrap_components as dbc

# import plotly_express as px
from dash.dependencies import Input, Output, State
import re
import os
import numpy as np
from sklearn.cluster import KMeans
import plotly.figure_factory as ff
import os
import base64
import datetime
from datetime import datetime as dt
import pathlib
import easy_prime
import io
import csv
from easy_prime.utils import get_parameters, print_parameters,vcf2fasta,fasta2vcf
from easy_prime import target_mutation
import subprocess
from joblib import Parallel, delayed
import urllib
from io import StringIO
import uuid
from Bio import SeqIO

import pandas as pd
import uuid
import subprocess
import matplotlib.pyplot as plt
npg_colors = ["#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F","#464d4f"]
my_colors = {}
my_colors['sgRNA'] = npg_colors[0]
my_colors['PBS'] = npg_colors[1]
my_colors['RTT'] = npg_colors[2]
my_colors['ngRNA'] = npg_colors[3]
my_colors['variant'] = "#e6fc3f"


def get_version():
	
	return "v%s"%(easy_prime.__version__)
	
def vis_pegRNA_png(df_file,jid):
	# bedtools 2.29.2
	# print ("call easy_prime vis")
	genome_fasta = "hg19.fa"
	if os.path.isfile(f"results/{jid}.fa"):
		subprocess.call(f"easy_prime_vis -f {df_file} -s results/{jid}.fa --output_file_name results/{jid}.png",shell=True)
	else:
		# print ("using genome fasta")
		cmd = f"easy_prime_vis -f {df_file} -s {genome_fasta} --output_file_name results/{jid}.png"
		# print (cmd)
		subprocess.call(cmd,shell=True)
	fig = f"results/{jid}.png"
	with open(fig, "rb") as image_file:
		img_string = base64.b64encode(image_file.read())
	return "data:image/png;base64,%s"%(img_string.decode("utf-8"))


def df2csv_string(df):
	csv_string = df.to_csv(index=False, encoding='utf-8')
	csv_string = "data:text/csv;charset=utf-8,%EF%BB%BF" + urllib.parse.quote(csv_string)
	return csv_string

def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "results/{}".format(urlquote(filename))
    return html.A(filename, href=location,target="_blank")

def get_current_pegRNA_table_title_and_download_links(df,jid):
	PE = "%.1f"%(df.predicted_efficiency[0])+"%"
	use_columns = ['chr','start','end','seq',"type",'strand']
	show_columns = ['#chr','start','end','seq',"type",'strand']
	df = df[use_columns]
	df.columns = show_columns
	table =  dash_table.DataTable(
		id='pegRNA-table',
		columns=[
			{'name': i, 'id': i, 'deletable': False} for i in show_columns
		],
		data=df.to_dict('records'),
	)
	
	header = dbc.FormGroup(
		[
			html.H5("Current pegRNA/ngRNA selection. Predicted efficiency: %s"%(PE),style={"margin-right":10,"margin-left":10}),
			html.A(html.Button('Download current selection',className="btn btn-dark"),href=df2csv_string(df[show_columns]),target="_blank",download="current_design.csv" ,style={"margin-right":10}),
			html.A(html.Button('Download all predictions',className="btn btn-dark"),href="results/{}".format(urlquote("%s_rawX_pegRNAs.csv.gz"%(jid))),target="_blank",style={"margin-right":10}),
		],
		row=True,
	)
	return [header,table]

def get_current_selection_table(df):
	show_columns = ['chr','start','end','seq',"predicted_efficiency",'strand']
	table =  dash_table.DataTable(
		id='pegRNA-table',
		columns=[
			{'name': i, 'id': i, 'deletable': False} for i in show_columns
		],
		data=df.to_dict('records'),
	)

	return [html.H5("Current pegRNA/ngRNA selection"),table]


download_current_selection_button = html.A(html.Button('Submit feedback!'),href='https://github.com/czbiohub/singlecell-dash/issues/new',target="_blank")


#--------------------------------- easy prime search ---------------------------
def write_fasta(file_name,myDict):
	out = open(file_name,"wt")
	for k in myDict:
		out.write(">"+k+"\n")
		out.write(myDict[k]+"\n")
	out.close()
def read_fasta(f):
	my_dict = {}
	for r in SeqIO.parse(f, "fasta"):
		my_dict[r.id] = str(r.seq).upper()
	return my_dict	
def run_steps(t,**kwargs):

	t.init(**kwargs)
	t.search(**kwargs)
	t.predict(**kwargs)

	return [t.topX,t.rawX,t.X_p,t.found_PE3b,t.found_PE3,t.found_dPAM,t.found_PE2,t.N_sgRNA_found]

def run_easy_prime_backend():

	error_message = ""
	return False,error_message

def run_easy_prime_backend(vcf,jid,parameters):
	# print (vcf.head())
	if vcf.shape[1]==5:
		vcf[5] = vcf2fasta(vcf,**parameters)
		vcf = vcf[list(range(6))]
	
	variant_list = vcf[2].tolist()
	my_targets = [target_mutation(*r) for i,r in vcf.iterrows()]	

	df_list = [run_steps(t,**parameters) for t in my_targets]

	summary = pd.DataFrame([x[3:8] for x in df_list]).astype(int)
	summary.columns = ['found_PE3b','found_PE3','found_dPAM','found_PE2',"N_sgRNA_found"]
	summary.index = variant_list
	summary.to_csv("results/%s_summary.csv"%(jid),index=True)

	df_top = pd.concat([x[0] for x in df_list])
	if df_top.shape[0]==0:
		# print ("no pegRNA (including PE2) were found for the input file")
		return summary,pd.DataFrame(),pd.DataFrame(),pd.DataFrame()
	df_top = df_top.sort_values("predicted_efficiency",ascending=False)
	df_top.to_csv("results/%s_topX_pegRNAs.csv"%(jid),index=False)
	df_all = pd.concat([x[1] for x in df_list])
	df_all = df_all.sort_values("predicted_efficiency",ascending=False)
	df_all.to_csv("results/%s_rawX_pegRNAs.csv.gz"%(jid),index=False,compression="gzip")
	
	X_p = pd.concat([x[2] for x in df_list])
	X_p = X_p.sort_values("predicted_efficiency",ascending=False)
	X_p.to_csv("results/%s_X_p_pegRNAs.csv.gz"%(jid),index=True,compression="gzip")
	return summary,df_top,df_all,X_p

def PD2fasta_dict(input_string):
	out = {}
	myDict = read_fasta(StringIO(input_string))
	for k in myDict:
		s = myDict[k]
		if s.count("(")>1:
			return 0
		before = s[:s.find("(")]
		after = s[s.find(")")+1:]
		mutation = s[s.find("(")+1:s.find(")")]
		if "+" in mutation:
			mutation = mutation.replace("+","")
			ref = before+after
			alt = before+mutation+after
		elif "-" in mutation:
			mutation = mutation.replace("-","")
			ref = before+mutation+after
			alt = before+after
		else:
			mutation = mutation.split("/")
			ref = before+mutation[0]+after
			alt = before+mutation[1]+after
		out["%s_ref"%(k)] = ref
		out["%s_alt"%(k)] = alt
	return out
main_chr_human = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chrX","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr20","chrY","chr19","chr22","chr21","chrM"]

def is_valid_DNA(sequence):
	valid_dna = "ACGT"
	return all(i in valid_dna for i in sequence)

def check_and_convert_input(input_type,chr,pos,variant_id,ref,alt,vcf_batch,fasta_input,PrimeDesign_input,jid):
	error_flag = False
	error_message = ""
	if input_type == "vcf_tab":
		if not "chr" in chr:
			chr = "chr"+chr
		if not chr in main_chr_human:
			return None, True,"Input chromosome %s not found in hg19 main chromosomes"%(chr)
		try:
			pos = int(pos)
		except:
			return None, True,"Input Position %s can't converted to integer"%(pos)
		name = variant_id
		if name == "":
			return None, True,"Input Name is empty!"
		ref = ref.upper()
		if ref == "":
			return None, True,"Input Name is empty!"
		alt = alt.upper()
		if alt == "":
			return None, True,"Input Name is empty!"
		if not is_valid_DNA(ref):
			return None, True,"Input reference allele is not valid, only A,C,G,T is allowed!"
		if not is_valid_DNA(alt):
			return None, True,"Input alternative allele is not valid, only A,C,G,T is allowed!"
		vcf = pd.DataFrame([[chr,pos,name,ref,alt]])

	if input_type == "vcf_batch_tab":
		try:
			vcf = pd.read_csv(StringIO(vcf_batch),comment="#",sep="\t",header=None)
			vcf[1] = vcf[1].astype(int)
			vcf =vcf.drop_duplicates(2) # remove duplicated names
			vcf[3] = [x.upper() for x in vcf[3]]
			vcf[4] = [x.upper() for x in vcf[4]]
		except Exception as e:
			# print (e)
			return None, True,"Input vcf_batch can't be parsed correctly! %s"%(e)
	if input_type == "fasta_tab":
		try:
			file_name = "results/%s.fa"%(jid)
			myDict = read_fasta(StringIO(fasta_input))
			write_fasta(file_name,myDict)
			vcf = fasta2vcf(file_name)
		except Exception as e:
			# print (e)
			return None,True,"Input fasta can't be parsed correctly! Make sure fasta length >= 50bp. Error: %s"%(e)
	if input_type == "PrimeDesign_tab":
		try:
			file_name = "results/%s.fa"%(jid)
			myDict = PD2fasta_dict(PrimeDesign_input)
			write_fasta(file_name,myDict)
			vcf = fasta2vcf(file_name)
		except Exception as e:
			# print (e)
			return None,True,"Input PrimeDesign sequences can't be parsed correctly! Currently, combinatorial editing in this format is not supported! %s"%(e)


	return vcf,False,error_message



def read_easy_prime_output(jid):
	rawX = pd.read_csv("results/%s_rawX_pegRNAs.csv.gz"%(jid))
	X_p = pd.read_csv("results/%s_X_p_pegRNAs.csv.gz"%(jid),index_col=0)
	X_p['sample_ID'] = X_p.index.tolist()
	rawX.index = rawX.sample_ID.tolist()
	
	rawX['location_name'] = rawX.chr+"_"+rawX.start.astype(str)+"_"+rawX.end.astype(str)+"_"+rawX.seq
	rawX['DeepSpCas9_score'] = X_p['cas9_score']
	rawX['target_pos'] = X_p['Target_pos']
	rawX['nick_pos'] = X_p['nick_to_pegRNA']
	rawX['PBS_length'] = X_p['PBS_length']
	rawX['RTT_length'] = X_p['RTT_length']
	return rawX,X_p

def get_options_dict(jid):
	df = pd.read_csv("results/%s_summary.csv"%(jid),index_col=0)
	# print (df.head())
	out = []
	first_valid = None
	for i,r in df.iterrows():
		if r.found_PE2 == 0:
			out.append({'label': i + " | no pegRNA found",'disabled': True, 'value': i})
		else:
			if not first_valid:
				first_valid = i
			out.append({'label': i, 'value': i})
	return out,first_valid

def get_annotation(x):
	if "dPAM" in x:
		return "PAM-disruption"
	if "PE3b" in x:
		return "PE3b"
	return ""
	
def to_sgRNA_table(rawX,sample_ID):
	rawX_df = rawX[rawX.sample_ID.str.contains(sample_ID)]
	rawX_df = rawX_df[rawX_df['type']=='sgRNA']
	# X_p_df = X_p[X_p.sample_ID.str.contains(sample_ID)]
	rawX_df = rawX_df.sort_values('predicted_efficiency',ascending=False)
	rawX_df = rawX_df.drop_duplicates(['chr','start','end'])
	columns = ['chr','start','end','seq','DeepSpCas9_score','strand','target_pos','annotation']
	
	rawX_df['annotation'] = rawX_df.sample_ID.apply(get_annotation)
	rawX_df = rawX_df.reset_index(drop=True)
	return rawX_df[columns]

def to_sgRNA_table2(rawX,X_p,sample_ID):
	rawX_df = rawX[rawX.sample_ID.str.contains(sample_ID)]
	rawX_df = rawX_df[rawX_df['type']=='sgRNA']
	X_p_df = X_p[X_p.sample_ID.str.contains(sample_ID)]
	rawX_df = rawX_df.sort_values('predicted_efficiency',ascending=False)
	rawX_df = rawX_df.drop_duplicates(['chr','start','end'])

	columns = ['chr','start','end','seq','DeepSpCas9_score','strand','target_pos','annotation']
	rawX_df['DeepSpCas9_score'] = X_p_df['cas9_score']
	rawX_df['target_pos'] = X_p_df['Target_pos']
	rawX_df['annotation'] = ""
	rawX_df = rawX_df.reset_index(drop=True)
	ID_list = rawX_df.sample_ID.tolist()
	return rawX_df[columns],ID_list


def to_ngRNA_table(rawX,sgRNA_location):
	sample_ID_list = rawX[rawX.location_name==sgRNA_location].sample_ID.tolist()
	rawX_df = rawX[rawX.sample_ID.isin(sample_ID_list)]
	rawX_df = rawX_df[rawX_df['type']=='ngRNA']
	rawX_df = rawX_df.drop_duplicates(['chr','start','end'])
	columns = ['chr','start','end','seq','nick_pos','strand']
	rawX_df = rawX_df.sort_values('nick_pos')
	rawX_df = rawX_df.reset_index(drop=True)
	return rawX_df[columns],rawX_df.predicted_efficiency.idxmax()


def to_ngRNA_table2(rawX,X_p,sample_ID_list,best_ID=None):
	if not best_ID:
		best_ID = sample_ID_list[0]
	rawX_df = rawX[rawX.sample_ID.isin(sample_ID_list)]
	X_p_df = X_p[X_p.sample_ID.isin(sample_ID_list)]
	rawX_df = rawX_df[rawX_df['type']=='ngRNA']
	rawX_df = rawX_df.drop_duplicates(['chr','start','end'])
	columns = ['chr','start','end','seq','nick_pos','strand']
	rawX_df['nick_pos'] = X_p_df['nick_to_pegRNA']
	rawX_df = rawX_df.reset_index(drop=True)
	return rawX_df[columns],rawX_df.predicted_efficiency.idxmax()



def to_PBS_table(rawX,sgRNA_location):
	sample_ID_list = rawX[rawX.location_name==sgRNA_location].sample_ID.tolist()
	rawX_df = rawX[rawX.sample_ID.isin(sample_ID_list)]
	rawX_df = rawX_df[rawX_df['type']=='PBS']
	columns = ['chr','start','end','seq','PBS_length','strand']
	rawX_df = rawX_df.drop_duplicates('seq')
	# rawX_df['PBS_length'] = rawX_df.seq.apply(len)
	rawX_df = rawX_df.sort_values('PBS_length')
	rawX_df = rawX_df.reset_index(drop=True)
	return rawX_df[columns],rawX_df.predicted_efficiency.idxmax()


def to_PBS_table2(rawX,sample_ID_list,best_ID=None):
	if not best_ID:
		best_ID = sample_ID_list[0]
	rawX_df = rawX[rawX.sample_ID.isin(sample_ID_list)]
	rawX_df = rawX_df[rawX_df['type']=='PBS']
	columns = ['chr','start','end','seq','PBS_length','strand']
	rawX_df = rawX_df.drop_duplicates('seq')
	rawX_df['PBS_length'] = rawX_df.seq.apply(len)
	rawX_df = rawX_df.sort_values('PBS_length')
	rawX_df = rawX_df.reset_index(drop=True)
	return rawX_df[columns],rawX_df[rawX_df.sample_ID == best_ID].index.tolist()

def to_RTT_table(rawX,sgRNA_location):
	sample_ID_list = rawX[rawX.location_name==sgRNA_location].sample_ID.tolist()
	rawX_df = rawX[rawX.sample_ID.isin(sample_ID_list)]
	rawX_df = rawX_df[rawX_df['type']=='RTT']
	columns = ['chr','start','end','seq','RTT_length','strand']
	rawX_df = rawX_df.drop_duplicates('seq')
	# rawX_df['RTT_length'] = rawX_df.seq.apply(len)
	rawX_df = rawX_df.sort_values('RTT_length')
	rawX_df = rawX_df.reset_index(drop=True)
	return rawX_df[columns],rawX_df.predicted_efficiency.idxmax()


def to_RTT_table2(rawX,sample_ID_list,best_ID=None):
	if not best_ID:
		best_ID = sample_ID_list[0]
	rawX_df = rawX[rawX.sample_ID.isin(sample_ID_list)]
	rawX_df = rawX_df[rawX_df['type']=='RTT']
	columns = ['chr','start','end','seq','RTT_length','strand']
	rawX_df = rawX_df.drop_duplicates('seq')
	rawX_df['RTT_length'] = rawX_df.seq.apply(len)
	rawX_df = rawX_df.sort_values('RTT_length')
	rawX_df = rawX_df.reset_index(drop=True)
	return rawX_df[columns],rawX_df[rawX_df.sample_ID == best_ID].index.tolist()


def get_uid():
	# return "easy_prime_yli11_2021-05-24_result_dir"
	return str(uuid.uuid4()).split("-")[-1]




#--------------------------------- dash app utils ---------------------------


def df2bedjs(df,output):
	# sample_ID,CHROM,POS,REF,ALT,type,seq,chr,start,end,strand,predicted_efficiency
	# FIG5G_HEK293T_HEK3_6XHIS_chr9_110184619_110184639_+_GGCCCAGACTGAGCACGTGA_candidate_1808,chr9,110184636,G,GCACCATCATCACCATCAT,ngRNA,GTCAACCAGTATCCCGGTGC,chr9,110184723,110184743,-,0.6275171041488647
	df['name'] = df.apply(lambda r:"""{"strand":"%s","name":"%s","color":"%s"}"""%(r.strand,r['type'],my_colors[r['type']]),axis=1)
	track_name = "pegRNA_design_%.1f"%(df.predicted_efficiency[0])
	df = df[['chr','start','end','name']]
	# print (df.head())
	df.sort_values('start').to_csv("results/%s.bed"%(output),sep="\t",header=False,index=False,quoting=csv.QUOTE_NONE)
	os.system("bgzip results/{0}.bed;tabix -p bed results/{0}.bed.gz".format(output))
	
	return '''{"type":"bedj","url":"http://easy-prime-test-dev.us-west-2.elasticbeanstalk.com/results/%s.bed.gz","stackheight":20,"stackspace":1,"name":"%s"},'''%(output,track_name)


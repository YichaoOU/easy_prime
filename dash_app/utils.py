from dna_features_viewer import GraphicFeature, GraphicRecord
import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import dash_table
from flask import Flask, send_from_directory
from urllib.parse import quote as urlquote

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
# TESTDATA = StringIO("""col1;col2;col3
#     1;4.4;99
#     2;4.5;200
#     3;4.7;65
#     4;3.2;140
#     """)

# df = pd.read_csv(TESTDATA, sep=";")

def get_version():
	
	return "v%s"%(easy_prime.__version__)

#---------- easy prime vis code ---------------------

def write_file(file_name,message):
	out = open(file_name,"wt")
	out.write(message)
	out.close()
def get_fasta_single(chr,start,end,genome_fasta=None):
	out_bed = str(uuid.uuid4()).split("-")[-1]
	out_fa = str(uuid.uuid4()).split("-")[-1]
	write_file(out_bed,"%s\t%s\t%s"%(chr,start,end))
	command = "bedtools getfasta -fi %s -bed %s -fo %s -tab"%(genome_fasta,out_bed,out_fa)
	subprocess.call(command,shell=True)
	lines = open(out_fa).readlines()[0]
	seq = lines.split()[-1]
	subprocess.call("rm %s;rm %s"%(out_bed,out_fa),shell=True)
	return seq
def get_strand(x):
	if x == "+":
		return 1
	else:
		return -1

def vis_pegRNA(df_file,genome_fasta=None,out_file_name=None,**kwargs):

	# print ("call easy_prime vis")
	subprocess.call(f"easy_prime_vis -f {df_file} -s {genome_fasta} --output_file_name results/{out_file_name}.png",shell=True)
	fig = f"results/{out_file_name}.png"
	with open(fig, "rb") as image_file:
		img_string = base64.b64encode(image_file.read())
	return "data:image/png;base64,%s"%(img_string.decode("utf-8"))


def vis_pegRNA2(df,genome_fasta=None,**kwargs):
	"""Given one instance of easy-prime prediction (rawX format), generate DNA visualization
	
	Input
	--------
	the data frame contains 4 rows: RTT, PBS, sgRNA, ngRNA
	
	"""
	pegRNA_id = df.index.tolist()[0]
	variant_id = pegRNA_id.split("_")[0]
	chr = df['CHROM'][0]
	start = df['start'].min()
	start -= start%10
	start -= 1
	end = df['end'].max()
	end -= end%10
	end += 10

	variant_pos = df.POS.min()
	ref = df.REF[0]
	alt = df.ALT[0]
	predicted_efficiency = df.predicted_efficiency[0]*100
	pos = variant_pos-start
	sequence = get_fasta_single(chr,start,end,genome_fasta).upper()
	fig,ax = plt.subplots()
	feature_list = []
	for s,r in df.iterrows():
		r_start = r.start-start
		r_end = r_start+(r.end-r.start)
		r_strand = get_strand(r.strand)
		gf = GraphicFeature(start=r_start, end=r_end, strand=r_strand, 
			color=my_colors[r.type],label=r.type)
		feature_list.append(gf)
	record = GraphicRecord(sequence=sequence, features=feature_list)

	# ax, _ = record.plot(figure_width=int(len(sequence)/5))
	record.plot(ax=ax,figure_width=int(len(sequence)/5))
	return 0
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

	my_stringIObytes = io.BytesIO()
	ax.figure.savefig(my_stringIObytes, format='png',bbox_inches='tight')
	my_stringIObytes.seek(0)
	img_string = base64.b64encode(my_stringIObytes.read())
	return "data:image/png;base64,%s"%(img_string.decode("utf-8"))


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


def easy_prime_main(input_data,jid,parameters):
	
	## read vcf
	# print (f"jid: {jid}")
	# vcf = pd.read_csv(input_data,comment="#",sep="\t",header=None)
	# vcf[1] = vcf[1].astype(int)
	# vcf =vcf.drop_duplicates(2) # remove duplicated names
	# vcf[3] = [x.upper() for x in vcf[3]]
	# vcf[4] = [x.upper() for x in vcf[4]]
	# vcf[5] = vcf2fasta(vcf,**parameters)
	# vcf = vcf[list(range(6))]


	## get a list of targets
	# from easy_prime import target_mutation
	# import pandas as pd
	## ## modified for fasta input
	try:
		vcf = pd.read_csv(StringIO(input_data),comment="#",sep="\t",header=None)
		vcf[1] = vcf[1].astype(int)
		vcf =vcf.drop_duplicates(2) # remove duplicated names
		vcf[3] = [x.upper() for x in vcf[3]]
		vcf[4] = [x.upper() for x in vcf[4]]
		vcf[5] = vcf2fasta(vcf,**parameters)
		vcf = vcf[list(range(6))]
		
		## for each target, create target mutation class
		
	except:
		try:
			# print (input_data)
			# print ("Reading fasta file: %s"%(input_data))
			
			# print (list(SeqIO.parse(StringIO(input_data),"fasta")))
			file_name = "results/%s.fa"%(jid)
			myDict = read_fasta(StringIO(input_data))
			write_fasta(file_name,myDict)
			vcf = fasta2vcf(file_name)
			# print (vcf)
		except Exception as e:
			print (e)
			# print ("Can't read %s as vcf or fasta. Please check input. Exit..."%(input_data))
			exit()

	variant_list = vcf[2].tolist()
	my_targets = [target_mutation(*r) for i,r in vcf.iterrows()]	

	## for each target, create target mutation class
	# my_targets = [target_mutation(*r) for i,r in vcf.iterrows()]
	# print (vcf)
	## find best pegRNAs
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

#--------------------------------- Right Column --------------------------------------

def init_fig():
	# print ("init figure")
	fig = "img/init.png"
	with open(fig, "rb") as image_file:
		img_string = base64.b64encode(image_file.read())
	return "data:image/png;base64,%s"%(img_string.decode("utf-8"))

def encode_fig(pegRNA_id=None):
	# print ("#######????????????#######",pegRNA_id)
	fig = "img/rs737092_chr20_55990370_55990390_TACCTCCTGGGCCTCCCGCA_candidate_162.top_pegRNA.png"
	with open(fig, "rb") as image_file:
		img_string = base64.b64encode(image_file.read())
	return "data:image/png;base64,%s"%(img_string.decode("utf-8"))

def vis_PE_design(app,pegRNA_list):
	## TODO add option when pegRNA_list=None
	# print (list(pegRNA_list[0].values())[0])
	dropdown = dcc.Dropdown(
		id="select_pegRNA_id",
		options=pegRNA_list,
		# value=list(pegRNA_list[0])[0],
		value="rs737092",
		style={"display": "inline-block","margin-left":"20px","width":"200px",'font-size':'13px'},
		clearable=False
	)
	## TODO complete the callback with real pegRNA prediction
	@app.callback(Output('PE_design_figure', 'src'),
		[Input('select_pegRNA_id', 'value')])
	def update_output(pegRNA_id):
		return encode_fig(pegRNA_id)


	return html.Div(
		id="vis_top_pegRNA",
		style={"background-color":"#FFFFFF","margin-top":"20px",'font-weight':'bold','font-size':'20px'},
		children=[
			html.Div(style={"vertical-align": "middle","display": "inline-block"},children=[
				# html.Div(style={"padding-bottom": "40px","display": "inline-block"},children=['PE design visualization for:'])
				'PE design visualization for:',
				dropdown
			]),
			html.Div([html.Img(id='PE_design_figure')]),		

		],
	)




def show_PE_table(app):
	## TODO complete this
	df = pd.read_csv("data/sankaran_2016_cell.top_pegRNAs.csv.top1.pegRNA.csv")
	df['ID'] = [x.split("_chr")[0] for x in df['sample_ID'].tolist()]
	show_columns = ['ID','type','seq','chr','start','end','strand']

	@app.callback(Output("pegRNA-table", "data"),
		[Input('select_pegRNA_id', 'value')])
	def update_output(pegRNA_id):
		# print ("update_output")
		return df.to_dict('records')
	@app.callback(Output('topX', 'href'),
		[Input('pegRNA-table', "page_current"),
		Input('pegRNA-table', "page_size"),
		Input('pegRNA-table', "sort_by"),
		Input('pegRNA-table', "filter_query"),
		])
	def download_topX(page_current, page_size, sort_by, filter):
		csv_string = df.to_csv(index=False, encoding='utf-8')
		# print (df.head())
		# print ("#### printing topX ####")
		# csv_string = "data:text/csv;charset=utf-8,%EF%BB%BF" + urllib.parse.quote(csv_string)
		csv_string = "data:text/csv;charset=utf-8,%EF%BB%BF" + urllib.parse.quote(csv_string)
		# print (csv_string)
		return csv_string

	# print (df.head())
	return dash_table.DataTable(
			id='pegRNA-table',
			columns=[
				{'name': i, 'id': i, 'deletable': False} for i in show_columns
			],
			style_cell={
				'minWidth': '0px', 'maxWidth': '20px',
				'width':'10%',
				'whiteSpace': 'normal'
			},			
			editable=True,
			filter_action="native",
			page_action="native",
			page_current= 0,
			page_size= 13,
			# data=df.to_dict('records'),
		)


# @app.callback(Output('topX', 'href')
#             , [Input('full-table', "page_current"),
#      Input('full-table', "page_size"),
#      Input('full-table', 'sort_by'),
#      Input('full-table', 'filter_query')])
# def update_table2(page_current, page_size, sort_by, filter):
#     filtering_expressions = filter.split(' && ')
#     dff = df
#     for filter_part in filtering_expressions:
#         col_name, operator, filter_value = functions.split_filter_part(filter_part)

#         if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
#             # these operators match pandas series operator method names
#             dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
#         elif operator == 'contains':
#             dff = dff.loc[dff[col_name].str.contains(filter_value)]
#         elif operator == 'datestartswith':
#             # this is a simplification of the front-end filtering logic,
#             # only works with complete fields in standard format
#             dff = dff.loc[dff[col_name].str.startswith(filter_value)]

#     if len(sort_by):
#         dff = dff.sort_values(
#             [col['column_id'] for col in sort_by],
#             ascending=[
#                 col['direction'] == 'asc'
#                 for col in sort_by
#             ],
#             inplace=False
#         )


#     csv_string = dff.to_csv(index=False, encoding='utf-8')
#     csv_string = "data:text/csv;charset=utf-8,%EF%BB%BF" + urllib.parse.quote(csv_string)
#     return csv_string





















#--------------------------------- Left Column --------------------------------------




def buttons():




	return html.Div([html.Button('Start', id='start_search',style={'margin-left':'8%'}),html.A(html.Button('Tutorial', id='help_menu',style={"float":"right",'margin-right':'8%'}),href="https://github.com/YichaoOU/easy_prime",target='_blank')])


def help():

	return html.Div([html.Button('Help', id='help_menu')])

intro="Easy-Prime is a machine learning based tool for prime editing gRNA design. Please input your desire edits in VCF format or as FASTA sequence and click start.  Additionally, you can play with the pegRNA searching parameters or choose a different genome in the bottom."

def welcome():
	return html.Div(
		style={"margin-bottom":"10%"},
		id="welcome",
		children=[
			html.H5("Welcome to Easy-Prime!"),intro,
		],
	)

def variant_input():
	variants = dcc.Textarea(
		id='variants',
		value='VCF format\nThe first 5 columns: chr, id, pos, ref, alt\n #comment line will be ignored\n FASTA format\nUse _ref to represent reference allele and _alt for alternative allele.\n>test_ref\nXXXXXXX\n>test_alt\nXXXXAXXXX',
		style={'width': '100%', 'height': 150,'padding':"0 0 0 0","font-size":"13px"},
	)
	return html.Div(
		id="variant_input",
		style={'font-weight':'bold',"font-size":"20px","margin-left":"3%","margin-right":"3%"},
		children=[			
			"Input Variants",
			variants,
		],
	)

def genome_input(app):
	dropdown = dcc.Dropdown(
		id="genome",
		options=[
			{'label': 'hg19', 'value': 'hg19'},
			{'label': 'hg38', 'value': 'hg38'},
			{'label': 'mm10', 'value': 'mm10'},
			{'label': 'mm9', 'value': 'mm9'},
			{'label': 'custom', 'value': 'custom'},

		],
		value="hg19",
		clearable=False
	)
	custom_geome_textarea = dcc.Textarea(
		id='user_genome',
		value='Please input fasta format sequences',
		style={'width': '100%', 'height': 300},
	)


	@app.callback(dash.dependencies.Output('triger_user_input', 'children'),
		[dash.dependencies.Input('genome', 'value')])
	def user_input_genome(value):
		if value == "custom":
			return custom_geome_textarea
		else:
			return None
	return html.Div(
		style={'font-weight':'bold',"margin-left":"8%","margin-right":"8%"},
		children=[			
			html.Div(
				style={'font-weight':'bold'},
				children=["Choose a genome (also support custome sequence)",dropdown],
			),
			html.Div(id='triger_user_input'),
		],
	)

def PBS_slider(app):
	labelDict={}
	labelDict[7] = {'label':'min 7','style': {'color': '#77b0b1'}}
	labelDict[12] = {'label':'12','style': {'color': '#77b0b1'}}
	labelDict[17] = {'label':'17','style': {'color': '#77b0b1'}}
	subtitle = 'Please select the maximum PBS length'
	subtitle = ''
	title = 'PBS length range'
	# return mySlider("pbs",7,17,labelDict,defaultValue=12,title=title,subtitle=subtitle,app=app)
	return myRangeSlider("pbs",7,17,labelDict,defaultValue=[10,14],title=title,subtitle=subtitle,app=app)

def ngRNA_slider(app):
	labelDict={}
	labelDict[0] = {'label':'min 0','style': {'color': '#77b0b1'}}
	labelDict[50] = {'label':'50','style': {'color': '#77b0b1'}}
	labelDict[100] = {'label':'100','style': {'color': '#77b0b1'}}
	labelDict[150] = {'label':'150','style': {'color': '#FF5500'}}
	labelDict[200] = {'label':'200','style': {'color': '#FF0000'}}
	subtitle = 'Please select the maximum ngRNA distance to pegRNA'
	subtitle = ''
	title = 'ngRNA to pegRNA distance range'
	return mySlider("ngRNA",0,200,labelDict,defaultValue=100,title=title,subtitle=subtitle,app=app)


def RTT_slider(app):
	labelDict={}
	labelDict[7] = {'label':'min 7','style': {'color': '#77b0b1'}}
	labelDict[20] = {'label':'20','style': {'color': '#77b0b1'}}
	labelDict[30] = {'label':'30','style': {'color': '#77b0b1'}}
	# labelDict[30] = {'label':'30','style': {'color': '#77b0b1'}}
	# labelDict[40] = {'label':'40','style': {'color': '#77b0b1'}}
	labelDict[50] = {'label':'50','style': {'color': '#FF5500'}}
	# labelDict[60] = {'label':'60','style': {'color': '#FF5500'}}
	# labelDict[70] = {'label':'70','style': {'color': '#FF0000'}}
	labelDict[80] = {'label':'80','style': {'color': '#FF0000'}}
	subtitle = 'Please select the maximum RTT length'
	subtitle = ''
	title = 'RTT length range'
	# return mySlider("rtt",7,80,labelDict,defaultValue=40,title=title,subtitle=subtitle,app=app)
	return myRangeSlider("rtt",7,80,labelDict,defaultValue=[10,20],title=title,subtitle=subtitle,app=app)

def mySlider(name,vmin,vmax,labelDict,defaultValue,title=None,subtitle=None,app=None):
	slider = dcc.Slider(
		id=f"slider-{name}",
		className="left-align-slider",
		min=vmin,
		max=vmax,
		value=defaultValue,
		marks=labelDict,
		included=True
	)
	@app.callback(dash.dependencies.Output(f'slider-output-{name}', 'children'),
		[dash.dependencies.Input(f"slider-{name}", 'value')])
	def update_output(value):
		return f"{title}: [{vmin}, {value}]"


	return html.Div(
		style={'font-weight':'bold',"margin-left":"8%","margin-right":"8%","margin-top":"8%","margin-bottom":"8%"},
		children=[
			html.Div(id=f'slider-output-{name}'),
			html.Div(
				style={"color": "grey"},
				children=[subtitle,slider,html.Div(id=f'{name}-output-container')],
			),
		],
	)

def myRangeSlider(name,vmin,vmax,labelDict,defaultValue,title=None,subtitle=None,app=None):
	slider = dcc.RangeSlider(
		id=f"slider-{name}",
		className="left-align-slider",
		min=vmin,
		max=vmax,
		value=defaultValue,
		marks=labelDict,
		included=True
	)
	@app.callback(dash.dependencies.Output(f'slider-output-{name}', 'children'),
		[dash.dependencies.Input(f"slider-{name}", 'value')])
	def update_output(value):
		return f"{title}: [{value[0]}, {value[1]}]"


	return html.Div(
		style={'font-weight':'bold',"margin-left":"8%","margin-right":"8%","margin-top":"8%","margin-bottom":"8%"},
		children=[
			html.Div(id=f'slider-output-{name}'),
			html.Div(
				style={"color": "grey"},
				children=[subtitle,slider,html.Div(id=f'{name}-output-container')],
			),
		],
	)

#--------------------------------- download data --------------------------------------



import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import dash_table
from urllib.parse import quote as urlquote
import dash_bootstrap_components as dbc


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
import io
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


#--------------------------------- welcome panel --------------------------------------


intro="Easy-Prime is a machine learning based tool for prime editing gRNA (pegRNA) design. Please input your desire edits in VCF format or as FASTA sequence and click start. Additionally, you can play with the pegRNA/ngRNA searching parameters below."

def welcome():
	return html.Div(
		style={"margin-bottom":10},
		id="welcome",
		children=[
			html.H5("Welcome to Easy-Prime!",style={"color":"#2c8cff"}),intro,
		],
	)

#--------------------------------- Output utils --------------------------------------

def display_table(df,table_name):

	return dash_table.DataTable(
					id = f'{table_name}-table',
					columns = [{'name': i, 'id': i} for i in df.columns.tolist()],
					data = df.to_dict('records'),
					style_cell={'textAlign': 'left', 'padding': '5px'},
					# style_as_list_view=True,
					style_header={
						'backgroundColor': 'white',
						# 'fontWeight': 'bold',
						'font-family':'HelveticaNeue','font-size':'14px'

					},
					style_table={
						'maxHeight': '600px',
						'overflowY': 'scroll'
					},
					sort_action = 'native',
					sort_mode = 'multi',
					row_selectable = 'single',
					# filter_action = 'native',
					style_data_conditional=[{
						'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PAM-disruption'},
						'backgroundColor': "#62c096",
						'color': 'white'
					},
					{
						'if': {'column_id': 'annotation', 'filter_query': '{annotation} eq PE3b'},
						'backgroundColor': "#62c096",
						'color': 'white'
					},
					]
			)


select_variant_to_show = dcc.Dropdown(
	options=[
		{'label': 'no input', 'value': 'no input'},
	],
	clearable=False,
	id="select_variant_to_show",
	placeholder ="select variant to show",
)


output_info = dbc.Row([
			dbc.Col(html.H5('Design Tables (Easy-Prime Output)', id='output_panel_header'),md=8),
			# dbc.Col(html.H6('Select variants to show:'),md=2,align="center"),
			dbc.Col(select_variant_to_show,md=4,align="center"),
])

init_sgRNA_table = pd.DataFrame([],columns=['chr','start','end','seq','DeepSpCas9_score','strand','target_pos','annotation'])
sgRNA_table = display_table(init_sgRNA_table,'sgRNA')


init_ngRNA_table = pd.DataFrame([],columns=['chr','start','end','seq','nick_pos','strand'])
ngRNA_table = display_table(init_ngRNA_table,'ngRNA')
init_RTT_table = pd.DataFrame([],columns=['chr','start','end','seq','RTT_length','strand'])
RTT_table = display_table(init_RTT_table,'RTT')
init_PBS_table = pd.DataFrame([],columns=['chr','start','end','seq','PBS_length','strand'])
PBS_table = display_table(init_PBS_table,'PBS')
















#--------------------------------- Input layouts --------------------------------------

chr_input = dbc.FormGroup(
	[
		dbc.Label("Chromosome:", html_for="chr_input", width=5,style ={"font-weight":"bold"}),
		dbc.Col(
			dbc.Input(
				type="text", id="chr_input", placeholder="Example: chr1"
			),
			width=6,
		),
	],
	row=True,style={"margin-top":10,"margin-left":10},
)

pos_input = dbc.FormGroup(
	[
		dbc.Label("Position:", html_for="pos_input", width=5,style ={"font-weight":"bold"}),
		dbc.Col(
			dbc.Input(
				type="text", id="pos_input", placeholder="Example: 158582552"
			),
			width=6,
		),
	],
	row=True,style={"margin-left":10},
)

variant_id_input = dbc.FormGroup(
	[
		dbc.Label("Variant ID:", html_for="variant_id_input", width=5,style ={"font-weight":"bold"}),
		dbc.Col(
			dbc.Input(
				type="text", id="variant_id_input", placeholder="Example: any_name"
			),
			width=6,
		),
	],
	row=True,style={"margin-left":10},
)


ref_input = dbc.FormGroup(
	[
		dbc.Label("Reference allele:", html_for="ref_input", width=5,style ={"font-weight":"bold"}),
		dbc.Col(
			dbc.Input(
				type="text", id="ref_input", placeholder="Example: G"
			),
			width=6,
		),
	],
	row=True,style={"margin-left":10},
)


alt_input = dbc.FormGroup(
	[
		dbc.Label("Alternative allele:", html_for="alt_input", width=5,style ={"font-weight":"bold"}),
		dbc.Col(
			dbc.Input(
				type="text", id="alt_input", placeholder="Example: A"
			),
			width=6,
		),
	],
	row=True,style={"margin-left":10},
)


VCF_input = dbc.Card(
	dbc.CardBody(
		[
			chr_input,
			pos_input,
			variant_id_input,
			ref_input,
			alt_input,
		],style={'width': '100%', 'height': 250,'padding':"0 0 0 0"},
	)
)

VCF_batch_input = dbc.Card(
	dbc.CardBody(
		[
			dcc.Textarea(
				id='vcf_batch_input',
				placeholder ='VCF format\nThe first 5 columns: chr, id, pos, ref, alt\n #comment line will be ignored',
				style={'width': '100%', 'height': 250,'padding':"0 0 0 0","font-size":"13px"},
			),
		],style={'width': '100%', 'height': 250,'padding':"0 0 0 0"},
	),
)

FASTA_input = dbc.Card(
	dbc.CardBody(
		[
			dcc.Textarea(
				id='fasta_input',
				placeholder ='FASTA format\nUse _ref to represent reference allele and _alt for alternative allele.\n>test_ref\nXXXXXXX\n>test_alt\nXXXXAXXXX',
				style={'width': '100%', 'height': 250,'padding':"0 0 0 0","font-size":"13px"},
			),
		],style={'width': '100%', 'height': 250,'padding':"0 0 0 0"},
	),
)

PrimeDesign_input = dbc.Card(
	dbc.CardBody(
		[
			dcc.Textarea(
				id='PrimeDesign_input',
				placeholder ='PrimeDesign input sequence format\nSee: https://github.com/pinellolab/PrimeDesign for more information',
				style={'width': '100%', 'height': 250,'padding':"0 0 0 0","font-size":"13px"},
			),
		],style={'width': '100%', 'height': 250,'padding':"0 0 0 0"},
	),
)

step1_info_button = html.Div(
	[
		dbc.Button("", id="step1_info_button",className="glyphicon glyphicon-info-sign btn-sm"),
		dbc.Modal(
			[
				dbc.ModalHeader("Input format selection"),
				dbc.ModalBody("Four input formats are provided. If you have a single desired edit, please input its chromosome, position, a variant name, the reference allele and the alternative allele (which is your edit). If you directly have a vcf file, you can just copy and paste it using VCF batch mode. If you only have a DNA sequence, you can input the sequence (at least 100bp) as the input_ref sequence and create a new DNA sequence that contains your desired edits as input_alt sequence. "),
				dbc.ModalFooter(
					dbc.Button("Close", id="step1_info_button_close")
				),
			],
			id="step1_info_modal",
		),
	]
)

step2_info_button = html.Div(
	[
		dbc.Button("", id="step2_info_button",className="glyphicon glyphicon-info-sign btn-sm"),
		dbc.Modal(
			[
				dbc.ModalHeader("Parameter definition"),
				dbc.ModalBody("In most cases, you can just use default searching parameters."),
				dbc.ModalFooter(
					dbc.Button("Close", id="step2_info_button_close")
				),
			],
			id="step2_info_modal",
		),
	]
)



modal = html.Div(
	[
		dbc.Button("Open", id="open-centered"),
		dbc.Modal(
			[
				dbc.ModalHeader("Header"),
				dbc.ModalBody("This modal is vertically centered"),
				dbc.ModalFooter(
					dbc.Button(
						"Close", id="close-centered", className="ml-auto"
					)
				),
			],
			id="modal-centered",
			centered=True,
		),
	]
)

def myRangeSlider(name,vmin,vmax,labelDict,defaultValue,subtitle=None):
	slider = dcc.RangeSlider(
		id=f"slider-{name}",
		# className="left-align-slider",
		min=vmin,
		max=vmax,
		value=defaultValue,
		marks=labelDict,
		step=1,
		included=True,
		allowCross=False
	)

	return html.Div(
		style={'font-weight':'bold'},
		children=[
			html.H6(subtitle,style={"margin-left":10}),
			html.Div(
				style={"color": "grey"},
				children=[html.Div(id=f'slider-output-{name}',style={'margin-bottom':10,"margin-left":10}),
				slider],
			),
			
		],
	)

def mySlider(name,vmin,vmax,labelDict,defaultValue,subtitle=None):
	slider = dcc.Slider(
		id=f"slider-{name}",
		# className="left-align-slider",
		min=vmin,
		max=vmax,
		value=defaultValue,
		marks=labelDict,
		step=1,
		included=True,
	)

	return html.Div(
		style={'font-weight':'bold'},
		children=[
			html.H6(subtitle,style={"margin-left":10}),
			html.Div(
				style={"color": "grey"},
				children=[html.Div(id=f'slider-output-{name}',style={'margin-bottom':10,"margin-left":10}),
				slider],
			),
			
		],
	)
	
def add_slider_callback(app,name,title):
	@app.callback(dash.dependencies.Output(f'slider-output-{name}', 'children'),
		[dash.dependencies.Input(f"slider-{name}", 'value')])
	def update_output(value):
		try:
			return f"{title}: [{value[0]}, {value[1]}]"
		except:
			return f"{title}: [0, {value}]"


def RTT_slider():
	labelDict={}
	labelDict[7] = {'label':'7','style': {'color': '#77b0b1'}}
	labelDict[10] = {'label':'10','style': {'color': '#77b0b1'}}
	labelDict[15] = {'label':'15','style': {'color': '#77b0b1'}}
	labelDict[20] = {'label':'20','style': {'color': '#77b0b1'}}
	labelDict[30] = {'label':'30','style': {'color': '#FF5500'}}
	labelDict[40] = {'label':'40','style': {'color': '#FF0000'}}
	labelDict[50] = {'label':'50','style': {'color': '#FF0000'}}
	labelDict[60] = {'label':'60','style': {'color': '#FF0000'}}
	subtitle = 'Reverse Transcription Template length'
	return myRangeSlider("rtt",7,60,labelDict,defaultValue=[10,20],subtitle=subtitle)


def PBS_slider():
	labelDict={}
	labelDict[7] = {'label':'7','style': {'color': '#77b0b1'}}
	labelDict[10] = {'label':'10','style': {'color': '#77b0b1'}}
	labelDict[12] = {'label':'12','style': {'color': '#77b0b1'}}
	labelDict[14] = {'label':'14','style': {'color': '#77b0b1'}}
	labelDict[17] = {'label':'17','style': {'color': '#77b0b1'}}
	# subtitle = 'Please select the maximum PBS length'
	subtitle = 'Primer Binding Sequence length'
	return myRangeSlider("pbs",7,17,labelDict,defaultValue=[10,14],subtitle=subtitle)

def ngRNA_slider():
	labelDict={}
	labelDict[0] = {'label':'0','style': {'color': '#77b0b1'}}
	labelDict[50] = {'label':'50','style': {'color': '#77b0b1'}}
	labelDict[100] = {'label':'100','style': {'color': '#77b0b1'}}
	labelDict[150] = {'label':'150','style': {'color': '#FF5500'}}
	labelDict[200] = {'label':'200','style': {'color': '#FF0000'}}
	# subtitle = 'Please select the maximum ngRNA distance to pegRNA'
	subtitle = 'nick-gRNA distance'
	# title = 'ngRNA to pegRNA distance range'
	return mySlider("ngRNA",0,200,labelDict,defaultValue=100,subtitle=subtitle)

RTT_input = dbc.Card(
	dbc.CardBody(
		[
			RTT_slider(),
		],style={'width': '100%', 'height': 150,'padding':"0 0 0 0"},
	)
)

PBS_input = dbc.Card(
	dbc.CardBody(
		[
			PBS_slider(),
		],style={'width': '100%', 'height': 150,'padding':"0 0 0 0"},
	)
)

ngRNA_input = dbc.Card(
	dbc.CardBody(
		[
			ngRNA_slider(),
		],style={'width': '100%', 'height': 150,'padding':"0 0 0 0"},
	)
)


step1_info = dbc.FormGroup(
	[
		html.H6("Step 1. Select the input format below."),
		# html.Span(className="glyphicon glyphicon-info-sign"),
		step1_info_button,
		
	],
	row=True,style={'margin-left': 0, 'margin-bottom': 0},
)

step2_info = dbc.FormGroup(
	[
		html.H6("Step 2. Choose searching parameters."),
		step2_info_button,
		
	],
	row=True,style={'margin-left': 0, 'margin-bottom': 0},
)


variant_input = html.Div([
	step1_info,
	dbc.Tabs(
		[
			dbc.Tab(VCF_input, label="VCF",tab_id='vcf_tab'),
			dbc.Tab(VCF_batch_input, label="VCF batch",tab_id='vcf_batch_tab'),
			dbc.Tab(FASTA_input, label="FASTA",tab_id='fasta_tab'),
			dbc.Tab(PrimeDesign_input, label="PrimeDesign",tab_id='PrimeDesign_tab'),
		],card=True,id="variant_input_tabs",className="nav-pills",style={"margin-bottom":10,"margin-left":10},
	),
])





parameter_input = html.Div([
	step2_info,
	dbc.Tabs(
		[
			dbc.Tab(RTT_input, label="RTT",id='rtt_input_tab'),
			dbc.Tab(PBS_input, label="PBS",id='pbs_input_tab'),
			dbc.Tab(ngRNA_input, label="ngRNA",id='ngRNA_input_tab'),
		],card=True,id="parameter_input_tabs",className="nav-pills",style={"margin-bottom":10,"margin-left":10},
	),
])

job_submission = dbc.FormGroup(
	[
		html.Button('Start', id='start_search'),
		html.Button('Examples', id='show_example',className="ml-auto mr-5"),
	],
	row=True,style={'margin-left':10,'margin-top':10,'margin-bottom':10},
)



#--------------------------------- Input panel (left) --------------------------------------

user_input = html.Div(
			id="user_input",
			children=[
				# welcome(),
				html.Div(id="parameters", style={"background-color":"white"},
					children=[variant_input,parameter_input,job_submission])
			],style={"background-color":"white","margin-top":10}
	
)

#--------------------------------- Output panel (right) --------------------------------------


Output_selection = html.Div([
	output_info,
	dbc.Tabs(
		[
			dbc.Tab(sgRNA_table, label="sgRNA table",id='sgRNA_table_tab'),
			dbc.Tab(PBS_table, label="PBS table",id='PBS_table_tab'),
			dbc.Tab(RTT_table, label="RTT table",id='RTT_table_tab'),
			dbc.Tab(ngRNA_table, label="ngRNA table",id='ngRNA_table_tab'),
			# dbc.Tab([], label="asd",id='predicted_eff',disabled=False),
			# dbc.Tab(PBS_table, label="VCF batch",id='vcf_batch_tab'),
			# dbc.Tab(RTT_table, label="FASTA",id='fasta_tab'),
			# dbc.Tab(ngRNA_table, label="PrimeDesign",id='PrimeDesign_tab'),
		],card=True,id="output_table_tabs",className="nav-pills",style={"margin-bottom":10,"margin-left":10},
	),
	html.Div(id="current_pegRNA_table",style={"margin-top":20}),
],style={"background-color":"white","margin-top":10})

#---------------------------- Visualization panel (bottom) --------------------------


Vis_loading = dcc.Loading(
	id="loading-2",
	children=[html.Div(id="vis_loading")],
	type="default",
)

vis_info = dbc.Row([
			dbc.Col(html.H5('Design Visualizations'),md=4),
			dbc.Col(Vis_loading,md=8),
])


vis_tracks = html.Div([
	vis_info,
	dbc.Tabs(
		[
			# dbc.Tab(sgRNA_table, label="sgRNA table",id='sgRNA_table_tab'),
			# dbc.Tab(PBS_table, label="PBS table",id='PBS_table_tab'),
			# dbc.Tab(RTT_table, label="RTT table",id='RTT_table_tab'),
			# dbc.Tab(ngRNA_table, label="ngRNA table",id='ngRNA_table_tab'),
			# dbc.Tab([], label="asd",id='predicted_eff',disabled=False),
			# dbc.Tab(PBS_table, label="VCF batch",id='vcf_batch_tab'),
			# dbc.Tab(RTT_table, label="FASTA",id='fasta_tab'),
			# dbc.Tab(ngRNA_table, label="PrimeDesign",id='PrimeDesign_tab'),
		],card=True,id="output_vis_tabs",className="nav-pills",style={"margin-bottom":10,"margin-left":10},
	),
],style={"background-color":"white","margin-top":10})

def add_vis_tab(name,flag,src,tab_id,track_src=None,view_location=None):
	image = html.Img(src=src)
	if flag == "iframe":
		iframe = html.Iframe(title=name,srcDoc=get_PE_vis_tracks(view_location,track_src),width="100%",height=400)
		# iframe = html.Iframe(title=name,srcDoc=get_PE_vis_tracks())
		return dbc.Tab([iframe,image],label=name,tab_id=tab_id)
		# return dbc.Tab(iframe,label=name,tab_id=tab_id)
	return dbc.Tab(image,label=name,tab_id=tab_id)

PE_vis_track_template = """
<html>
<body>

<script src="https://proteinpaint.stjude.org/bin/proteinpaint.js" charset="utf-8"></script>

<div id=a style="margin:10px"></div>

<script>
runproteinpaint({
	host:'https://proteinpaint.stjude.org',
	holder:document.getElementById('a'),
	parseurl:true,
	block:true,
	nobox:1,
	noheader:1,
	genome:'hg19',
	position:'{view_location}',
	nativetracks:'RefGene',
	tracks:[   
		{
		type:"bigwig",
	  "url":"http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw",
			"name":"UCSC phyloP 100ways",
			"height":100
		},
		{track_src}
	]
}).then(_=>{_.block.highlight_1basedcoordinate('{view_location}')})
</script>

</body>
</html>

"""


def get_PE_vis_tracks(view_location,track_src):
	print (view_location,track_src)
	PE_vis_tracks = PE_vis_track_template.replace("{track_src}",track_src)
	PE_vis_tracks = PE_vis_tracks.replace("{view_location}",view_location)

	return PE_vis_tracks


import dash
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
import plotly.graph_objs as go
import dash_table
import plotly_express as px
from dash.dependencies import Input, Output, State
import re
import numpy as np
from sklearn.cluster import KMeans
import plotly.figure_factory as ff
import os
import base64
import datetime
from datetime import datetime as dt
import pathlib

def get_version():
	return "v1.0"

#--------------------------------- Right Column --------------------------------------


def encode_fig(pegRNA_id=None):
	print ("#######????????????#######",pegRNA_id)
	fig = "img/rs737092_chr20_55990370_55990390_TACCTCCTGGGCCTCCCGCA_candidate_162.top_pegRNA.png"
	with open(fig, "rb") as image_file:
		img_string = base64.b64encode(image_file.read())
	return "data:image/png;base64,%s"%(img_string.decode("utf-8"))

def vis_PE_design(app,pegRNA_list):
	## TODO add option when pegRNA_list=None
	print (list(pegRNA_list[0].values())[0])
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
		return df.to_dict('records')
	print (df.head())
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
























#--------------------------------- Left Column --------------------------------------


def buttons():

	return html.Div([html.Button('Start', id='start_search',style={'margin-left':'8%'}),html.Button('Tutorial', id='help_menu',style={"float":"right",'margin-right':'8%'})])


def help():

	return html.Div([html.Button('Help', id='help_menu')])

intro="Easy-Prime is a machine learning based tool for prime editing gRNA design. Please input your desire edits in VCF format and click start.  Additionally, you can play with the pegRNA searching parameters or choose a different genome in the bottom."

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
		value='VCF format\nThe first 5 columns: chr, id, pos, ref, alt\n #comment line will be ignored',
		style={'width': '100%', 'height': 100,'padding':"0 0 0 0","font-size":"13px"},
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
	return mySlider("pbs",7,17,labelDict,defaultValue=12,title=title,subtitle=subtitle,app=app)

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
	labelDict[30] = {'label':'30','style': {'color': '#77b0b1'}}
	labelDict[40] = {'label':'40','style': {'color': '#77b0b1'}}
	labelDict[50] = {'label':'50','style': {'color': '#FF5500'}}
	labelDict[60] = {'label':'60','style': {'color': '#FF5500'}}
	labelDict[70] = {'label':'70','style': {'color': '#FF0000'}}
	labelDict[80] = {'label':'80','style': {'color': '#FF0000'}}
	subtitle = 'Please select the maximum RTT length'
	subtitle = ''
	title = 'RTT length range'
	return mySlider("rtt",7,80,labelDict,defaultValue=40,title=title,subtitle=subtitle,app=app)

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



import utils

# exec(open("test_dash.py").read())
import imp  
imp.reload(utils)
from utils import *

app = dash.Dash(
	__name__,
	meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

server = app.server
app.config.suppress_callback_exceptions = True

## a fake list
pegRNA_list=[
{'label': 'rs737092', 'value': 'rs737092'}
]

app.layout = html.Div(
	id="app-container",
	children=[
		# Banner
		html.Div(
			id="banner",
			className="banner",
			children=[html.Img(src=app.get_asset_url("logo.png")),html.H4("Easy-Prime %s"%(get_version()))],
		),
		# Left column
		html.Div(
			id="left-column",
			className="four columns",
			children=[
				# welcome
				welcome(),
				html.Div([
					dcc.Input(
					id = 'jid',
					value = '',
					)
				], style= {'display': 'none'} # <-- This is the line that will be changed by the dropdown callback
				),
				# input
				# variant_input(),

				# parameter
				html.Div(id="parameters", style={"background-color":"white"},
					children=[variant_input(),buttons(),RTT_slider(app),PBS_slider(app),ngRNA_slider(app),genome_input(app)])




			],
	
		),
		# Right column
		html.Div(
			id="right-column",
			className="eight columns",
			children=[
				# vis pe design
				html.Div(
					id="vis_top_pegRNA",
					style={"background-color":"#FFFFFF","margin-top":"20px",'font-weight':'bold','font-size':'20px'},
					children=[
						html.Div(style={"vertical-align": "middle","display": "inline-block"},children=[
							'PE design visualization for:',
							dcc.Input(
								id = 'select_pegRNA_id',
								value = '',
								debounce=True,
							)
						]),
						html.Div([html.Img(id='PE_design_figure')]),		

					],
				),
				# vis_PE_design(app,pegRNA_list=pegRNA_list),
				# PE table title
				html.Div(style={"margin-top":"20px","margin-bottom":"10px",'font-weight':'bold','font-size':'20px'},children=['Sequences of top pegRNA and ngRNA combinations are shown below.']),
				# show PE sequence table
				# show_PE_table(app),
				dash_table.DataTable(
					id='pegRNA-table',
					columns=[
						{'name': i, 'id': i, 'deletable': False} for i in ['sample_ID','type','seq','chr','start','end','strand']
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
				),

				html.Div([
					html.Div(html.A('Download top-pegRNAs.csv', id='topX', download="top_pegRNAs.csv",href="",target="_blank"),style={'width': '20%','display': 'inline-block'}),
					html.Div(html.A('Download all-pegRNAs.csv', id='X_p', download="X_p_pegRNAs.csv",href="",target="_blank"),style={'width': '20%','display': 'inline-block'}),
					html.Div(html.A('Download raw-pegRNAs.csv', id='rawX', download="rawX_pegRNAs.csv",href="",target="_blank"),style={'width': '20%','display': 'inline-block'})
				]),


			],
		),
	],
)

## main easy prime

@app.callback(
	[
	Output("jid", "value"),
	Output("pegRNA-table", "data"),

	], 
	[
	Input('start_search', 'n_clicks'),
	Input('slider-pbs', 'value'), # PBS
	Input('slider-rtt', 'value'), # RTT
	Input('slider-ngRNA', 'value'), # ngRNA
	Input('genome', 'value'), # genome
	Input('variants', 'value'), # variants


	]
)
def search_pegRNA(n_clicks,pbs_value,rtt_value,ngRNA_value,genome,input_variants):
	"""search pegRNA, show table, and vis the top one

	users can change visualization by searching different
	
	"""

	## get parameters
	print ("search_pegRNA")
	if n_clicks == None:
		return [None,None]
	parameters = get_parameters("config.yaml")
	print_parameters(parameters)



	jid=str(uuid.uuid4()).split("-")[-1]
	input_variants = StringIO(input_variants)
	parameters['min_PBS_length'] = pbs_value[0]
	parameters['max_PBS_length'] = pbs_value[1]
	parameters['min_RTT_length'] = rtt_value[0]
	parameters['max_RTT_length'] = rtt_value[1]
	parameters['max_ngRNA_distance'] = ngRNA_value
	summary,df_top,df_all,X_p = easy_prime_main(input_variants,jid,parameters)

	## initiate data table



	return [jid,df_top.to_dict('records')]




@app.callback(
	
	Output("PE_design_figure", "src"),
	[
	Input('jid', 'value'),
	Input('select_pegRNA_id', 'value'),
	]
)
def update_vis_pegRNA(jid,select_pegRNA_id):
	"""search pegRNA, show table, and vis the top one

	users can change visualization by searching different
	
	"""
	print ("my values")
	print ("jid is ",jid)
	print ("select_pegRNA_id is ",select_pegRNA_id)
	if jid == None:
		print ("jid is ",jid)
		return init_fig()
	# if select_pegRNA_id == "":
	# 	print ("select_pegRNA_id is ",select_pegRNA_id)
	# 	return init_fig()

	parameters = get_parameters("config.yaml")
	df_all = pd.read_csv("results/%s_rawX_pegRNAs.csv.gz"%(jid),index_col=0)
	df_top = pd.read_csv("results/%s_topX_pegRNAs.csv"%(jid),index_col=0)
	print (df_all)
	print (df_top)
	if select_pegRNA_id == "":
		select_pegRNA_id = df_top.index.tolist()[0]
	tmp_df_name = jid+".vis.df"
	df_all.loc[select_pegRNA_id].to_csv("results/%s"%(tmp_df_name))

	# img = vis_pegRNA(,parameters['genome_fasta'])
	img = vis_pegRNA("results/%s"%(tmp_df_name),genome_fasta=parameters['genome_fasta'],out_file_name=tmp_df_name)
	return img

## download data table (need jid)



if __name__ == '__main__':
	app.run_server(debug=False,host='0.0.0.0',port=8051)

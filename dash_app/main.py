import utils

# exec(open("test_dash.py").read())
import imp  
imp.reload(utils)
from utils import *
server = Flask(__name__)
app = dash.Dash(
	__name__,
	meta_tags=[{"name": "viewport", "content": "width=800px, initial-scale=1"}],
	server=server,
)

# server = app.server
app.config.suppress_callback_exceptions = True
@server.route("/results/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory("results", path, as_attachment=True)


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
						html.Div(children=[
							'PE design visualization for:',
							dcc.Input(
								id = 'select_pegRNA_id',
								value = '',
								style = {"width":"500px"}
							)
						]),
						html.Div([html.Img(id='PE_design_figure')]),		

					],
				),
				# vis_PE_design(app,pegRNA_list=pegRNA_list),
				# PE table title
				html.Div(style={"margin-top":"20px","margin-bottom":"10px",'font-weight':'bold','font-size':'20px'},children=['Sequences of top pegRNA and ngRNA combinations are shown below.']),
				# show PE sequence table
				dash_table.DataTable(
					id='pegRNA-table',
					columns=[
						{'name': i, 'id': i, 'deletable': False} for i in ['sample_ID','type','seq','chr','start','end','strand','predicted efficiency']
					],
					style_cell={
						'minWidth': '0px', 'maxWidth': '20px',
						# 'width':'12%',
						'whiteSpace': 'normal',
						'height': 'auto',
						# 'whiteSpace': 'normal',
						# 'text-overflow': 'ellipsis',
					},	
					style_data={
						'whiteSpace': 'normal',
						'height': 'auto',
					},	
					style_cell_conditional=[
						{'if': {'column_id': 'type'},
						 'width': '7%'},
						{'if': {'column_id': 'strand'},
						 'width': '7%'},
						{'if': {'column_id': 'start'},
						 'width': '10%'},
						{'if': {'column_id': 'end'},
						 'width': '10%'},
						{'if': {'column_id': 'chr'},
						 'width': '10%'},
						{'if': {'column_id': 'predicted efficiency'},
						 'width': '10%'},
					],
					style_header={
						'text-align':'center',
						'font-weight': 'bold',


					},						
					editable=True,
					filter_action="native",
					page_action="native",
					page_current= 0,
					page_size= 13,
				),


				html.Div(id='download_data'),
				

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
	# print ("search_pegRNA")
	if n_clicks == None:
		return [None,None]
	parameters = get_parameters("config.yaml")
	# print_parameters(parameters)



	jid=str(uuid.uuid4()).split("-")[-1]
	# input_variants = StringIO(input_variants)
	parameters['min_PBS_length'] = pbs_value[0]
	parameters['max_PBS_length'] = pbs_value[1]
	parameters['min_RTT_length'] = rtt_value[0]
	parameters['max_RTT_length'] = rtt_value[1]
	parameters['max_ngRNA_distance'] = ngRNA_value
	summary,df_top,df_all,X_p = easy_prime_main(input_variants,jid,parameters)
	df_top['predicted_efficiency'] = ["{0:.2f}%".format(x * 100) for x in df_top['predicted_efficiency']]
	df_top=df_top.rename(columns = {'predicted_efficiency':'predicted efficiency'})
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
	# print ("my values")
	# print ("jid is ",jid)
	# print ("select_pegRNA_id is ",select_pegRNA_id)
	if jid == None:
		# print ("jid is ",jid)
		return init_fig()
	# if select_pegRNA_id == "":
	# 	print ("select_pegRNA_id is ",select_pegRNA_id)
	# 	return init_fig()

	parameters = get_parameters("config.yaml")
	df_all = pd.read_csv("results/%s_rawX_pegRNAs.csv.gz"%(jid),index_col=0)
	df_top = pd.read_csv("results/%s_topX_pegRNAs.csv"%(jid),index_col=0)
	# print (df_all)
	# print (df_top)
	if select_pegRNA_id == "":
		select_pegRNA_id = df_top.index.tolist()[0]
	tmp_df_name = jid+".vis.df"
	df_all.loc[select_pegRNA_id].to_csv("results/%s"%(tmp_df_name))

	# img = vis_pegRNA(,parameters['genome_fasta'])
	file_name = "results/%s.fa"%(jid)
	if os.path.isfile(file_name):
		genome_fasta = file_name
	else:
		genome_fasta = parameters['genome_fasta']
	img = vis_pegRNA("results/%s"%(tmp_df_name),genome_fasta=genome_fasta,out_file_name=tmp_df_name)
	return img

## download data table (need jid)
def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "results/{}".format(urlquote(filename))
    return html.A(filename, href=location,target="_blank")
@app.callback(Output('download_data', 'children'), [Input('jid', 'value')])
def download_predict_data(jid):
	if jid == None:
		return None
		
	return [
				"Download files:",
				html.Ul([
					html.Li(file_download_link("%s_rawX_pegRNAs.csv.gz"%(jid))),
					html.Li(file_download_link("%s_topX_pegRNAs.csv"%(jid))),
					html.Li(file_download_link("%s_X_p_pegRNAs.csv.gz"%(jid)))
				]),
			]
'''
@app.callback(Output('topX', 'href'),
	[Input('jid', 'value')]
	)
def download_topX(jid):
	df = pd.read_csv("results/%s_topX_pegRNAs.csv"%(jid),index_col=0)
	csv_string = df.to_csv(index=False, encoding='utf-8')
	# print (df.head())
	# print ("#### printing topX ####")
	# csv_string = "data:text/csv;charset=utf-8,%EF%BB%BF" + urllib.parse.quote(csv_string)
	# csv_string = "data:text/csv;charset=utf-8,%EF%BB%BF" + urllib.parse.quote(csv_string)
	csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(csv_string)
	# print (csv_string)
	return csv_string
'''

if __name__ == '__main__':
	app.run_server(debug=False,host='0.0.0.0',port=8051)

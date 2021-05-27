import utils

# exec(open("track_view_version.py").read())
import imp  
imp.reload(utils)
from utils import *
server = Flask(__name__)
# external_scripts =['https://proteinpaint.stjude.org/bin/proteinpaint.js']
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(
	__name__,
	# meta_tags=[{"name": "viewport", "content": "width=800px, initial-scale=1"}],
	meta_tags=[{"name": "viewport"}],
	server=server,
	serve_locally=True,
	# external_scripts=external_scripts,
	external_stylesheets=external_stylesheets,
)


app.config.suppress_callback_exceptions = True

@server.route("/results/<path:path>")
def download(path):
	"""Serve a file from the upload directory."""
	return send_from_directory("results", path, as_attachment=True)

Banner = html.Div(
			id="banner",
			className="banner",
			children=[html.Img(src=app.get_asset_url("logo.png")),html.H4("Easy-Prime %s"%(get_version()))],
)

Parameter_selection2 = dbc.Row([
	dbc.Col()
])




Parameter_selection = html.Div(
			id="Parameter_selection",
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
	
)

Output_selection = html.Div(
			id="sdf",
			className="sdf",
			children=[],
)

Track_view = html.Div(id="PE_track_view")

app.layout = dbc.Container(
	[
		dbc.Row([dbc.Col(Banner)]),
		dbc.Row( 
			[
				dbc.Col(Parameter_selection, md=4),
				dbc.Col(Output_selection, md=8),

			],
			align="top",
		),
		dbc.Row(Track_view),
	],
	fluid=True,
)

test_ttt = """
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
	position:'chr12:25334648-25426944',
	nativetracks:'RefGene',
	tracks:[   
        {
	    type:"bigwig",
      "url":"http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw",
			"name":"UCSC phyloP 100ways",
			"height":100
		},
	]
})
</script>

</body>
</html>

"""

@app.callback(Output('PE_track_view', 'children'), Input('start_search', 'n_clicks'))
def embed_iframe(value):
	if not value:
		return value

	return html.Iframe(id="ifm",title="ttttt",srcDoc=test_ttt,width=1600,height=800)

'''
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
# if __name__ == '__main__':
app.run_server(debug=False,host='0.0.0.0',port=9866)

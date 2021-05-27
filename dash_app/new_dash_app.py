import utils2
import layout

# exec(open("new_dash_app.py").read())
import imp  
imp.reload(utils2)
imp.reload(layout)
from utils2 import *
from layout import *
server = Flask(__name__)
external_stylesheets = [dbc.themes.BOOTSTRAP]
app = dash.Dash(
	__name__,
	meta_tags=[{"name": "viewport"}],
	server=server,
	serve_locally=True,
	external_stylesheets=external_stylesheets,
)


app.config.suppress_callback_exceptions = True

@server.route("/results/<path:path>")
def download(path):
	"""Serve a file from the upload directory."""
	return send_from_directory("results", path, as_attachment=True)


#--------------------------------- Header panel (top) --------------------------------------

Banner = html.Div(
			id="banner",
			className="banner",
			children=[html.Img(src=app.get_asset_url("logo.png")),html.H4("Easy-Prime %s"%(get_version()))],
)

Status = html.Div(
	[
		dbc.Button("Check Running Status", id="error_info_button"),
		dbc.Modal(
			[
				dbc.ModalHeader("No status found!",id="status_header"),
				dbc.ModalBody("Possibly you just opened this url",id="status_body"),
				dbc.ModalFooter(
					dbc.Button("Close", id="error_info_button_close")
				),
			],
			id="error_info_modal",
		),
	]
)
Status_loading = dcc.Loading(
	id="loading-1",
	children=[Status],
	type="circle",
)




#----------------------- MAIN layout ---------------------------

app.layout = dbc.Container(
	[
		# header, showing logo and web page running status
		dbc.Row([
			dbc.Col(Banner,md=10),
			dbc.Col(Status_loading,md=2,align="center"),
		],style={"background-color":"white"}),
		
		# main content
		# left is input parameter selection
		# right is table output
		
		dbc.Row([
			dbc.Col(user_input, md=4),
			dbc.Col(Output_selection, md=8),
			],
		align="top",
		),
		
		# bottom is visualization
		dbc.Row([
			dbc.Col(vis_tracks,md=12),
			
		]),
		
		# hidden elements
		html.Div(id='store-rawX', style={'display': 'none'}),
		html.Div(id='store-X_p', style={'display': 'none'}),
		html.Div(id='store-sgRNA-table', style={'display': 'none'}),
		html.Div(id='store-PBS-table', style={'display': 'none'}),
		html.Div(id='store-RTT-table', style={'display': 'none'}),
		html.Div(id='store-ngRNA-table', style={'display': 'none'}),
		html.H6("",id="current_jid",hidden=True),
	],
	fluid=True,
)

#----------------------- helpful callbacks ---------------------------

# load examples
@app.callback(
	[
	Output('chr_input', 'value'), 
	Output('pos_input', 'value'), 
	Output('variant_id_input', 'value'), 
	Output('ref_input', 'value'), 
	Output('alt_input', 'value'), 
	Output('vcf_batch_input', 'value'), 
	Output('fasta_input', 'value'), 
	Output('PrimeDesign_input', 'value'), 
	],
	[
	Input('show_example','n_clicks'),
	]
)
def load_examples(value):
	if not value:
		return ""
	return "chr1","158582552",'rs2251964',"G","A","## comment line, will be ignored\nchr9	110184636	FIG5G_HEK293T_HEK3_6XHIS	G	GCACCATCATCACCATCAT\nchr1	185056772	FIG5E_U2OS_RNF2_1CG	G	C\nchr1	173878832	rs5878	T	C\nchr11	22647331	FIG3C_FANCF_7AC_PE3B	T	G\nchr19	10244324	EDFIG5B_DNMT1_dPAM	G	T\n",">rs2251964_ref\nGTTACCAAAGCAAATGACATCTTGTGAAAGGGGAGGTCTGAAAAAAAAAAACAAGTGGGTGGGTTTTTTCAAAGTAGGCCACCGGGCCTGAGATGACCAGAATTCAAATTAGGATGACAGTGTAGTAGGGGAAGCAACCAGAATCGGACCT\n>rs2251964_alt\nGTTACCAAAGCAAATGACATCTTGTGAAAGGGGAGGTCTGAAAAAAAAAAACAAGTGGGTGGGTTTTTTCAAAGTAGGCCACCGGGCCTGAGATAACCAGAATTCAAATTAGGATGACAGTGTAGTAGGGGAAGCAACCAGAATCGGACCT\n",">test_SNV\nGCCTGTGACTAACTGCGCCAAAACGGCCTGTGACTAACTGCGCCAGCCTGTGACTAACTGCGCCAAAACGAAACG(T/A)GCCTGGCCTGTGACTAACTGCGCCAAAACGTGACTAACTGCGCCAAAACGCTTCCAATCCCCTTATCCAATTTA\n>test_insertion\nGCCTGTGCCTGTGACTAACTGCGCCAAAACGGAGCCTGTGACTAACTGCGCCAAAACGCTAACTGCGCCAAAACGT(+CTT)CTTCCGCCTGGCCTGTGACTAACTGCGCCAAAACGTGACTAACTGCGCCAAAACGAATCCCCTTATCCAATTTA\n>test_deletion\nGCCTGTGACTAGCCTGTGACTAACTGCGCCAAAACGACTGCGCGCCTGTGACTAACTGCGCCAAAACGCAAAAC(-GTCT)TCCAATCGCCTGTGACTAACTGCGCCAAAACGCCCTTATCCGCCTGTGACTAACTGCGCCAAAACGAATTTA\n"


# select variants and update sgRNA table
@app.callback(
	[
	Output('sgRNA-table', 'data'),
	Output('sgRNA-table', 'selected_rows'),
	Output('store-sgRNA-table', 'children'),
	],
	[
	Input('select_variant_to_show','value'),
	],
	state = [
	State('store-rawX', 'children'), 
	State('store-X_p', 'children'), 
	# add parameter states
	]
)
def update_sgRNA_table(sample_ID,rawX,X_p):
	if not sample_ID:
		# return init_sgRNA_table.to_dict('records'),None,None
		return None
	# print ("inside update_sgRNA_table")
	rawX = pd.read_json(rawX, orient='split')
	# print (rawX.head())
	X_p = pd.read_json(X_p, orient='split')
	# print (X_p.head())
	# sgRNA_df,ID_list = to_sgRNA_table(rawX,X_p,sample_ID)
	sgRNA_df = to_sgRNA_table(rawX,sample_ID)
	selected_sgRNA = [0]
	sgRNA_df_stored = sgRNA_df.copy()
	
	sgRNA_df_stored['name'] = sgRNA_df_stored.chr+"_"+sgRNA_df_stored.start.astype(str)+"_"+sgRNA_df_stored.end.astype(str)+"_"+sgRNA_df_stored.seq
	return sgRNA_df.to_dict('records'),selected_sgRNA,sgRNA_df_stored.to_json(date_format='iso', orient='split')


# select sgRNA and update PBS,RTT,ngRNA tables
@app.callback(
	[
	# Output('sgRNA-table', 'data'),
	Output('PBS-table', 'data'),
	Output('RTT-table', 'data'),
	Output('ngRNA-table', 'data'),
	Output('store-PBS-table', 'children'),
	Output('store-RTT-table', 'children'),
	Output('store-ngRNA-table', 'children'),
	Output('PBS-table', 'selected_rows'),
	Output('RTT-table', 'selected_rows'),
	Output('ngRNA-table', 'selected_rows'),
	],
	[
	Input('sgRNA-table','selected_rows'),
	],
	state = [
	State('store-rawX', 'children'), 
	State('store-sgRNA-table', 'children'), 
	# add parameter states
	]
)
def update_other_3_tables(sgRNA_table_index,rawX,sgRNA_df):
	if not sgRNA_table_index:
		# return init_PBS_table.to_dict('records'),init_RTT_table.to_dict('records'),init_ngRNA_table.to_dict('records'),None,None,None,None,None,None
		return None
	sgRNA_df = pd.read_json(sgRNA_df, orient='split')
	# print (sgRNA_df.head())
	sgRNA_location = sgRNA_df.at[sgRNA_table_index[0],'name']
	rawX = pd.read_json(rawX, orient='split')
	RTT_df,selected_RTT = to_RTT_table(rawX,sgRNA_location)
	# print (RTT_df,selected_RTT)
	ngRNA_df,selected_ngRNA = to_ngRNA_table(rawX,sgRNA_location)
	# print (ngRNA_df,selected_ngRNA)
	PBS_df,selected_PBS = to_PBS_table(rawX,sgRNA_location)
	# print (PBS_df,selected_PBS)

	ngRNA_df_stored = ngRNA_df.copy()
	ngRNA_df_stored['name'] = ngRNA_df_stored.chr+"_"+ngRNA_df_stored.start.astype(str)+"_"+ngRNA_df_stored.end.astype(str)+"_"+ngRNA_df_stored.seq

	PBS_df_stored = PBS_df.copy()
	PBS_df_stored['name'] = PBS_df_stored.chr+"_"+PBS_df_stored.start.astype(str)+"_"+PBS_df_stored.end.astype(str)+"_"+PBS_df_stored.seq

	RTT_df_stored = RTT_df.copy()
	RTT_df_stored['name'] = RTT_df_stored.chr+"_"+RTT_df_stored.start.astype(str)+"_"+RTT_df_stored.end.astype(str)+"_"+RTT_df_stored.seq

	return PBS_df.to_dict('records'),RTT_df.to_dict('records'),ngRNA_df.to_dict('records'),PBS_df_stored.to_json(date_format='iso', orient='split'),RTT_df_stored.to_json(date_format='iso', orient='split'),ngRNA_df_stored.to_json(date_format='iso', orient='split'),[selected_PBS],[selected_RTT],[selected_ngRNA]

# given sgRNA,PBS,RTT,ngRNA selection, update track view
@app.callback(
	[
	Output('output_vis_tabs', 'children'),
	Output('output_vis_tabs', 'active_tab'),
	Output('current_pegRNA_table', 'children'),
	Output('vis_loading', 'children'),
	],
	[
	Input('sgRNA-table', 'selected_rows'),
	Input('PBS-table', 'selected_rows'),
	Input('RTT-table', 'selected_rows'),
	Input('ngRNA-table', 'selected_rows'),
	],
	state = [
	State('store-rawX', 'children'), 
	State('store-sgRNA-table', 'children'),
	State('store-PBS-table', 'children'),
	State('store-RTT-table', 'children'),
	State('store-ngRNA-table', 'children'),
	State('current_jid','children'),
	State('variant_input_tabs', 'active_tab'), 
	State('output_vis_tabs', 'children'), 
	State('select_variant_to_show', 'value'), 
	# add parameter states
	]
)
def update_track_view(sgRNA_table_index,PBS_table_index,RTT_table_index,ngRNA_table_index,rawX,sgRNA_df,PBS_df,RTT_df,ngRNA_df,jid,active_tab,vis_tab,variant_id):
	if not PBS_table_index:
		return None
	rawX = pd.read_json(rawX, orient='split')
	sgRNA_df = pd.read_json(sgRNA_df, orient='split')
	PBS_df = pd.read_json(PBS_df, orient='split')
	RTT_df = pd.read_json(RTT_df, orient='split')
	ngRNA_df = pd.read_json(ngRNA_df, orient='split')
	sgRNA_location = sgRNA_df.at[sgRNA_table_index[0],'name']
	PBS_location = PBS_df.at[PBS_table_index[0],'name']
	RTT_location = RTT_df.at[RTT_table_index[0],'name']
	ngRNA_location = ngRNA_df.at[ngRNA_table_index[0],'name']
	
	ID1 = rawX[rawX.location_name==sgRNA_location].sample_ID.tolist()
	ID2 = rawX[rawX.location_name==PBS_location].sample_ID.tolist()
	ID3 = rawX[rawX.location_name==RTT_location].sample_ID.tolist()
	ID4 = rawX[rawX.location_name==ngRNA_location].sample_ID.tolist()
	# print (ID1)
	# print (ID2)
	# print (ID3)
	# print (ID4)
	ID = set(ID1).intersection(ID2,ID3,ID4)
	ID = list(ID)[0]
	# print ("vis:",ID)
	vis_df = rawX[rawX.sample_ID==ID]
	
	vis_df['vis_name'] = vis_df.target_pos.astype(str)+"_"+vis_df.PBS_length.astype(str)+"_"+vis_df.RTT_length.astype(str)+"_"+vis_df.nick_pos.astype(str)
	view_location = "%s:%s-%s"%(rawX.CHROM[0],rawX.POS[0],rawX.POS[0]+len(rawX.REF[0]))
	
	vis_name = variant_id+"_"+vis_df.vis_name[0]
	tab_id = "tab-%s"%(len(vis_tab))
	

	filename = "results/%s.vis.csv"%(jid)
	vis_df_copy = vis_df.copy()
	vis_df_copy['predicted_efficiency'] = vis_df_copy['predicted_efficiency']/100
	vis_df_copy.to_csv(filename,index=False)
	img_src = vis_pegRNA_png(filename,jid)
	# vis_tab.append(add_vis_tab(vis_name,"img",src,tab_id))
	# print (vis_tab.label)
	# print (vis_tab[0].label)
	if active_tab in ['vcf_tab','vcf_batch_tab']:
		vis_bed = "%s_%s"%(jid,len(vis_tab))
		print (vis_bed)
		track_src = df2bedjs(vis_df,vis_bed)
		# print ("showing both iframe and image")
		vis_tab.append(add_vis_tab(vis_name,"iframe",img_src,tab_id,track_src = track_src,view_location=view_location))
	else:
		vis_tab.append(add_vis_tab(vis_name,"",img_src,tab_id))
	return vis_tab,tab_id,get_current_selection_table(vis_df),None
#----------------------- main callbacks ---------------------------

# first create a new jid, then jid trigger easy_prime
@app.callback(
	Output("current_jid", "children"),
	Input('start_search','n_clicks')
)
def set_jid(value):
	if not value:
		return ""
	return get_uid()

@app.callback(
	[
	Output('store-rawX', 'children'),
	Output('store-X_p', 'children'),
	Output('select_variant_to_show', 'options'),
	Output('select_variant_to_show', 'value'),
	# Output('sgRNA-table', 'data'),
	# Output('PBS-table', 'data'),
	# Output('RTT-table', 'data'),
	# Output('ngRNA-table', 'data'),
	# Output('sgRNA-table', 'selected_rows'),
	# Output('PBS-table', 'selected_rows'),
	# Output('RTT-table', 'selected_rows'),
	# Output('ngRNA-table', 'selected_rows'),
	Output("status_header", "children"),
	Output("status_body", "children"),
	],
	[
	Input('current_jid','children'),
	],
	state = [
	State('variant_input_tabs', 'active_tab'), 
	State('chr_input', 'value'), 
	State('pos_input', 'value'), 
	State('variant_id_input', 'value'), 
	State('ref_input', 'value'), 
	State('alt_input', 'value'), 
	State('vcf_batch_input', 'value'), 
	State('fasta_input', 'value'), 
	State('PrimeDesign_input', 'value'), 
	State('slider-pbs', 'value'), 
	State('slider-rtt', 'value'), 
	State('slider-ngRNA', 'value'), 
	# add parameter states
	]
)
def start_easy_prime(jid,active_tab,chr,pos,variant_id,ref,alt,vcf_batch,fasta_input,PrimeDesign_input,pbs_value,rtt_value,ngRNA_value):
	if jid == "":
		return None
	# print (jid,active_tab)

	status_header = "No errors found."
	finished_message = "Your result should be finished. Job ID is: %s. If you didn't see any results in the table or track visualization below, please email Yichao.Li@stjude.org."%(jid)
	# check input
	vcf,input_error_flag,error_message = check_and_convert_input(active_tab,chr,pos,variant_id,ref,alt,vcf_batch,fasta_input,PrimeDesign_input,jid)
	print (vcf)
	if input_error_flag:
		status_header = "Input format error!"
		print (error_message)
		# return None,None,None,None,None,None,status_header,error_message
		return None
	
	
	# start easy_prime
	try:
		parameters = get_parameters("config.yaml")
		parameters['min_PBS_length'] = pbs_value[0]
		parameters['max_PBS_length'] = pbs_value[1]
		parameters['min_RTT_length'] = rtt_value[0]
		parameters['max_RTT_length'] = rtt_value[1]
		parameters['max_ngRNA_distance'] = ngRNA_value
		summary,df_top,df_all,X_p = run_easy_prime_backend(vcf,jid,parameters)
	except Exception as e:
		status_header = "Easy-prime running error!"
		print (e)
		return None
		# return None,None,None,None,None,None,status_header,error_message
	
	# read easy_prime output
	easy_prime_error_flag,error_message,rawX,X_p = read_easy_prime_output(jid)
	options,sample_ID =  get_options_dict(jid)
	# print (options)
	# print (sample_ID)
	if easy_prime_error_flag:
		status_header = "Parsing Easy-Prime output error!"
		return None
		# return None,None,None,None,None,None,status_header,error_message

	# print ("updating Main")
	return rawX.to_json(date_format='iso', orient='split'),X_p.to_json(date_format='iso', orient='split'),options,sample_ID,status_header,finished_message




#--------------------------------- slider range callbacks --------------------------------------

add_slider_callback(app,"rtt",'RTT length range')
add_slider_callback(app,"pbs",'PBS length range')
add_slider_callback(app,"ngRNA",'ngRNA distance range')

#--------------------------------- slider range callbacks (end) --------------------------------------

#--------------------------------- Info button callbacks --------------------------------------


@app.callback(
	Output("step1_info_modal", "is_open"),
	[Input("step1_info_button", "n_clicks"), Input("step1_info_button_close", "n_clicks")],
	[State("step1_info_modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
	if n1 or n2:
		return not is_open
	return is_open


@app.callback(
	Output("step2_info_modal", "is_open"),
	[Input("step2_info_button", "n_clicks"), Input("step2_info_button_close", "n_clicks")],
	[State("step2_info_modal", "is_open")],
)
def toggle_modal(n1, n2, is_open):
	if n1 or n2:
		return not is_open
	return is_open

@app.callback(
	Output("error_info_modal", "is_open"),
	[Input("error_info_button", "n_clicks"), Input("error_info_button_close", "n_clicks")],
	[State("error_info_modal", "is_open")],
)
def toggle_modal(n1,n2, is_open):
	if n1 or n2:
		return not is_open
	return is_open

#--------------------------------- Info button callbacks (end) --------------------------------------

if __name__ == '__main__':
	app.run_server(debug=True,host='0.0.0.0',port=9866,use_reloader=False)

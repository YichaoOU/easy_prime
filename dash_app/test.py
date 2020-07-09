import utils

# exec(open("test.py").read())
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
			className="three columns",
			children=[
				# welcome
				welcome(),

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
			className="nine columns",
			children=[
				# vis pe design
				vis_PE_design(app,pegRNA_list=pegRNA_list),
				# PE table title
				html.Div(style={"margin-top":"20px","margin-bottom":"10px",'font-weight':'bold','font-size':'20px'},children=['Sequences of top pegRNA and ngRNA combinations are shown below.']),
				# show PE sequence table
				show_PE_table(app),
			],
		),
	],
)

if __name__ == '__main__':
	app.run_server(debug=False,host='0.0.0.0',port=8051)

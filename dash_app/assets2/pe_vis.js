runproteinpaint({
	host:'https://proteinpaint.stjude.org',
	holder:document.getElementById('PE_track_view'),
	parseurl:true,
	block:true,
	nobox:1,
	noheader:1,
	genome:'hg19',
	position:'chr12:25334648-25426944',
	nativetracks:'RefGene',
	tracks:[   
		{
			"type":"bigwig",
			"url":"http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw",
			// alternatively, load from a local file on the ProteinPaint server
			//"file":"hg19/hg19.100way.phastCons.bw",
			"name":"UCSC phyloP 100ways",
			"height":100
		},
	]
})

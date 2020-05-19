
from .utils import *
	
class sgRNA:
	def __init__(self,chr,start,end,seq,strand,cut_position,target_pos,ref,alt,user_defined_name_variant_name,dist_dict,strand_dict,candidate_pegRNA_df):
		"""
		
		search key is chr_start_end_seq because seq can be duplicated, this is unique
		
		"""
		self.chr = chr
		self.start = start
		self.end = end
		self.seq = seq
		self.user_defined_name_variant_name = user_defined_name_variant_name
		self.name = "_".join([chr,str(start),str(end),seq])
		self.strand = strand
		self.cut_position = cut_position
		# print (cut_position)
		# exit()
		self.target_pos = target_pos
		self.ref = ref
		self.alt = alt
		self.dist_dict = dist_dict ## other gRNAs in search space
		self.strand_dict = strand_dict ## other gRNAs in search space
		self.candidate_pegRNA_df = candidate_pegRNA_df	 ## other gRNAs in search space
		self.nick_gRNA_list = []
		self.PBS_df = pd.DataFrame()
		self.RTS_df = pd.DataFrame()
		self.rawX = pd.DataFrame()
	def find_PBS(self,min_PBS_length=5,max_PBS_length=18,PBS_length_step=3,**kwargs):
		if len(self.nick_gRNA_list) == 0:
			return False
	
		out = []
		chr = self.chr
		
		if self.strand=="+":
			end = self.cut_position
			for l in range(min_PBS_length,max_PBS_length+1,PBS_length_step):
				start = end-l
				out.append([chr,start,end])
		if self.strand=="-":
			start = self.cut_position - 1
			for l in range(min_PBS_length,max_PBS_length+1,PBS_length_step):
				end = start+l
				out.append([chr,start,end])		
		out = pd.DataFrame(out)	
		out[4] = "PBS"
		self.PBS_df = out.copy()
		# print (self.PBS_df.head())
		# exit()
		# print (self.name,self.PBS_df.shape)
	
	def find_RTS(self,min_RTS_length=10,max_RTS_length=50,RTS_length_step=3,**kwargs):
		if len(self.nick_gRNA_list) == 0:
			# print (self.seq,len(self.nick_gRNA_list))
			return False
		out = []
		chr = self.chr
		if self.strand=="+":
			if len(self.ref)!=len(self.alt):
				target_pos = self.target_pos+1
			else:
				target_pos = self.target_pos
			start = self.cut_position
			for l in range(min_RTS_length,max_RTS_length+1,RTS_length_step):
				end = start+l
				if start+1<=target_pos <=end:
					out.append([chr,start,end])
		if self.strand=="-":
			end = self.cut_position - 1
			for l in range(min_RTS_length,max_RTS_length+1,RTS_length_step):
				start = end-l
				if end>=self.target_pos >=start+1:
					out.append([chr,start,end])		
		out = pd.DataFrame(out)	
		out[4] = "RTS"
		self.RTS_df = out.copy()
		# if out.shape[0] == 0:
			# print (self.chr,self.start,self.end,self.seq,self.strand)
		# print (self.RTS_df.head())
		# exit()
	
	
	def find_nick_gRNA(self,gRNA_search_space=500,**kwargs):
		max_nick_distance = gRNA_search_space
		opposite_strand_gRNAs = [g for g in self.strand_dict if self.strand_dict[g] != self.strand]
		
		for g in opposite_strand_gRNAs:
			# print (g,self.name)
			if self.dist_dict[self.name][g]<=max_nick_distance:
				self.nick_gRNA_list.append(g)	
		if len(self.nick_gRNA_list) == 0:
			print ("no valid nick-gRNAs found for query gRNA %s"%(self.seq))
		# print (self.seq,len(self.nick_gRNA_list))

	
	
	def add_variant(self,genome_fasta=None,N_nick_gRNA_per_sgRNA=5,**kwargs):
		if len(self.nick_gRNA_list) == 0 or self.RTS_df.shape[0]==0:
			# print (self.seq,len(self.nick_gRNA_list))
			return False
	
		## get PBS, RTS sequences
		df = pd.concat([self.PBS_df,self.RTS_df])
		# print (df.shape)
		df[1] = df[1].astype(int)
		df[2] = df[2].astype(int)
		# df = df.dropna()
		# print (df.shape)
		df[3] = self.name
		df[5] = get_opposite_strand(self.strand)
		nick_gRNA = self.candidate_pegRNA_df.loc[self.nick_gRNA_list][[0,1,2,3,4,5,'cut']]
		nick_gRNA[4]="nick-gRNA"
		nick_gRNA[6] = nick_gRNA['cut']-self.target_pos-50
		nick_gRNA[6] = nick_gRNA[6].abs()
		nick_gRNA = nick_gRNA.sort_values(6)
		nick_gRNA = nick_gRNA.head(n=N_nick_gRNA_per_sgRNA)
		nick_gRNA = nick_gRNA.drop(['cut',6],axis=1)
		# N_nick_gRNA_per_sgRNA, pick top N nick gRNA w.r.t 50bp
		
		# print (nick_gRNA)
		df = pd.concat([df,nick_gRNA])
		# if self.seq == "GTCATCTTAGTCATTACCTG":
			# print (df)
		# exit()
		outfile = str(uuid.uuid4()).split("-")[-1]
		df.to_csv(outfile,sep="\t",header=False,index=False)
		# print (df.head())
		df['dummy_index'] = df[3]+"::"+df[0]+":"+df[1].astype(str)+"-"+df[2].astype(str)+"("+df[5]+")"
		seqs = get_seq_from_bed(outfile,genome_fasta)
		seqs[1] = seqs[1].str.upper()
		# print (df)
		# print (seqs.head())
		
		# exit()
		df[3] = seqs.loc[df['dummy_index'].tolist()][1].tolist()
		df = df.drop(['dummy_index'],axis=1)
		# df.to_csv("%s.bed"%(self.name),sep="\t",header=False,index=False)
		# print (df)
		## add variant
		out = []
		# remove_index_list = []
		for line in df.values.tolist():
			RTS = line[3]
			if line[4] != "RTS":
				out.append(line)
				continue
			if line[5] == "-": 
				RTS = revcomp(RTS)
			relative_pos = self.target_pos - int(line[1] ) -1	 ## 0-index
			get_ref = RTS[relative_pos:relative_pos+len(self.ref)]
			if relative_pos == -1:
				## its OK for indels to have -1, because the ref is on PBS sequence, not RTS sequence
				get_ref = self.ref
				# relative_pos = 0
				# get_ref = RTS[relative_pos:relative_pos+len(self.ref)]
				# get_ref = self.ref
				new_RTS = RTS+ self.alt 
			else:
				new_RTS = RTS[:relative_pos]	+ self.alt + RTS[relative_pos+len(self.ref):]
			if get_ref !=self.ref:
				# remove_index_list.append(line[0])
				print ("Something is wrong %s"%(line))
				print ("relative position:",relative_pos)
				print ("RTS:",RTS)
				print ("ref:",self.ref)
				print ("get_ref:",get_ref)	
				print (self.chr,self.start,self.end,self.seq,self.strand)
				if len(self.ref) >1:
					print ("Could be this particular RTS is not long enough to cover the entire reference allel %s"%(self.ref))	
				continue
			
			if line[5] == "-":
				line[3] = revcomp(new_RTS)
			else:
				line[3] = new_RTS
			out.append(line)
		df = pd.DataFrame(out)
		# self.pre_rawX = df.copy()
		self.get_rawX(df)		
		delete_files([outfile])
		
		pass
		
	def get_rawX(self,pegRNA,**kwargs):
	
		# pegRNA = self.pre_rawX
		PBS = pegRNA[pegRNA[4]=="PBS"]
		RTS = pegRNA[pegRNA[4]=="RTS"]
		if RTS.shape[0] == 0:			
			return False
		nick_gRNA = pegRNA[pegRNA[4]=="nick-gRNA"]
		output_index = []
		selected_rows = []
		count = 0
		
		all_list = [PBS.index.tolist(),RTS.index.tolist(),nick_gRNA.index.tolist()]
		# print (all_list)
		temp = list(itertools.product(*all_list)) 
		# print (temp)
		for s in temp:
			count += 1
			current_index = "%s_%s_candidate_%s"%(self.user_defined_name_variant_name,self.name,count)
			for j in s:
				selected_rows.append(j)
				output_index.append(current_index)
			
		out = 	pegRNA.loc[selected_rows]
		out['sample_ID'] = output_index
		out['CHROM'] = self.chr
		out['POS'] = self.target_pos
		out['REF'] = self.ref
		out['ALT'] = self.alt
		out['type'] = out[4]
		out['seq'] = out[3]
		out['chr'] = self.chr
		out['start'] = out[1]
		out['end'] = out[2]
		out['strand'] = out[5]
		out = out[["sample_ID","CHROM","POS","REF","ALT","type","seq","chr","start","end","strand"]]

		# gRNA.columns = [str(x) for x in gRNA.columns]

		out2 = pd.DataFrame()
		out2['sample_ID'] = list(set(output_index))
		out2['CHROM'] = self.chr
		out2['POS'] = self.target_pos
		out2['REF'] = self.ref
		out2['ALT'] = self.alt
		out2['type'] = "pegRNA"

		out2['seq'] = self.seq
		out2['chr'] = self.chr
		out2['start'] = self.start
		out2['end'] = self.end
		out2['strand'] = self.strand

		out = pd.concat([out,out2])
		out = out.sort_values("sample_ID")
		# print (out.head())
		# out.to_csv("%s.bed"%(self.name),sep="\t",header=False,index=False)
		self.rawX = out.copy()	
	
	
	

# from .utils import *
from utils import *
	
class sgRNA:
	def __init__(self,chr=None,start=None,end=None,seq=None,strand=None,cut_position=None,mutation_pos = None,mutation_ref = None,mutation_alt = None,variant_id=None,dist_dict=None,opposite_strand_sgRNAs=None,all_sgRNA_df=None,**kwargs):
		"""
		
		search key is chr_start_end_seq because seq can be duplicated, this is unique
		
		
		Searching steps
		--------------------
		
		1. RTT
		2. PBS
		3. ngRNA
		
		
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
		self.target_pos = target_pos ## this is the corrected, actual position of variant
		self.max_nick_distance = max_nick_distance
		self.ref = ref
		self.alt = alt
		self.dist_dict = dist_dict ## other gRNAs in search space
		self.strand_dict = strand_dict ## other gRNAs in search space
		self.candidate_pegRNA_df = candidate_pegRNA_df	 ## other gRNAs in search space
		self.nick_gRNA_list = []
		
		#for rawX
		self.PBS_df = pd.DataFrame()
		self.RTT_df = pd.DataFrame()
		self.ngRNA_df = pd.DataFrame()
		# for X
		self.PBS_features = pd.DataFrame() 
		self.RTT_features = pd.DataFrame()
		self.ngRNA_features = pd.DataFrame()
		
		
		self.rawX = pd.DataFrame()
		# sgRNA name = chr,start,end,seq,strand
		
		## flags
		self.no_RTT = False ## alwasy assume we can find valid RTT
		
		
		
		
		
	def find_PBS(self,min_PBS_length=7,max_PBS_length=17,**kwargs):
		"""find PBS sequences as bed format df
		
		Output
		--------
		
		PBS coordinates, as a dataframe

		chr,start,end, strand 
		
		"""
		if self.no_RTT:
			return 0
	
		out = []
		chr = self.chr
		
		if self.strand=="+":
			end = self.cut_position
			for l in range(min_PBS_length,max_PBS_length+1):
				start = end-l
				out.append([chr,start,end])
		if self.strand=="-":
			start = self.cut_position - 1
			for l in range(min_PBS_length,max_PBS_length+1):
				end = start+l
				out.append([chr,start,end])		
		self.PBS_df = pd.DataFrame(out)	
		self.PBS_df[3] = get_opposite_strand(self.strand)
		
		
		## make features
		
		

	def find_RTT(self,min_RTT_length=7,max_RTT_length=40,min_distance_RTT5=5,**kwargs):
		"""find RTT sequences as bed format df
		
		Output
		--------
		
		RTT coordinates, as a dataframe

		chr,start,end, strand 
		
		
		
		"""
		out = []
		chr = self.chr
		if self.strand=="+":
			start = self.cut_position # remember out cut position, the actual nucleotide, we use -4
			for l in range(min_RTT_length,max_RTT_length+1):
				end = start+l
				if start+1<=self.target_pos <=end-min_distance_RTT5:
					out.append([chr,start,end])
		if self.strand=="-":
			end = self.cut_position - 1
			for l in range(min_RTT_length,max_RTT_length+1):
				start = end-l
				if end>=self.target_pos >=start+1+min_distance_RTT5:
					out.append([chr,start,end])	
		if len(out) == 0:
			## valid RTT not found
			self.no_RTT = True
			return 0
		self.RTT_df = pd.DataFrame(out)
		self.RTT_df.columns = ['chr','start','end']
		self.RTT_df["strand"] = get_opposite_strand(self.strand)
		
		temp = get_fasta_simple(self.target_fa,self.RTT_df, self.target_pos)
		self.RTT_df['old_seq'] = temp[3].tolist()
		## add variant
		self.RTT_df['seq'] = [add_variant(r['old_seq'],self.mutation_pos-r['start'],self.mutation_ref,self.mutation_alt) for i,r in self.RTT_df.iterrows()]
		
		## make features
		feature_list = []
		


	def find_nick_gRNA(self,max_nick_distance=150,debug=0,**kwargs):
		"""find all valid ngRNAs given the sgRNA name
		
		Input
		------
		
		all possible ngRNA list is a list of ngRNA names that are on the opposite strand of the given gRNA
		
		Then we just need to check the distance threshold
		
		
		"""

		for g in self.ngRNA_list:
			if abs(self.dist_dict[self.name][g])<=max_nick_distance:
				self.nick_gRNA_list.append(g)
		if debug > 5:
			if len(self.nick_gRNA_list) == 0:
				print ("no valid nick-gRNAs found for query gRNA %s"%(self.seq))
		# print (self.seq,len(self.nick_gRNA_list))
		
		
		## make features

	
	
	def add_variant(self,RTT_seq,relative_pos,ref,alt,**kwargs):
		"""
		
		Steps
		-------
		
		Assume all input is on the positive strand
		
		ref alt is the corrected version
		
		e.g., GC - > C becomes C to "" (empty)
		
		relative_pos 0-index
		
		"""
		

		
		## function check
		if len(ref)>0:
			get_ref = RTT_seq[relative_pos:relative_pos+len(ref)]
			if get_ref != ref:
				print (self.name,"references do not match",RTT_seq,relative_pos,ref,alt)
		if relative_pos < 0:
			print (self.name,"relative position < 0, result will not be correct!")
			relative_pos = 0
		new_RTT = RTT_seq[:relative_pos]	+ self.alt + RTT_seq[relative_pos+len(self.ref):]
		return new_RTT

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
		out2['type'] = "sgRNA"

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
	
	
	
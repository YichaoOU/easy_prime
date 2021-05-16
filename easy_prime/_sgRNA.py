
from .utils import *
# from utils import *
	
class sgRNA:
	def __init__(self,chr=None,start=None,end=None,seq=None,sgRNA_name=None,strand=None,cut_position=None,mutation_pos = None,mutation_ref = None,mutation_alt = None,variant_id=None,dist_dict=None,opposite_strand_sgRNAs=None,all_sgRNA_df=None,target_fa=None,scaffold_seq=None,user_target_pos=None,user_ref=None,user_alt=None,offset=None,target_to_sgRNA=None,PAM=None,DeepSpCas9=None,**kwargs):
	
		"""
		
		search key is chr_start_end_seq because seq can be duplicated, this chr_start_end_seq is unique. This name doesn't show strand, but seq is alwasy 5-3 direction
		
		
		Searching steps
		--------------------
		
		1. RTT
		2. PBS
		3. ngRNA
		
							   chr     start       end  target_to_RTT5 strand  ...        11        12        13    RTT_GC  RTT_length
		01d0f4a564dc_RTT_0   chr20  55990403  55990410               5      -  ...  0.224020  0.256546  0.297711  0.857143           7
		01d0f4a564dc_RTT_33  chr20  55990403  55990443              38      -  ...  0.435531  0.151318  0.172131  0.675000          40

		[34 rows x 24 columns]
							   chr     start       end strand                 seq    PBS_GC  PBS_length
		01d0f4a564dc_PBS_0   chr20  55990396  55990403      -             ACTTATC  0.285714           7
		01d0f4a564dc_PBS_11  chr20  55990385  55990403      -  ACTTATCGGCCCCTGCGG  0.666667          18
								chr     start       end                   seq                                      sgRNA_name strand  sgRNA_distance_to_ngRNA
		01d0f4a564dc_ngRNA_0  chr20  55990397  55990417  AAGCAGAGCCAGACACTTAT  chr20_55990397_55990417_-_AAGCAGAGCCAGACACTTAT      -                       -2
		01d0f4a564dc_ngRNA_1  chr20  55990388  55990408  CAGACACTTATCGGCCCCTG  chr20_55990388_55990408_-_CAGACACTTATCGGCCCCTG      -                      -11
		01d0f4a564dc_ngRNA_2  chr20  55990387  55990407  AGACACTTATCGGCCCCTGC  chr20_55990387_55990407_-_AGACACTTATCGGCCCCTGC      -                      -12
		01d0f4a564dc_ngRNA_3  chr20  55990384  55990404  CACTTATCGGCCCCTGCGGG  chr20_55990384_55990404_-_CACTTATCGGCCCCTGCGGG      -                      -15
		01d0f4a564dc_ngRNA_4  chr20  55990378  55990398  TCGGCCCCTGCGGGAGGCCC  chr20_55990378_55990398_-_TCGGCCCCTGCGGGAGGCCC      -                      -21
				
		
		"""
		self.chr = chr
		self.start = start
		self.end = end
		self.seq = seq
		self.variant_id = variant_id
		self.target_to_sgRNA = target_to_sgRNA
		self.DeepSpCas9 = DeepSpCas9
		self.is_dPAM = 0
		self.PAM = PAM
		self.offset = offset
		self.sgRNA_name = sgRNA_name
		# print ("init",sgRNA_name)
		self.uid = str(uuid.uuid4()).split("-")[-1]
		# self.name = "_".join([chr,str(start),str(end),seq])
		self.strand = strand
		self.cut_position = cut_position
		# print (cut_position)
		# exit()
		self.target_pos = mutation_pos ## this is the corrected, actual position of variant
		self.target_fa = target_fa 
		self.scaffold_seq = scaffold_seq 
		# self.max_nick_distance = max_nick_distance
		self.ref = mutation_ref
		self.alt = mutation_alt
		self.dist_dict = dist_dict ## other gRNAs in search space
		self.opposite_strand_sgRNAs = opposite_strand_sgRNAs 
		# self.candidate_pegRNA_df = candidate_pegRNA_df	 ## other gRNAs in search space
		self.nick_gRNA_list = []
		self.user_target_pos = user_target_pos
		self.user_ref = user_ref
		self.user_alt = user_alt

		#for rawX
		self.PBS_df = pd.DataFrame()
		self.RTT_df = pd.DataFrame()
		self.ngRNA_df = pd.DataFrame()
		# for X
		# self.PBS_feature_list = []
		# self.RTT_feature_list = []
		# self.ngRNA_feature_list = []
		self.PBS_feature_list = ["PBS_GC",'PBS_length']
		# self.RTT_feature_list = ["target_to_RTT5","RTT_GC","RTT_length","0","1","2","3","4","5","6","7","8","9","10","11","12","13"]
		self.RTT_feature_list = ["target_to_RTT5","RTT_GC","RTT_length","0","1","2","3","4","5","6","7","8","9"]
		self.ngRNA_feature_list = ['sgRNA_distance_to_ngRNA','is_PE3b']
		
		self.rawX = pd.DataFrame()
		self.X = pd.DataFrame()
		self.X_p = pd.DataFrame()
		# sgRNA name = chr,start,end,seq,strand
		
		## flags
		self.no_RTT = False ## alwasy assume we can find valid RTT
		self.no_ngRNA = False ## alwasy assume we can find valid ngRNA
		
		
	def find_PBS(self,min_PBS_length=7,max_PBS_length=17,**kwargs):
		"""find PBS sequences as bed format df
		
		Output
		--------
		
		PBS_df
		
		PBS_feature_list
		
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
		self.PBS_df.columns = ['chr','start','end']
		self.PBS_df["strand"] = get_opposite_strand(self.strand)
		self.PBS_df.index = ["%s_PBS_%s"%(self.uid,i) for i in range(self.PBS_df.shape[0])]
		temp = get_fasta_simple(self.target_fa,self.PBS_df, self.user_target_pos)
		self.PBS_df['seq'] = temp[3].tolist()
		
		if self.strand == "+": ## when sgRNA is positive strand, RTT should use the negative strand
			self.PBS_df['seq'] = [revcomp(x) for x in self.PBS_df['seq']]

		self.PBS_df['PBS_GC'] = [GC_content(x) for x in self.PBS_df['seq'] ]
		self.PBS_df['PBS_length'] = [len(x) for x in self.PBS_df['seq'] ]
		

	def find_RTT(self,debug=0,min_RTT_length=10,max_RTT_length=20,max_max_RTT_length=40,min_distance_RTT5=5,**kwargs):
		"""find RTT sequences as bed format df
		
		Output
		--------
		
		RTT_df
		
		RTT_feature_list
		
		"""
		out = []
		chr = self.chr
		pbs_start = None
		pbs_end = None
		user_max_RTT_length = max_RTT_length
		large_deletion_flag=False
		if len(self.ref)+min_distance_RTT5 > max_RTT_length: ## in case of large deletion
			large_deletion_flag = True
			max_RTT_length = max_RTT_length+len(self.ref)
		# target_to_RTT5_feature=[]
		deletion_length = len(self.ref)-len(self.alt)
		if self.strand=="+":
			start = self.cut_position # remember out cut position, the actual nucleotide, we use -4
			pbs_end = start
			pbs_start = pbs_end - 14 
			for l in range(min_RTT_length,max_RTT_length+1+max(0,deletion_length)):
				end = start+l
				if start+1<=self.target_pos <=end-min_distance_RTT5:
					out.append([chr,start,end,end-self.target_pos-max(0,deletion_length)])
		if self.strand=="-":
			end = self.cut_position - 1
			pbs_start = end
			pbs_end = pbs_start + 14 
			for l in range(min_RTT_length,max_RTT_length+1+max(0,deletion_length)):
				start = end-l
				if end>=self.target_pos >=start+1+min_distance_RTT5:
					out.append([chr,start,end,self.target_pos-start-1-max(0,deletion_length)])	
					
		current_max_RTT_length = max_RTT_length+5
		while len(out) == 0:
			if current_max_RTT_length > max_max_RTT_length:
				break
			out,_,no_RTT_flag = self.find_longer_RTT(min_RTT_length=min_RTT_length,max_RTT_length=current_max_RTT_length,min_distance_RTT5=min_distance_RTT5,**kwargs)
			if len(out) > 0:
				print ("max_RTT_length increased from %s to %s"%(max_RTT_length,current_max_RTT_length))
				break
			current_max_RTT_length += 5
					
		if len(out) == 0:
			## valid RTT not found
			self.no_RTT = True
			print ("No valid RTT found given current max length:%s"%(max_max_RTT_length))
			return 0
		self.RTT_df = pd.DataFrame(out)
		self.RTT_df.columns = ['chr','start','end','target_to_RTT5']
		self.RTT_df["strand"] = get_opposite_strand(self.strand)
		self.RTT_df.index = ["%s_RTT_%s"%(self.uid,i) for i in range(self.RTT_df.shape[0])]
		if debug>10:
			print (self.target_fa)
			print (self.user_target_pos)
		temp = get_fasta_simple(self.target_fa,self.RTT_df, self.user_target_pos)
		self.RTT_df['old_seq'] = temp[3].tolist()
	
		## add variant
		# relative_pos = self.target_pos-r['start']-1 # start is 0-index
		if debug>10:
			pd.set_option('display.max_columns', None)
			print (self.RTT_df)
		self.RTT_df['seq'] = [self.add_variant(r['old_seq'],self.target_pos-r['start']-1,self.ref,self.alt) for i,r in self.RTT_df.iterrows()]
		self.RTT_df = self.RTT_df[self.RTT_df['seq']!=0]
		if self.strand == "+": ## when sgRNA is positive strand, RTT should use the negative strand
			self.RTT_df['seq'] = [revcomp(x) for x in self.RTT_df['seq']]
		# print (self.RTT_df)
		self.RTT_df['RTT_length'] = [len(x) for x in self.RTT_df['seq'] ]
		if large_deletion_flag:
			self.RTT_df = self.RTT_df[self.RTT_df['RTT_length']<=user_max_RTT_length]
		self.RTT_df = self.RTT_df[self.RTT_df['RTT_length']>=min_RTT_length]
		# filter out first C
		self.RTT_df['firstC'] = [x[0] for x in self.RTT_df['seq']]
		self.RTT_df = self.RTT_df[self.RTT_df['firstC']!="C"]
		self.RTT_df = self.RTT_df.drop(['firstC'],axis=1)	

		
		current_max_RTT_length = max_RTT_length+5
		while self.RTT_df.shape[0] == 0:
			if debug>=10:
				print ("increasing max_RTT_length to:", current_max_RTT_length)
			if current_max_RTT_length > max_max_RTT_length:
				break
			_,RTT_df,no_RTT_flag = self.find_longer_RTT(min_RTT_length=min_RTT_length,max_RTT_length=current_max_RTT_length,min_distance_RTT5=min_distance_RTT5,**kwargs)
			if RTT_df.shape[0] > 0:
				print ("max_RTT_length increased from %s to %s"%(max_RTT_length,current_max_RTT_length))
				self.RTT_df = RTT_df
				break
			current_max_RTT_length += 5		

		
		
		if self.RTT_df.shape[0] == 0:
			## valid RTT not found
			self.no_RTT = True
			print ("NO RTT found: new RTT sequence length is longer than: %s"%(max_max_RTT_length))
			return 0

		
		## make features
		## modify it to fit fasta input
		attached_minimal_PBS = sub_fasta_single(self.target_fa,self.user_target_pos, pbs_start,pbs_end)
		if self.strand == "+": ## when sgRNA is positive strand, RTT should use the negative strand
			attached_minimal_PBS = revcomp(attached_minimal_PBS)
		# print ("attached_minimal_PBS",attached_minimal_PBS)
		# print (self.RTT_df)

		
		self.RTT_df['RNAfold_seq']= [(self.scaffold_seq+x+attached_minimal_PBS).replace("T","U") for x in self.RTT_df['seq']]
		RNAfold_features_df = pd.DataFrame([call_RNAplfold(x,len(self.scaffold_seq))for x in self.RTT_df['RNAfold_seq']])
		# print (RNAfold_features_df)
		RNAfold_features_df.index = self.RTT_df.index.tolist()
		self.RTT_df = pd.concat([self.RTT_df,RNAfold_features_df],axis=1)
		self.RTT_df['RTT_GC'] = [GC_content(x) for x in self.RTT_df['seq'] ]
		
		self.RTT_df.columns = [str(x) for x in self.RTT_df.columns]
		self.is_dPAM = is_dPAM(self.PAM, self.RTT_df['seq'][0], self.offset)
		# print (self.RTT_df)
		
	def find_longer_RTT(self,min_RTT_length=10,max_RTT_length=20,min_distance_RTT5=5,**kwargs):
		"""find RTT sequences as bed format df
		
		Output
		--------
		
		RTT_df
		
		RTT_feature_list
		
		"""
		out = []
		chr = self.chr
		pbs_start = None
		pbs_end = None
		user_max_RTT_length = max_RTT_length
		large_deletion_flag=False
		if len(self.ref)+min_distance_RTT5 > max_RTT_length: ## in case of large deletion
			large_deletion_flag = True
			max_RTT_length = max_RTT_length+len(self.ref)
		# target_to_RTT5_feature=[]
		if self.strand=="+":
			start = self.cut_position # remember out cut position, the actual nucleotide, we use -4
			pbs_end = start
			pbs_start = pbs_end - 14 
			for l in range(min_RTT_length,max_RTT_length+1):
				end = start+l
				if start+1<=self.target_pos <=end-min_distance_RTT5:
					out.append([chr,start,end,end-self.target_pos])
		if self.strand=="-":
			end = self.cut_position - 1
			pbs_start = end
			pbs_end = pbs_start + 14 
			for l in range(min_RTT_length,max_RTT_length+1):
				start = end-l
				if end>=self.target_pos >=start+1+min_distance_RTT5:
					out.append([chr,start,end,self.target_pos-start-1])	
		if len(out) == 0:
			return out,self.RTT_df,True
		self.RTT_df = pd.DataFrame(out)
		self.RTT_df.columns = ['chr','start','end','target_to_RTT5']
		self.RTT_df["strand"] = get_opposite_strand(self.strand)
		self.RTT_df.index = ["%s_RTT_%s"%(self.uid,i) for i in range(self.RTT_df.shape[0])]
		temp = get_fasta_simple(self.target_fa,self.RTT_df, self.user_target_pos)
		self.RTT_df['old_seq'] = temp[3].tolist()
		
		## add variant
		# relative_pos = self.target_pos-r['start']-1 # start is 0-index
		self.RTT_df['seq'] = [self.add_variant(r['old_seq'],self.target_pos-r['start']-1,self.ref,self.alt) for i,r in self.RTT_df.iterrows()]
		self.RTT_df = self.RTT_df[self.RTT_df['seq']!=0]
		if self.strand == "+": ## when sgRNA is positive strand, RTT should use the negative strand
			self.RTT_df['seq'] = [revcomp(x) for x in self.RTT_df['seq']]
		# print (self.RTT_df)
		self.RTT_df['RTT_length'] = [len(x) for x in self.RTT_df['seq'] ]
		if large_deletion_flag:
			self.RTT_df = self.RTT_df[self.RTT_df['RTT_length']<=user_max_RTT_length]
		self.RTT_df = self.RTT_df[self.RTT_df['RTT_length']>=min_RTT_length]
		# filter out first C
		self.RTT_df['firstC'] = [x[0] for x in self.RTT_df['seq']]
		self.RTT_df = self.RTT_df[self.RTT_df['firstC']!="C"]
		self.RTT_df = self.RTT_df.drop(['firstC'],axis=1)	
		if self.RTT_df.shape[0] == 0:
			return out,self.RTT_df,True
		return out,self.RTT_df,False
			

	def make_PE3b_ngRNA(self,myIndex,start,end,seq,strand):
		"""check if target mutation overlaps with ngRNA and update its sequence
		
		only work for len(ref) = len(alt), namely substitutions
		
		return 
		
		sequence
		
		"""

		# check if overlaps
		
		target_pos_list = list(range(self.target_pos,self.target_pos+len(self.ref)))
		seq_pos_list = list(range(start+1,end+1))
		overlaps = set(seq_pos_list).intersection(target_pos_list)
		
		if len(overlaps) == 0:
			return [myIndex,seq,0]
		# print (seq,"contain overlaps",overlaps)
		if strand == "-":
			new_seq = list(revcomp(seq))
		else:
			new_seq = list(seq)
		for i in overlaps:
			relative_ngRNA_pos = i - start - 1
			relative_target_pos = i - self.target_pos
			if relative_ngRNA_pos < 0 or relative_target_pos<0:
				print ("Error: make_PE3b_ngRNA",relative_target_pos,relative_ngRNA_pos,self.sgRNA_name,start,end,new_seq,strand,self.variant_id)
				return [myIndex,seq,0]
				relative_ngRNA_pos = 0
			ngRNA_ref = new_seq[relative_ngRNA_pos]
			target_ref = self.ref[relative_target_pos]
			target_alt = self.alt[relative_target_pos]
			if target_ref != ngRNA_ref:
				print ("Error: target_ref != ngRNA_ref","%s != %s"%(target_ref,ngRNA_ref),relative_target_pos,relative_ngRNA_pos,start,end,new_seq,strand,self.variant_id)
				return [myIndex,seq,0]
			new_seq[relative_ngRNA_pos] = target_alt
		if strand == "-":
			new_seq = revcomp("".join(new_seq))
			return [myIndex,new_seq,1]
		else:
			new_seq = "".join(new_seq)
			return [myIndex,new_seq,1]
		
	def find_nick_gRNA(self,max_ngRNA_distance=100,max_max_ngRNA_distance=200,debug=0,**kwargs):
		"""find all valid ngRNAs given the sgRNA name
		
		Input
		------
		
		ngRNA_df
		
		ngRNA_feature_list
		
		"""

		self.ngRNA_df = self.opposite_strand_sgRNAs.copy()
		# print (self.ngRNA_df)
		if self.ngRNA_df.shape[0] == 0:
			self.no_ngRNA = True
			print ("no ngRNA for %s"%(self.sgRNA_name))
			return 0
		self.ngRNA_df.index = ["%s_ngRNA_%s"%(self.uid,i) for i in range(self.ngRNA_df.shape[0])]
		self.ngRNA_df.columns = ['chr','start','end','seq','sgRNA_name','strand']
		if self.strand == "-":
			self.ngRNA_df['sgRNA_distance_to_ngRNA'] = [-self.dist_dict[x][self.sgRNA_name] for x in self.ngRNA_df["sgRNA_name"]]
		if self.strand == "+":
			self.ngRNA_df['sgRNA_distance_to_ngRNA'] = [self.dist_dict[x][self.sgRNA_name] for x in self.ngRNA_df['sgRNA_name']]
		# print (max_ngRNA_distance,"max_ngRNA_distance")
		current_ngRNA_df = self.ngRNA_df[self.ngRNA_df['sgRNA_distance_to_ngRNA'].abs()<=max_ngRNA_distance]
		current_max_ngRNA_distance = max_ngRNA_distance + 20
		while current_ngRNA_df.shape[0] == 0:
			if current_max_ngRNA_distance > max_max_ngRNA_distance:
				self.ngRNA_df = current_ngRNA_df
				break
			current_ngRNA_df = self.ngRNA_df[self.ngRNA_df['sgRNA_distance_to_ngRNA'].abs()<=current_max_ngRNA_distance]
			if current_ngRNA_df.shape[0] > 0:
				self.ngRNA_df = current_ngRNA_df
				print ("max ngRNA distance increased from %s to %s:"%(max_ngRNA_distance,current_max_ngRNA_distance))
				break
			current_max_ngRNA_distance += 20
		self.ngRNA_df = current_ngRNA_df
		if self.ngRNA_df.shape[0] == 0:
			print ("no ngRNA for max distance: %s"%(max_ngRNA_distance))
			self.no_ngRNA = True
			return 0
		
		self.ngRNA_df['is_PE3b'] = 0
		
		if len(self.ref) == len(self.alt): # check PE3b
			pe3b = pd.DataFrame([self.make_PE3b_ngRNA(		i,
																														r['start'],
																														r['end'],
																														r['seq'],
																														r['strand']
																													) for i,r in self.ngRNA_df.iterrows()])
			# print ("-"*20,self.variant_id,"-"*20)
			# print (pe3b)
			pe3b = pe3b.set_index(0)
			self.ngRNA_df[['seq','is_PE3b']] = pe3b[[1,2]]
		
		# print (self.ngRNA_df)

	
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
			# print (self.sgRNA_name,"relative_pos",relative_pos,"get_ref",get_ref)
			if get_ref != ref:
				if len(get_ref)<len(ref):
					## large deletion, given RTT is not long enough
					return 0
				print (self.variant_id,self.sgRNA_name,"references do not match",RTT_seq,relative_pos,ref,alt)
				return 0
		if relative_pos < 0:
			print (self.variant_id,self.sgRNA_name,"relative position < 0, result will not be correct!")
			return 0
			relative_pos = 0
		new_RTT = RTT_seq[:relative_pos]	+ self.alt + RTT_seq[relative_pos+len(self.ref):]
		return new_RTT

	def get_rawX_and_X(self,debug=0,**kwargs):
		"""get rawX and X formated dataframe
		
		given PBS, RTT, ngRNA sequences and features
		
		we can use the itertools to take all combinations of the index
		
		then we concat them as rows, -> rawX
		
		we concat them as columns -> X
		
		return
		------
		
		self.rawX
		
		self.X

		
		"""
		
		rawX_columns = ["seq","chr","start","end","strand"]
		
		if self.no_RTT:
			return 0
		## not allow PE2 cases
		if self.no_ngRNA:
			print ("no_ngRNA found",self.sgRNA_name)
			return 0
			self.ngRNA_df=pd.DataFrame([np.nan]*(len(self.ngRNA_feature_list)+len(rawX_columns))).T
			self.ngRNA_df.columns = rawX_columns + self.ngRNA_feature_list
		if debug>=10:
			print (self.variant_id,"showing PBS_df, RTT_df, and ngRNA df")
			print (self.PBS_df.head())
			print (self.RTT_df.head())
			print (self.ngRNA_df.head())
		X_index = []
		PBS_selected_rows =[]
		RTT_selected_rows =[]
		ngRNA_selected_rows =[]
		all_list = [self.PBS_df.index.tolist(),self.RTT_df.index.tolist(),self.ngRNA_df.index.tolist()]
		temp = list(itertools.product(*all_list)) 
		count = 0
		for s in temp:
			count += 1
			current_index = "%s_%s_candidate_%s"%(self.variant_id,self.sgRNA_name,count)
			# current_index = "%s_%s"%(self.variant_id,count) # shorter name
			if count == 1  and debug>10:
				print (current_index)
			PBS_selected_rows.append(s[0])
			RTT_selected_rows.append(s[1])
			ngRNA_selected_rows.append(s[2])
			if self.ngRNA_df.at[s[2],'is_PE3b']==1:
				current_index+="_PE3b"
				# print (s[2],"PE3B")
			if self.is_dPAM:
				current_index+="_dPAM"	
				# print (self.sgRNA_name,"dPAM")		
			if self.no_ngRNA:
				current_index+="_PE2"	
			X_index.append(current_index)
			
			
		# ------------------------------------------------------------  rawX  ------------------------------------------------------------
		rawX_PBS = self.PBS_df.loc[PBS_selected_rows][rawX_columns]
		rawX_RTT = self.RTT_df.loc[RTT_selected_rows][rawX_columns]
		rawX_ngRNA = self.ngRNA_df.loc[ngRNA_selected_rows][rawX_columns]
		rawX_sgRNA = pd.DataFrame([self.seq,self.chr,self.start,self.end,self.strand]).T
		rawX_sgRNA.columns = rawX_columns
		rawX_sgRNA = rawX_sgRNA.loc[[0]*len(ngRNA_selected_rows)][rawX_columns]
		rawX_PBS['sample_ID']  = X_index
		rawX_RTT['sample_ID']  = X_index
		rawX_ngRNA['sample_ID']  = X_index
		rawX_sgRNA['sample_ID']  = X_index
		rawX = pd.concat([rawX_PBS,rawX_RTT,rawX_ngRNA,rawX_sgRNA])
		rawX['CHROM'] = self.chr
		rawX['POS'] = self.user_target_pos
		rawX['REF'] = self.user_ref
		rawX['ALT'] = self.user_alt
		rawX['type'] = ['PBS']*rawX_PBS.shape[0]+['RTT']*rawX_RTT.shape[0]+['ngRNA']*rawX_ngRNA.shape[0]+['sgRNA']*rawX_ngRNA.shape[0]
		self.rawX = rawX[["sample_ID","CHROM","POS","REF","ALT","type","seq","chr","start","end","strand"]]
		self.rawX = self.rawX.sort_values("sample_ID")
		self.rawX.index = self.rawX['sample_ID'].tolist()
		


		# ------------------------------------------------------------  X  ------------------------------------------------------------
		
		X_PBS = self.PBS_df.loc[PBS_selected_rows][self.PBS_feature_list]
		X_PBS = X_PBS.reset_index(drop=True)
		X_RTT = self.RTT_df.loc[RTT_selected_rows][self.RTT_feature_list]
		X_RTT = X_RTT.reset_index(drop=True)
		X_ngRNA = self.ngRNA_df.loc[ngRNA_selected_rows][self.ngRNA_feature_list]
		X_ngRNA = X_ngRNA.reset_index(drop=True)
		self.X = pd.concat([X_PBS,X_RTT,X_ngRNA],axis=1)
		self.X.index = X_index	
		self.X['is_dPAM'] = self.is_dPAM
		self.X['target_to_sgRNA'] = self.target_to_sgRNA
		self.X['DeepSpCas9'] = self.DeepSpCas9


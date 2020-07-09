
from .utils import *
from ._sgRNA import sgRNA

# import _sgRNA
# import utils
# import imp
# imp.reload(utils)
# imp.reload(_sgRNA)
# from utils import *
# from _sgRNA import sgRNA

"""
'nick_to_pegRNA', 'target_to_pegRNA', 'target_to_RTT5','aln_ref_alt_mis', 'aln_ref_alt_del', 'aln_ref_alt_ins',  'PBS_GC', 'RTS_GC','PBS_length', 'RTS_length', 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 'is_dPAM'


nick_to_pegRNA,target_to_pegRNA,'is_dPAM': defined here
'aln_ref_alt_mis', 'aln_ref_alt_del', 'aln_ref_alt_ins': defined here


target_to_RTT5, RTS_GC, RTS_length,0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 : defined in sgRNA class find_RTT

PBS_GC,PBS_length:  defined in sgRNA class find_PBS
"""


def find_mutation_pos(pos,ref,alt):
	"""solve the problem when user specification of ref, alt contain redundancy
	ATTTT-> ATTT, should be T -> ""
	
	GC-> GGC should be G -> GC
	
	G - > GC will just be G GC
	
	ref can't be empty, otherwise we can't check if add variant is correct or not
	"""
	count=0
	for i in range(min(len(ref),len(alt))):
		# if len(ref)-i <= 1:
			# return pos,ref[i:],alt[i:]
		x=ref[i]
		y=alt[i]
		if x != y:
			return pos,ref[i:],alt[i:]
		else:
			pos+=1
			count+=1
	return pos,ref[count:],alt[count:]


class target_mutation:
	def __init__(self,chr,pos,name,ref,alt,target_fa,**kwargs):
		"""target mutation class
		
		
		sgRNA name: chr_start_end_strand_seq
		
		target_mutation name: id_chr_pos_ref_alt
		
		pos is corrected, and the corrected pos, ref, alt is used 
		
		target_fa is the +-1000 extended sequences
		
		"""
		self.chr = chr
		self.target_pos = pos
		self.name = name.replace("/","_").replace(",","_")
		self.ref = ref
		self.alt = alt
		self.target_fa = target_fa
		self.debug_folder = "easy_prime_debug_files"
		self.dist_dict = {}
		self.strand_dict = {}
		self.rawX = pd.DataFrame()
		self.X = pd.DataFrame()
		self.X_p = pd.DataFrame()
		self.topX = pd.DataFrame()
		self.allX = pd.DataFrame()
		self.pegRNA_flag=True
		## flags
		self.found_PE3b = False
		self.found_PE3 = False
		self.found_PE2 = False
		self.found_dPAM = False
		self.N_sgRNA_found = 0
		

		self.feature_for_prediction = ["sgRNA_distance_to_ngRNA","target_to_sgRNA","target_to_RTT5","N_subsitution","N_deletion","N_insertions","PBS_GC","RTT_GC","PBS_length","RTT_length",'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13',"is_dPAM"] # match the order of training features
		
		
		
		self.mutation_pos,self.mutation_ref,self.mutation_alt = find_mutation_pos(pos,ref,alt)
		# print (self.target_pos,self.ref,self.alt)
		# print ("%s minimal mutation format: "%(self.name),self.mutation_pos,self.mutation_ref,self.mutation_alt)
		
		## helper structures
		self.sgRNA_strand_df={}
		self.sgRNA_strand_df["+"]=pd.DataFrame()
		self.sgRNA_strand_df["-"]=pd.DataFrame()
		self.valid_init_sgRNA = pd.DataFrame()
		self.all_sgRNA = pd.DataFrame() # used to find ngRNA
		
		#---------------- features ------------
		# target mutation feature
		self.ref_alt = global_alignments(self.ref,self.alt)
		# >>> t.ref_alt
		# [0, 1, 0, 0]
		# sgRNA distance to target mutation
		self.sgRNA_target_distance_dict = {} ## contain valid and invalid sgRNA in the key, but the latter with  distance <0
		self.sgRNA_target_dPAM_dict = {} ## contain binary values of whether the target affect this sgRNA PAM
		# sgRNA distance to ngRNA 
		self.dist_dict = {} ## sgRNA_ngRNA_distance_dict

		
	def init(self,gRNA_search_space=200,search_iteration=1,sgRNA_length=20,PAM="NGG",offset=-3,debug=0,genome_fasta=None,max_RTT_length=40,min_distance_RTT5=5,max_target_to_sgRNA=10,**kwargs):
		"""First step: search for candidate sgRNAs around target mutation
		
		Input
		-----
		
		gRNA_search_space: extend pos by +- gRNA_search_space
		
		search_iteration: if in the search space defined by gRNA_search_space, we fail to find sgRNAs, we will extend the gRNA_search_space further to find at least one sgRNA. (no need to increase it)
		
		Output
		------
		
		we have the following columns
		
		chr, start, end, sgRNA name, seq, strand, cut_position, valid
		
		self.dist_dict  2d dict
		
		self.sgRNA_target_distance_dict
		
		self.sgRNA_strand_df={}
		self.sgRNA_strand_df["+"]=None
		self.sgRNA_strand_df["-"]=None
		self.valid_init_sgRNA = None
		self.all_sgRNA = None # used to find ngRNA

		These will be used later
		------
		
		self.offset
		self.PAM
		
		"""
		## debug folder
		if debug>0:
			subprocess.call("mkdir -p %s"%(self.debug_folder),shell=True)
		## setup
		self.offset = offset
		self.PAM = PAM
		## other sequences, coordinates are all computed from this sequence

		### find all sgRNA given a sequence
		for i in range(search_iteration):
			extend = gRNA_search_space*(i+1)
			if i >=1:
				print ("No sgRNA were found using %s gRNA_search_space"%(extend))
			start = self.mutation_pos-extend
			end = self.mutation_pos+extend
			search_fa = sub_fasta_single(self.target_fa,self.target_pos, start,end)
			df = run_pam_finder(search_fa,"N"*sgRNA_length,self.PAM,start,self.chr)
			## this above df contains all sgRNAs
			self.N_sgRNA_found = df.shape[0]
			
			try:
				df[1] = df[1].astype(int)
				df[2] = df[2].astype(int)
				## sgRNA name
				df[4] = df[0]+"_"+df[1].astype(str)+"_"+df[2].astype(str)+"_"+df[5].astype(str)+"_"+df[3].astype(str)
				df.index = df[4].to_list()
				df['cut'] = [get_gRNA_cut_site(x[1],x[2],x[5],self.offset) for i,x in df.iterrows()]
				## distance is the coordiantes used in D Liu paper
				# ooo|abcNGG
				# -3-2-1|123NGG
				# in the above example, b has a distance of 2 to the cut site
				df['target_distance'] = [is_gRNA_valid([r[0],r['cut']],[self.chr,self.mutation_pos],r[5]) for i,r in df.iterrows()]
				# print (df)
																	# 0         1         2  ...         		           5    cut                 target_distance
				# chr20_55990386_55990406_+_CGCAGGGGCCGATAAGTGTC  chr20  55990386  55990406  ...  +  55990403               2
				# chr20_55990397_55990417_-_AAGCAGAGCCAGACACTTAT  chr20  55990397  55990417  ...  -  55990401              -1
				# chr20_55990388_55990408_-_CAGACACTTATCGGCCCCTG  chr20  55990388  55990408  ...  -  55990392              -1
				# chr20_55990387_55990407_-_AGACACTTATCGGCCCCTGC  chr20  55990387  55990407  ...  -  55990391              -1
				# chr20_55990384_55990404_-_CACTTATCGGCCCCTGCGGG  chr20  55990384  55990404  ...  -  55990388              -1
				# chr20_55990378_55990398_-_TCGGCCCCTGCGGGAGGCCC  chr20  55990378  55990398  ...  -  55990382              -1		
				
				## gRNA validation given target mutation
				if debug > 5:
					df.to_csv("%s/%s.init.all_sgRNAs.bed"%(self.debug_folder,self.name),sep="\t",header=False,index=False)

				self.valid_init_sgRNA = df[df.target_distance.between(0,max_target_to_sgRNA)][[0,1,2,3,4,5,'cut']]
				
				
				## sgRNA features
				self.sgRNA_target_distance_dict = df['target_distance'].to_dict()
				if debug > 5:
					print (df[df.target_distance.between(0,max_target_to_sgRNA)])
				df = df.drop(['target_distance'],axis=1)
				if self.valid_init_sgRNA.shape[0] == 0:
					print ("No sgRNA was found for %s using %s gRNA_search_space"%(self.name,extend))
					continue
				else:
					self.found_PE2 = True
					print ("%s valid sgRNAs found for  %s"%(self.valid_init_sgRNA.shape[0],self.name))
					self.dist_dict = distance_matrix(df.values.tolist())
					self.sgRNA_strand_df['+'] = df[df[5]=="+"][[0,1,2,3,4,5]]
					self.sgRNA_strand_df['-'] = df[df[5]=="-"][[0,1,2,3,4,5]]
					self.all_sgRNA = df.copy()
					self.sgRNA_target_dPAM_dict = {i: is_dPAM(self.PAM, self.target_pos,self.ref,self.alt,r[0:4].tolist()) for i, r in self.valid_init_sgRNA.iterrows()}
					
					
					break
					
			except Exception as e:
				print (e)
				print ("Error or No sgRNA was found for %s using %s gRNA_search_space"%(self.name,extend))

		if debug > 5:
			print ("Target name: ",self.name)
			print (self.valid_init_sgRNA.head().to_string(index=False))


	def search(self,debug=0,scaffold=None,**kwargs):
		"""Second step: search for all possible PBS, RTS, pegRNA, nick-gRNA combos
		
		Input
		-----
		
		length min and max to define search space
		
		Output
		------
		
		1. valid sgRNA list
			- PBS dataframe
			- RTT dataframe
			- ngRNA dataframe
		
		
		"""
		if not self.found_PE2:
			return 0
			

		self.sgRNA_list = [sgRNA(
								chr = x[0],
								start = x[1],
								end = x[2],
								seq = x[3],
								sgRNA_name = x[4],
								strand = x[5],
								cut_position = x[6],
								mutation_pos = self.mutation_pos,mutation_ref = self.mutation_ref,mutation_alt = self.mutation_alt,
								user_target_pos = self.target_pos,user_ref = self.ref,user_alt = self.alt,
								is_dPAM = self.sgRNA_target_dPAM_dict[x[4]],target_to_sgRNA = self.sgRNA_target_distance_dict[x[4]],
								variant_id = self.name,
								dist_dict = self.dist_dict,
								opposite_strand_sgRNAs = self.sgRNA_strand_df[get_opposite_strand(x[5])],
								all_sgRNA_df = self.all_sgRNA,
								target_fa = self.target_fa,
								scaffold_seq = scaffold
								) 
						for x in self.valid_init_sgRNA.values.tolist()]

		[run_sgRNA_search(s,**kwargs) for s in self.sgRNA_list]
		
		self.rawX = pd.concat([s.rawX for s in self.sgRNA_list])
		if self.rawX.shape[0]==0:
			self.found_PE2=False
			return 0
		self.X = pd.concat([s.X for s in self.sgRNA_list])
		no_ngRNA = sum([s.no_ngRNA for s in self.sgRNA_list])
		if no_ngRNA == len(self.sgRNA_list):
			print ("%s only PE2 found"%(self.name))
		else:
			self.found_PE3 = True
			
		
		self.X['N_insertions'] = self.ref_alt[2]
		self.X['N_subsitution'] = self.ref_alt[1]
		self.X['N_deletion'] = self.ref_alt[3]
		

		self.found_PE3b = (self.X['is_PE3b']==1).any()
		self.found_dPAM = (self.X['is_dPAM']==1).any()		


	def predict(self,debug=0,ML_model=None,**kwargs):
		if not self.found_PE2:
			return 0
		with open(ML_model, 'rb') as file:  
			xgb_model = pickle.load(file)		

		pred_y = xgb_model.predict(self.X[self.feature_for_prediction])

		myPred = pd.DataFrame()
		myPred['predicted_efficiency'] = pred_y.tolist()
		myPred.index = self.X.index.tolist()
		self.X_p = pd.concat([self.X,myPred],axis=1)
		self.rawX['predicted_efficiency'] = myPred.loc[self.rawX.index]['predicted_efficiency']
		
		self.X_p = self.X_p.sort_values("predicted_efficiency",ascending=False)
		self.rawX = self.rawX.sort_values("predicted_efficiency",ascending=False)
		self.topX = self.rawX .loc[self.rawX.index[0]]
		
	

def run_sgRNA_search(s,**kwargs):
	s.find_RTT(**kwargs)
	s.find_PBS(**kwargs)
	s.find_nick_gRNA(**kwargs)
	s.get_rawX_and_X(**kwargs)


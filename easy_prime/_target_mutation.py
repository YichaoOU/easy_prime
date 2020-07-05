
# from .utils import *
# from ._sgRNA import sgRNA
# from .featurize import get_X_feature
# from .featurize import sequence_kmer
# from .featurize import global_alignments

from utils import *
from _sgRNA import sgRNA
from featurize import get_X_feature
from featurize import sequence_kmer
from featurize import global_alignments

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
	"""
	count=0
	for i in range(min(len(ref),len(alt))):
		x=ref[i]
		y=alt[i]
		if x != y:
			return pos,ref[i:],alt[i:]
		else:
			pos+=1
			count+=1
	return pos,ref[count:],alt[count:]


class target_mutation:
	def __init__(self,chr,pos,name,ref,alt,**kwargs):
		"""target mutation class
		
		
		sgRNA name: chr_start_end_strand_seq
		
		target_mutation name: id_chr_pos_ref_alt
		
		pos is corrected, and the corrected pos, ref, alt is used 
		
		
		
		"""
		self.chr = chr
		self.target_pos = pos
		self.name = name.replace("/","_").replace(",","_")
		self.ref = ref
		self.alt = alt
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
		

		
		
		
		
		self.mutation_pos,self.mutation_ref,self.mutation_alt = find_mutation_pos(pos,ref,alt)
		
		
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

		
	def init(self,gRNA_search_space=200,search_iteration=1,sgRNA_length=20,PAM="NGG",offset=-3,debug=0,genome_fasta=None,max_RTT_length=40,min_distance_RTT5=5,**kwargs):
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
		
		self.target_fa
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
		self.target_extend_space = 1000
		self.target_fa = get_fasta_single(self.chr,self.mutation_pos-self.target_extend_space,self.mutation_pos+self.target_extend_space,genome_fasta=genome_fasta).upper() # only call bedtools once
		
		
		### find all sgRNA given a sequence
		for i in range(search_iteration):
			extend = gRNA_search_space*(i+1)
			if i >=1:
				print ("No sgRNA were found using %s gRNA_search_space"%(extend))
			start = self.mutation_pos-extend
			end = self.mutation_pos+extend
			search_fa = sub_fasta_single(self.target_fa,self.mutation_pos, start,end)
			df = run_pam_finder(search_fa,"N"*sgRNA_length,self.PAM,start,self.chr)
			## this above df contains all sgRNAs
			
			
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

				self.valid_init_sgRNA = df[df.target_distance.between(0,max_RTT_length-min_distance_RTT5)][[0,1,2,3,4,5]]
				
				
				## sgRNA features
				self.sgRNA_target_distance_dict = df['target_distance'].to_dict()
				
				if self.valid_init_sgRNA.shape[0] == 0:
					print ("No sgRNA was found for %s using %s gRNA_search_space"%(self.name,extend))
					continue
				else:
					self.found_PE2 = True
					print ("%s valid sgRNAs found for  %s"%(self.valid_init_sgRNA.shape[0],self.name))
					self.dist_dict = distance_matrix(df.values.tolist())
					self.sgRNA_strand_df['+'] = df[df[5]=="+"][[0,1,2,3,4,5]]
					self.sgRNA_strand_df['-'] = df[df[5]=="-"][[0,1,2,3,4,5]]
					self.all_sgRNA = df[[0,1,2,3,4,5]]
					self.sgRNA_target_dPAM_dict = {i: is_dPAM(self.PAM, self.target_pos,self.ref,self.alt,r[0:4].tolist()) for i, r in self.valid_init_sgRNA.iterrows()}
					
					
					break
					
			except Exception as e:
				print (e)
				print ("Error or No sgRNA was found for %s using %s gRNA_search_space"%(self.name,extend))

		if debug > 5:
			print ("Target name: ",self.name)
			print (self.valid_init_sgRNA.head().to_string(index=False))


	def search(self,debug=0,**kwargs):
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
		sgRNA_list = [sgRNA(x[0],x[1],x[2],x[3],x[5],x[-1],self.pos,self.ref,self.alt,self.name,self.dist_dict,self.strand_dict,self.candidate_pegRNA_df,self.max_nick_distance) for x in self.candidate_pegRNA_df.values.tolist()]

		[run_sgRNA_search(s,self.target_fa,**kwargs) for s in sgRNA_list]
		
		self.rawX = pd.concat([s.rawX for s in sgRNA_list])
		self.X = pd.concat([s.X for s in sgRNA_list])

	
	def predict(self,debug=0,model,**kwargs):
		if self.rawX.shape[0] == 0:
			# print ("No valid pegRNAs were found for %s"%(self.name))
			return 0
		with open(model, 'rb') as file:  
			xgb_model = pickle.load(file)		

	
		pred_y = xgb_model.predict(self.X)

		myPred = pd.DataFrame()
		myPred['predicted_efficiency'] = pred_y.tolist()
		myPred.index = self.X.index.tolist()
		self.X_p = pd.concat([self.X,myPred],axis=1)
		self.X_p = self.X_p.sort_values("predicted_efficiency",ascending=False)

	def featurize(self,debug=0,N_combinations_for_optimization=800,**kwargs):
		"""
		
		input
		------
		
		a list of sgRNA object
		
		1. perform featurize of each PBS, RTT, ngRNA features
		2. make combinations
		
		
		Output
		--------
		
		1. rawX format
		
		2. X format
		
		"""
		if self.rawX.shape[0] == 0:
			# print ("No valid pegRNAs were found for %s"%(self.name))
			return False
		rawX = self.rawX.copy()
		rawX = rawX.drop_duplicates('sample_ID')
		rawX['unique_name'] = [x.split("_candidate_")[0].replace(self.name,"") for x in rawX['sample_ID'].tolist()]
		rawX['distance'] = rawX['unique_name'].map(self.sgRNA_target_distance_dict)
		n_sgRNA = rawX['unique_name'].nunique()	
		if self.rawX['sample_ID'].nunique() > N_combinations_for_optimization:
			rawX = rawX.sort_values('distance')
			select_list = rawX.head(n=N_combinations_for_optimization)['sample_ID'].tolist()
			print ("Total number of combinations is above threshold. \n Selecting top %s combinations based on sgRNA to target mutation distance"%(N_combinations_for_optimization))
			self.rawX = self.rawX[self.rawX['sample_ID'].isin(select_list)]
		print ("There are total %s combinations (%s unique sgRNAs) for %s"%(self.rawX['sample_ID'].nunique(),n_sgRNA,self.name))	
		## use for loop
		# feature_sample_list = Parallel(n_jobs=1)(delayed(get_X_feature)(s,d,self.seq_kmer,self.ref_alt,**kwargs) for s,d in self.rawX.groupby("sample_ID"))
		feature_sample_list = [get_X_feature(s,d,self.seq_kmer,self.ref_alt,**kwargs) for s,d in self.rawX.groupby("sample_ID")]
		out = pd.concat(feature_sample_list,axis=1)
		out = pd.DataFrame(out.T)	
		# print (out.head())
		self.X = out.copy()
		if debug>0:
			out.to_csv("%s/%s.X.csv"%(self.debug_folder,self.name))
	
	
		pass
		
	



def run_sgRNA_search(s,**kwargs):
	s.find_RTT(**kwargs)
	s.find_PBS(**kwargs)
	s.find_nick_gRNA(**kwargs)


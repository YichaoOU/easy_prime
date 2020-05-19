
from .utils import *
from ._sgRNA import sgRNA
from .featurize import get_X_feature
from .featurize import sequence_kmer
from .featurize import global_alignments


class target_mutation:
	def __init__(self,chr,pos,name,ref,alt):
		self.chr = chr
		self.pos = pos

		self.name = name.replace("/","_")
		self.ref = ref
		self.alt = alt
		self.candidate_pegRNA_df = pd.DataFrame()
		self.debug_folder = "easy_prime_debug_files"
		self.dist_dict = {}
		self.strand_dict = {}
		self.rawX = pd.DataFrame()
		self.X = pd.DataFrame()
		self.X_p = pd.DataFrame()
		self.topX = pd.DataFrame()
		self.allX = pd.DataFrame()
		
		
	def init(self,genome_fasta=None,gRNA_search_space=500,search_iteration=1,sgRNA_length=20,PAM="NGG",offset=-3,debug=0,PCA_model=None,**kwargs):
		"""First step: search for candidate sgRNAs around target mutation
		
		Input
		-----
		
		gRNA_search_space: extend pos by +- gRNA_search_space
		
		search_iteration: if in the search space defined by gRNA_search_space, we fail to find sgRNAs, we will extend the gRNA_search_space further to find at least one sgRNA.
		
		Output
		------
		
		candidate_pegRNA_df
		
		
		"""

		pos = self.pos
		chr = self.chr
		temp_file_list = []
		self.PAM = PAM
		self.offset = offset
		
		# for featurize later
		with open(PCA_model, 'rb') as file:  
			pca_3mer = pickle.load(file)		
		local_seq = get_fasta_single(chr,pos-50,pos+50,genome_fasta)
		kmer = sequence_kmer(3,local_seq.upper())
		self.seq_kmer = list(pca_3mer.transform(kmer.T)[0])	
		self.ref_alt = global_alignments(self.ref,self.alt)
		# print (self.seq_kmer)
		# print (self.ref_alt)
		# exit()
		## debug folder
		if debug>0:
			os.system("mkdir -p %s"%(self.debug_folder))

		for i in range(search_iteration):
			if i >=1:
				print ("No sgRNA were found using %s gRNA_search_space"%(gRNA_search_space))
			extend = gRNA_search_space*(i+1)
			start = pos-extend
			end = pos+extend
			
			extended_file = to_bed3(chr,start,end)
			temp_file_list.append(extended_file)
			
			extended_fa = get_fasta_given_bed(genome_fasta,extended_file)	
			temp_file_list.append(extended_fa)
			
			cas_output = run_casOFFinder(extended_fa,PAM,["N"*sgRNA_length])
			temp_file_list.append(cas_output)
			try:
				df = cas_local_to_df(cas_output,PAM,sgRNA_length)
				df.columns = [0,1,2,3,4,5]
				df[1] = df[1].astype(int)
				df[2] = df[2].astype(int)
				df['cut'] = [get_gRNA_cut_site(x[1],x[2],x[5],self.offset) for i,x in df.iterrows()]
				df['flag'] = [is_gRNA_valid([r[0],r['cut']],[self.chr,self.pos],r[5]) for i,r in df.iterrows()]
				## gRNA validation given target mutation
				if debug > 5:
					df.to_csv("%s/%s.init.all_sgRNAs.bed"%(self.debug_folder,self.name),sep="\t",header=False,index=False)
				
				df = df[df['flag']==True]
				df = df.drop(['flag'],axis=1)
				df[4] = df[0]+"_"+df[1].astype(str)+"_"+df[2].astype(str)+"_"+df[3].astype(str)
				df.index = df[4].to_list()

				if df.shape[0] == 0:
					print ("No sgRNA was found for %s using %s gRNA_search_space"%(self.name,gRNA_search_space))
					continue
				else:
					

					self.candidate_pegRNA_df = df.copy()
					self.dist_dict = distance_matrix(df.values.tolist(),self.offset)
					self.strand_dict = df[5].to_dict()

					break
			except Exception as e:
				print (e)
				print ("Error or No sgRNA was found for %s using %s gRNA_search_space"%(self.name,gRNA_search_space))
				
		

		if debug > 5:
			print ("Target name: ",self.name)
			print (self.candidate_pegRNA_df.head())
	
		delete_files(temp_file_list)
		
	def search(self,debug=0,**kwargs):
		"""Second step: search for all possible PBS, RTS, pegRNA, nick-gRNA combos
		
		Input
		-----
		
		length min and max to define search space
		
		Output
		------
		
		rawX format
		
		
		"""
		sgRNA_list = [sgRNA(x[0],x[1],x[2],x[3],x[5],x[-1],self.pos,self.ref,self.alt,self.name,self.dist_dict,self.strand_dict,self.candidate_pegRNA_df) for x in self.candidate_pegRNA_df.values.tolist()]
		Parallel(n_jobs=1)(delayed(run_sgRNA_search)(s,**kwargs) for s in sgRNA_list)
		df = pd.concat([s.rawX for s in sgRNA_list])
		self.rawX = df.copy()
		if debug>1:
			df.to_csv("%s/%s.rawX.csv"%(self.debug_folder,self.name))

		pass
	
	def predict(self,debug=0,model1=None,model2=None,w1=None,w2=None,**kwargs):
		if self.rawX.shape[0] == 0:
			# print ("No valid pegRNAs were found for %s"%(self.name))
			return False
		with open(model1, 'rb') as file:  
			m1 = pickle.load(file)		
		with open(model2, 'rb') as file:  
			m2 = pickle.load(file)		
	
		pred_y1 = m1.predict(self.X)
		pred_y2 = m2.predict(self.X)
		myPred = pd.DataFrame()
		myPred['model1_predction'] = pred_y1.tolist()
		myPred['model2_predction'] = pred_y2.tolist()
		myPred['predicted_efficiency'] = myPred['model1_predction']*w1+myPred['model2_predction']*w2
		myPred.index = self.X.index.tolist()
		X_p = pd.concat([self.X,myPred],axis=1)
		X_p = X_p.sort_values("predicted_efficiency",ascending=False)
		self.X_p = X_p
		if debug>0:
			X_p.to_csv("%s/%s.X_p.csv"%(self.debug_folder,self.name))
		
			
		pass
		
	def featurize(self,debug=0,**kwargs):
		if self.rawX.shape[0] == 0:
			print ("No valid pegRNAs were found for %s"%(self.name))
			return False
		print ("There are total %s combinations for %s"%(self.rawX['sample_ID'].nunique(),self.name))
		feature_sample_list = Parallel(n_jobs=1)(delayed(get_X_feature)(s,d,self.seq_kmer,self.ref_alt,**kwargs) for s,d in self.rawX.groupby("sample_ID"))
		out = pd.concat(feature_sample_list,axis=1)
		out = pd.DataFrame(out.T)	
		# print (out.head())
		self.X = out.copy()
		if debug>0:
			out.to_csv("%s/%s.X.csv"%(self.debug_folder,self.name))
	
	
		pass
		
	

	def output(self,N_top_pegRNAs=3,debug=0,**kwargs):
		if self.rawX.shape[0] == 0:
			print ("No valid pegRNAs were found for %s"%(self.name))
			return False	
		df = self.X_p.copy()
		df['unique_name'] = [x.split("_candidate_")[0] for x in df.index.tolist()]
		top_list=[]
		for s,d in df.groupby('unique_name'):
			top_list+=d.sort_values('predicted_efficiency',ascending=False).head(n=N_top_pegRNAs).index.tolist()
		# top_list = self.X_p.head(n=N_top_pegRNAs).index.tolist()
		topX = self.rawX[self.rawX['sample_ID'].isin(top_list)]
		allX = self.rawX.copy()
		# print (self.rawX.head())
		
		# print (top_list)
		topX['predicted_efficiency'] = self.X_p.loc[topX['sample_ID'].tolist()]['predicted_efficiency'].tolist()
		allX['predicted_efficiency'] = self.X_p.loc[allX['sample_ID'].tolist()]['predicted_efficiency'].tolist()
		self.topX = topX
		self.allX = allX
		if debug>0:
			topX.to_csv("%s/%s.topX.csv"%(self.debug_folder,self.name))
			
		# return topX
		# print (self.topX.head())
		
		pass


def run_sgRNA_search(s,**kwargs):
	s.find_nick_gRNA(**kwargs)
	s.find_RTS(**kwargs)
	s.find_PBS(**kwargs)
	s.add_variant(**kwargs)

	
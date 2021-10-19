import shap
import warnings
import re
warnings.filterwarnings("ignore")
import os 
import argparse
import shutil
import datetime
import getpass
import uuid
import pandas as pd
import yaml
import itertools
import pickle
import subprocess
from skbio.alignment import global_pairwise_align_nucleotide
from skbio import DNA
import RNA
from copy import deepcopy as dp
import numpy as np
import sys
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.pylab as plt
import scipy
import glob

from scipy import interp
import matplotlib
matplotlib.use('agg')
import xgboost as xgb
from sklearn.model_selection import KFold,StratifiedKFold
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression,RidgeClassifier,SGDClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB 
from sklearn.ensemble import RandomForestClassifier
from mlxtend.classifier import StackingCVClassifier
import umap
from sklearn.metrics import roc_curve,roc_auc_score,average_precision_score
from sklearn.datasets import load_iris
from mlxtend.feature_selection import ColumnSelector
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
from sklearn.ensemble import RandomForestRegressor,GradientBoostingRegressor,RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
from sklearn.metrics.scorer import make_scorer
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestRegressor,GradientBoostingRegressor
from sklearn.base import TransformerMixin
from sklearn.datasets import make_regression
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import PolynomialFeatures
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_absolute_error
from sklearn import linear_model
from sklearn.kernel_ridge import KernelRidge
from sklearn.svm import SVR
from mlxtend.regressor import StackingCVRegressor
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge,Lars,BayesianRidge
from xgboost import XGBRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.ensemble import RandomForestRegressor,GradientBoostingClassifier
from sklearn.svm import SVR,LinearSVC
from sklearn.neighbors import KNeighborsClassifier,RadiusNeighborsClassifier
from xgboost import XGBClassifier
from decimal import Decimal
from joblib import Parallel, delayed
from sklearn.model_selection import ParameterGrid
from sklearn.metrics import mean_squared_error, r2_score,max_error,mean_squared_error,mean_squared_log_error



def plot_correlation(plot_df,output):
	df2 = plot_df.copy()
	df2['group'] = [x.split("_REP")[0] for x in df2.index.tolist()]
	df2 = df2.groupby('group').mean()
	# print ("===54756=")
	# return 0
	plt.figure()
	# sns.regplot(x=plot_df['true'],y=plot_df['pred'],x_bins =20)
	sns.regplot(x=plot_df['true'],y=plot_df['pred'],scatter_kws={"color":"black","alpha":0.7,"s":3},line_kws={"color":"orange","alpha":0.8})
	r,p = scipy.stats.pearsonr(plot_df["pred"],plot_df["true"])
	plt.xlabel("True value")
	plt.ylabel("Predicted value")
	max_value = plot_df.max().max()
	min_value = plot_df.min().min()
	
	plt.text(max_value*0.1, max_value*0.85, 'r=%f'%(r))
	plt.text(max_value*0.1, max_value*0.8, 'p=%.2E'%(Decimal(p)))
	plt.xlim(0,min(max_value+0.1,1))
	plt.ylim(0,min(max_value+0.1,1))
	
	plt.savefig("%s_correlation_plot.pdf"%(output), bbox_inches='tight')
	return r
	
def get_correlation(plot_df):
	# print ("===54756=")
	# return 0
	# plt.figure()
	# sns.regplot(x=plot_df['true'],y=plot_df['pred'],x_bins =20)
	# sns.regplot(x=plot_df['true'],y=plot_df['pred'],scatter_kws={"color":"black","alpha":0.7,"s":3},line_kws={"color":"orange","alpha":0.8})
	r,p = scipy.stats.pearsonr(plot_df["pred"],plot_df["true"])
	# plt.xlabel("True value")
	# plt.ylabel("Predicted value")
	# max_value = plot_df.max().max()
	# min_value = plot_df.min().min()
	
	# plt.text(max_value*0.1, max_value*0.85, 'r=%f'%(r))
	# plt.text(max_value*0.1, max_value*0.8, 'p=%.2E'%(Decimal(p)))
	# plt.xlim(0,min(max_value+0.1,1))
	# plt.ylim(0,min(max_value+0.1,1))
	
	# plt.savefig("%s_correlation_plot.pdf"%(output), bbox_inches='tight')
	return r
	

def plot_top_features(reg,X,y,output,n):
	
	current_feature_df = pd.DataFrame()
	current_feature_df['features'] = X.columns.tolist()

	reg.fit(X,y)
	try:
		current_feature_df['score'] = list(reg.feature_importances_) 
	except:
		try:
			current_feature_df['score'] = list(reg.coef_) 
		except:
			current_feature_df['score'] = list(reg.coef_[0]) 
				
	current_feature_df = current_feature_df.sort_values('score',ascending=False)
	if output != None:
		current_feature_df.to_csv("%s_feature_ranking.tsv"%(output),sep="\t")
	if n != None:
		current_feature_df = current_feature_df.head(n=n)
	

	# plt.figure(figsize=(max(5,n/2),5))
	# sns.barplot(x=current_feature_df['features'],y=current_feature_df['score'] )
	# plt.xticks(rotation=90)
	# plt.xlabel("")
	# plt.ylabel("Feature importance")
	# plt.savefig("%s_feature_importance.pdf"%(output), bbox_inches='tight')
	return current_feature_df
	# return current_feature_df.features.tolist()

def plot_top_features1(reg1,reg2,X,y,output,n):
	
	current_feature_df = pd.DataFrame()
	current_feature_df['features'] = X.columns.tolist()

	reg1.fit(X,y)
	try:
		current_feature_df['score1'] = list(reg1.feature_importances_) 
	except:
		try:
			current_feature_df['score1'] = list(reg1.coef_) 
		except:
			current_feature_df['score1'] = list(reg1.coef_[0]) 

	reg2.fit(X,y)
	try:
		current_feature_df['score2'] = list(reg2.feature_importances_) 
	except:
		try:
			current_feature_df['score2'] = list(reg2.coef_) 
		except:
			current_feature_df['score2'] = list(reg2.coef_[0]) 

				
	current_feature_df = current_feature_df.sort_values('score',ascending=False)
	current_feature_df.to_csv("%s_feature_ranking.tsv"%(output),sep="\t")
	current_feature_df = current_feature_df.head(n=n)
	

	plt.figure(figsize=(max(5,n/2),5))
	sns.barplot(x=current_feature_df['features'],y=current_feature_df['score'] )
	plt.xticks(rotation=90)
	plt.xlabel("")
	plt.ylabel("Feature importance")
	plt.savefig("%s_feature_importance.pdf"%(output), bbox_inches='tight')

def get_importance(reg,features):
	
	current_feature_df = pd.DataFrame()
	try:
		current_feature_df['score'] = list(reg.feature_importances_) 
	except:
		try:
			current_feature_df['score'] = list(reg.coef_) 
		except:
			current_feature_df['score'] = list(reg.coef_[0]) 
	current_feature_df.index = features
	current_feature_df = current_feature_df.sort_values('score',ascending=False)
	return current_feature_df




def MAE_grid(estimator,parameterDict,X,y):
    myModel = GridSearchCV(estimator,parameterDict,scoring="neg_mean_absolute_error",cv=3,verbose=10,n_jobs=10,refit=True)
    myModel.fit(X,y)
    return myModel


def sklearn_GBT_reg(par=False):
	est = GradientBoostingRegressor(n_estimators=100)
	if par:
		est = GradientBoostingRegressor(**par)
	myDict = {}
	myDict['loss']=['ls','lad']
	myDict['learning_rate']=[0.05,0.1]
	myDict['max_depth']=[5,8]
	return est, myDict

def sklearn_RF_reg(par=False):
	est = RandomForestRegressor(n_estimators=300)
	if par:
		est = RandomForestRegressor(**par)
	myDict = {}
	myDict['max_depth']=[2,5,9]
	myDict['ccp_alpha']=[0,0.01]
	myDict['max_features']=['log2','auto']
	return est, myDict
	
	
def get_2_weight():
	myList = []
	for i in range(11):
		w1 = i/10.0
		for j in range(11):
			w2 = j/10.0
			if w1+w2 == 1:
				myList.append([w1,w2])
	return myList

def get_model_list():
	model_list = {}
	model_list['sklearn_RF_reg'] = sklearn_RF_reg
	model_list['sklearn_GBT_reg'] = sklearn_GBT_reg
	model_list['xgb_reg'] = xgb_reg
	# model_list['ada_RF_reg'] = ada_RF_reg 
	return model_list

class linear_blending2:
	def __init__(self,m1,m2,p1,p2,w1,w2,existing_model=False):
		if existing_model == False:
			model_list = get_model_list()
			self.m1,_ = model_list[m1](p1)
			self.m2,_ = model_list[m2](p2)
			self.w1 = w1
			self.w2 = w2
		else:
			self.m1,self.m2,self.w1,self.w2 = self.read_model(existing_model)
	def read_model(self,input):
	
		pass
		
	def fit(self,X,y):
		self.m1.fit(X)
		self.m2.fit(X)
		self.features = X.columns.tolist()

	def predict(self,X,Y=None):
		pred1 = self.m1.predict(X)
		pred2 = self.m2.predict(X)
		pred = pred1*self.w1+pred2*self.w2
		return pred.tolist()

	def compute_importance(self):
		df1 = get_importance(self.m1,self.features)
		df2 = get_importance(self.m2,self.features)
		df = pd.concat([df1,df2])
		df.columns = ['m1','m2']
		df['importance'] = df.m1*self.w1 + df.m2*self.w2
		df = df.sort_values('importance',ascending=False)
		return df



	
def xgb_reg(par=False):
	est = XGBRegressor(n_estimators=300,objective="reg:squarederror")
	if par:
		est = XGBRegressor(seed=100,**par)
	myDict = {}
	myDict['max_depth']=[2,5,9]
	myDict['learning_rate'] = [0.01,0.1]
	myDict['min_child_weight']=[1,5,10]
	myDict['colsample_bylevel']=[0.2,0.6,1]
	myDict['colsample_bytree']=[0.2,0.6,1]
	myDict['subsample']=[0.2,0.6,1]
	myDict['reg_alpha']=[0,0.1,1]
	myDict['reg_lambda']=[0,1,2]

	return est, myDict

def ridge_reg(par=False):
	est = Ridge()
	if par:
		est = Ridge(**par)	
	myDict = {}
	myDict['alpha']=[0.1,0.5,1,2]
	return est, myDict
def lasso_reg(par=False):
	est = Lasso()
	if par:
		est = Lasso(**par)	
	myDict = {}
	myDict['alpha']=[0.1,0.5,1,2]
	return est, myDict
def Lars_reg(par=False):
	est = Lars()
	if par:
		est = Lars(**par)	
	myDict = {}
	myDict['n_nonzero_coefs']=[1,5,10,50,100,200,300]
	return est, myDict
def BayesianRidge_reg(par=False):
	est = BayesianRidge()
	if par:
		est = BayesianRidge(**par)	
	myDict = {}
	myDict['alpha_1']=[1e-10,1e-5,1e-3,0.1]
	myDict['alpha_2']=[1e-10,1e-5,1e-3,0.1]
	myDict['lambda_1']=[1e-10,1e-5,1e-3,0.1]
	myDict['lambda_2']=[1e-10,1e-5,1e-3,0.1]
	return est, myDict
def KernelRidge_reg(par=False):
	## dual_coef_ is not working
	est = KernelRidge()
	if par:
		est = KernelRidge(**par)	
	myDict = {}
	myDict['alpha']=[1e-5,1e-3,0.1,0.5,1,2]
	myDict['gamma']=[1e-5,1e-3,0.1,0.5,1,2]
	myDict['kernel']=['linear','rbf']

	return est, myDict
def SVM_reg(par=False):
	est = SVR(kernel="linear")
	if par:
		est = SVR(**par)	
	myDict = {}
	myDict['C']=[1e-5,1e-3,0.1,0.5,1,2,10,100]
	myDict['gamma']=[1e-5,1e-3,0.1,0.5,1,2]
	myDict['kernel']=['linear']

	return est, myDict

def linear_stacking_reg(X,y):
	RANDOM_SEED=0
	reg1,dict1 = SVM_reg()
	reg2,dict2 = sklearn_GBT_reg()
	reg3,dict3 = xgb_reg()
	reg4,dict4 = ridge_reg()
	reg5,dict5 = lasso_reg()
	reg6,dict6 = Lars_reg()
	reg7,dict7 = KNeighborsRegressor_reg()
	reg8,dict8 = sklearn_RF_reg()
	reg9,dict9 = KernelRidge_reg()
	reg10,dict10 = BayesianRidge_reg()
	params ={}
	my_dict_list = [dict1,dict2,dict3,dict4,dict5,dict6,dict7,dict8,dict9,dict10]
	name_list = ["svr","gradientboostingregressor","xgbregressor",
				"ridge","lasso","lars","kneighborsregressor","randomforestregressor",
				"kernelridge","meta_regressor"]
	for i in range(len(my_dict_list)):
		current_dict = my_dict_list[i]
		for xx in current_dict:
			if i == 9:
				continue
			params['%s__%s'%(name_list[i],xx)]=current_dict[xx]
	linear_stacker = StackingCVRegressor(regressors=(reg1,reg2,reg3,reg4,reg5,reg6,reg7,reg8,reg9),
					meta_regressor=reg10,use_features_in_secondary=True)

	grid = RandomizedSearchCV(linear_stacker,params,n_iter=10,n_jobs=10,scoring="neg_mean_absolute_error",cv=3,verbose=10)

	linear_stacker.fit(X.values,y)
	# grid.fit(X.values,y)
	
	# print("Best: %f using %s" % (grid.best_score_, grid.best_params_))
	# best_model = grid.best_estimator_ 
	
	return linear_stacker
	
	

def KNeighborsRegressor_reg(par=False):
	est = KNeighborsRegressor()
	if par:
		est = KNeighborsRegressor(**par)	
	myDict = {}
	myDict['n_neighbors']=[1,2,3,5,10,20]
	myDict['leaf_size']=[1,5,10,20]
	myDict['p']=[1,2]
	myDict['weights']=['uniform','distance']

	return est, myDict
											
def recursive_feature_elimination_reg(reg,X,y):
	features = X.columns
	score_list = []
	my_feature_score_df = pd.DataFrame()
	feature_list = []
	for i in range(len(features-1)):
		X = X[features]
		feature_list.append(features)
		current_feature_df = pd.DataFrame()
		current_feature_df['features'] = features
		# print (current_feature_df.shape)
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)
		reg.fit(X_train,y_train)
		y_pred = reg.predict(X_test)
		MAE_score = mean_absolute_error(y_test, y_pred)
		score_list.append(MAE_score)
		# print (i,MAE_score)
		# print (len(reg.feature_importances_ ))
		try:
			current_feature_df['score'] = list(reg.feature_importances_) 
		except:
			try:
				current_feature_df['score'] = list(reg.coef_) 
			except:
				# print (reg.coef_[0])
				# current_feature_df['score'] = list(reg.dual_coef_) 
				current_feature_df['score'] = list(reg.coef_[0]) 
				
		current_feature_df = current_feature_df.sort_values('score')
		current_feature_df = current_feature_df.reset_index(drop=True)
		current_feature_df = current_feature_df.drop(0)
		features = current_feature_df['features'].tolist()
	my_feature_score_df['score'] = 	score_list
	my_feature_score_df['features'] = 	feature_list
	my_feature_score_df = my_feature_score_df.sort_values('score')
	best_MAE = my_feature_score_df['score'].tolist()[0]
	best_featureSet = my_feature_score_df['features'].tolist()[0]
	print ("Best MAE is %s, using %s features"%(best_MAE,len(best_featureSet)))
	return best_featureSet



def loo_reg(reg,X,y):
	my_pred = []
	my_true = []   
	index_list = []
	for i in X.index:
		print (i)
		train_X = X.drop([i])
		train_y = y.drop([i])
		test_X = pd.DataFrame(X.loc[i]).T
		test_y = y.loc[i]
		reg.fit(train_X.values,train_y)
		pred_y = reg.predict(test_X.values)
		my_pred.append(pred_y[0])
		my_true.append(test_y)
		# print (pred_y,test_y)
		index_list.append(i)
	return my_true,my_pred,index_list

def get_top_features(reg,X,y,top_n):

	current_feature_df = pd.DataFrame()
	current_feature_df['features'] = X.columns.tolist()

	reg.fit(X,y)
	try:
		current_feature_df['score'] = list(reg.feature_importances_) 
	except:
		try:
			current_feature_df['score'] = list(reg.coef_) 
		except:
			current_feature_df['score'] = list(reg.coef_[0]) 
				
	current_feature_df = current_feature_df.sort_values('score',ascending=False)
	return current_feature_df['features'].tolist()[:top_n]
	
def GC_content_gRNA(seq):
	GC=["G","C"]
	count=0
	for i in seq:
		if i in GC:
			count+=1
	return float(count+0.5)/20
	pass
	
def GC_content_raw(seq):
	GC=["G","C"]
	count=0
	for i in seq:
		if i in GC:
			count+=1
	return float(count)/len(seq)
	pass
	
    
tab = str.maketrans("ACTG", "TGAC")
def revcomp(seq):
	return seq.translate(tab)[::-1]
		

ulength = 31
md = RNA.md()
md.max_bp_span = 70
md.window_size = 70

def call_RNAplfold(seq,scaffold_length=76):
    """RNA fold python binding
    """
    # scaffold - RTT - PBS
    # we only care about RTT+PBS (3' extension) pairs with scaffold

    data = []
    fc = RNA.fold_compound(seq, md, RNA.OPTION_WINDOW)
    fc.probs_window(ulength, RNA.PROBS_WINDOW_BPP | RNA.PROBS_WINDOW_UP, pf_window_callback2, data)
    
    # parse output
    df = pd.DataFrame(data)
    df2 = df.copy()
    df2[0] = df[1]
    df2[1] = df[0]
    df = pd.concat([df,df2])
    seq_length = 16
    RTT_start = scaffold_length + 1
    RTT_end = scaffold_length + seq_length
    scaffold_start = 1
    scaffold_end =  scaffold_length
    
    ## subset df
    df = df[df[0]>=RTT_start]
    df = df[df[0]<=RTT_end]
    df = df[df[1]>=scaffold_start]
    df = df[df[1]<=scaffold_end]
    df = df[[0,2]]
    df = df.groupby(0).max()
    myDict = df[2].to_dict()
    
    value_list = []
    for i in range(seq_length):
        RTT_pos = scaffold_length +i+ 1
        if RTT_pos in myDict:
            value_list.append(myDict[RTT_pos])
        else:
            value_list.append(0)
    return value_list


def pf_window_callback2(v, v_size, i, maxsize, what, data=None):
    if what & RNA.PROBS_WINDOW_UP:
        pass
    else:
        data+=[[i,j,p] for j, p in enumerate(v) if (p is not None) and (p >= 0.01)]

def local_alignments(ref,q):
    query = StripedSmithWaterman(ref)
    alignment = query(q)
    return alignment['optimal_alignment_score']


def alignments_to_cigar(ref,q):
    Match = 0
    Mis = 0
    D = 0
    I = 0
    for i in range(len(ref)):
        a=ref[i]
        b=q[i]
        if a=="-":
            I=I+1
        elif b=="-":
            D=D+1
        elif a==b:
            Match+=1
        else:
            Mis+=1
    return [Match,Mis,D,I]    

def global_alignments(ref,q):

    s1 = DNA(ref)
    s2 = DNA(q)
    alignment, score, start_end_positions = global_pairwise_align_nucleotide(s1,s2,match_score=4,mismatch_score=1)
    return alignments_to_cigar(alignment[0]._string.decode("utf-8"),alignment[1]._string.decode("utf-8"))


def leave_N_target_mutation_out(y,nfold=10):
    """
    read raw X
    """

    kf = KFold(n_splits=nfold,shuffle=True,random_state=0)



    df = pd.DataFrame(y.copy())
    df2 = pd.read_csv("sample_ID_desired_edit_names.csv",index_col=0)
    df2['group'] = [x.split("_REP")[0] for x in df2.index.tolist()]

    df2 = df2.drop_duplicates('group')
    df2.index = df2['group'].tolist()
    df['desired_edit'] = df2['name']
    df['group'] = [x.split("_REP")[0] for x in df.index.tolist()]
    df = df.drop_duplicates('group')
    df.index = df['group'].tolist()
    index_set = set(df.index.tolist())
    out = []
    group_list = df['desired_edit'].unique().tolist()
    print (group_list)
    # for s,d in df.groupby('desired_edit'):
    for i,j in kf.split(group_list):

        # print ("inside",)
        train_group = [group_list[x] for x in i]
        test_group = [group_list[x] for x in j]
        train = df[df['desired_edit'].isin(train_group)].index.tolist()
        # print ("inside",train)
        test = df[df['desired_edit'].isin(test_group)].index.tolist()


        out.append([train,test])			

    return out	

	
def define_cell(x):
	cells=['K562','U2OS','HEK293T','HELA']
	for c in cells:
		if c in x:
			return c
	return "HEK293T"
	


# cross validation find best parameter and correlation

def train_eval2(X,y,model_name,par):	
    """return best model for the given single parameter combination"""
    my_pred=[]
    my_true=[]
    index_list = []

    for train_index, test_index in leave_N_target_mutation_out(y,3):
        X_train, X_test = X.loc[train_index], X.loc[test_index]
        y_train, y_test = y.loc[train_index], y.loc[test_index]
#         model,_ = model_list[model_name](par)
        model = XGBRegressor(seed=100,**par)
        model.fit(X_train,y_train)
        myPred = model.predict(X_test)
        my_pred += myPred.tolist()
        my_true += y_test.tolist()
        index_list += test_index
    r,p = scipy.stats.pearsonr(my_true,my_pred)
    # r,p = scipy.stats.spearmanr(my_true,my_pred)
    r2 = r2_score(my_true,my_pred)
    mae = mean_absolute_error(my_true,my_pred)	
    return [r,r2,mae]	
	

def nested_cv(X,y,model_name,p_list,nfold=10,debug=0):
	"""return best model for the given single parameter combination"""
	my_pred=[]
	my_true=[]
	index_list = []
	inner_best_row = []
	## outer loop
	for train_index, test_index in leave_N_target_mutation_out(y,nfold):
	# for train_index, test_index in leave_N_sample_out(y,nfold):
		if debug > 10:
			print ("train size: %s test size: %s"%(len(train_index),len(test_index)))
			print (test_index[:2])
		X_train, X_test = X.loc[train_index], X.loc[test_index]
		y_train, y_test = y.loc[train_index], y.loc[test_index]
		
		## inner loop find the best parameter
		result_list = Parallel(n_jobs=-1,verbose=0)(delayed(train_eval2)(X_train,y_train,model_name,p) for p in p_list)
		result_df = pd.DataFrame(result_list)
		result_df['parameters'] = p_list
		# result_df = result_df.sort_values(0,ascending=False) # sort cause error
		result_df = result_df.loc[result_df[1].idxmax()]
		if debug > 2:
			print ("-------inner loop best result---------\n",result_df.to_string())		
		best_line = result_df.values.tolist()
		best_parameter = best_line[-1]
		inner_best_row.append(best_line)
		
		## outer loop testing
		model,_ = model_list[model_name](best_parameter)
		model.fit(X_train,y_train)
		myPred = model.predict(X_test)
		my_pred += myPred.tolist()
		my_true += y_test.tolist()
		index_list += test_index
	outer_df = pd.DataFrame()
	outer_df['true'] = my_true
	outer_df['pred'] = my_pred
	outer_df.index = index_list
	r,p = scipy.stats.pearsonr(my_true,my_pred)
	# r,p = scipy.stats.spearmanr(my_true,my_pred)
	r2 = r2_score(my_true,my_pred)
	n=X.shape[0]
	p=X.shape[1]
	ar2 = 1-(1-r2)*(n-1)/(n-p-1)
	mae = mean_absolute_error(my_true,my_pred)	
	print ("Nested %s-fold CV result. See below"%nfold)
	print ("Pearson correlation: %s"%(r))
	print ("adjusted R2: %s"%(ar2))
	inner_df = pd.DataFrame(inner_best_row)
	inner_df.columns = ['r','r2','mae',"best_parameter"]
	return outer_df,inner_df



def is_dPAM2(PAM_seq, target_pos,ref,alt,pegRNA_loc):
    # currently accept N as ambiguos letter, R will cause error
    # report a bug in Biopython
    # https://github.com/biopython/biopython/issues/3023
    # currently, assuming ref do not contain IUPAC letter
    # get PAM location
    PAM_abs_pos=[]
    PAM_nucleotide=[]
    flag = 0
    if pegRNA_loc[3] == "-":
        for i in range(len(PAM_seq)):
            if PAM_seq[i]!= "N":
                PAM_nucleotide.append(PAM_seq[i])
                PAM_abs_pos.append(pegRNA_loc[1]-i)
    else:
        for i in range(len(PAM_seq)):
            if PAM_seq[i]!= "N":
                PAM_nucleotide.append(PAM_seq[i])
                PAM_abs_pos.append(pegRNA_loc[2]+1+i)
    if len(PAM_abs_pos) == 0:
        # out = pd.DataFrame([flag])
        # out.index = ['is_dPAM']
        return flag
    
    ## same length
    diff = len(ref)-len(alt)
    if diff==0:
        for i in range(len(ref)):
            current_pos = target_pos + i
            if current_pos in PAM_abs_pos:
                flag = 1
                break
    else: # indel
        for i in range(len(alt)):
            current_pos = target_pos + i
            for x in range(len(PAM_abs_pos)):
                PAM_pos = PAM_abs_pos[x]
                PAM_nuc = PAM_nucleotide[x]
                if PAM_pos == current_pos:
                    if alt[i] != PAM_nuc:
                        flag = 1
    # out = pd.DataFrame([flag])
    # out.index = ['is_dPAM']
    # print ("flag",flag)
    return flag
def is_dPAM(PAM_seq, RTT, cut_offset=-3):
    # Assuming no N is RTT, which should be true
    # match PAM seq to RTT, should be abs(cut_offset)
    # print (PAM_seq, RTT)
    # will need to do revcomp no matter what, because RTT is always xxxxxxxPAM

    seq = revcomp(RTT)
    fwd_search = SeqUtils.nt_search(seq, PAM_seq)
    flag = 1
    if len(fwd_search) > 1:
        if abs(cut_offset) in fwd_search:
            flag = 0

    return flag    
    
def target_to_RTT5_feature(pegRNA,nick_gRNA,target_loc,RTS_length,alt_length):
    # 1. target mutation distance to cut 1 (pegRNA)
    # 2. target mutation distance to cut 2 (nick-gRNA)
    # target_loc [chr,pos]
    # 3. cut 1 to cut 2
    cut1 = get_gRNA_cut_site(pegRNA[1],pegRNA[2],pegRNA[3])
    ## if nick_gRNA not exist return a big number
    if nick_gRNA[0] == "0":
        cut2=0
    else:
        cut2 = get_gRNA_cut_site(nick_gRNA[1],nick_gRNA[2],nick_gRNA[3])
    
    # cut2 to cut1
    a=cut2-cut1-1


    # target to cut1
    b=target_loc[1]-cut1
    
    if pegRNA[3]=="-":
        a+=2 # match to coordinate system
        a=-a
        b=-b
    c=cut2-target_loc[1]    
    if nick_gRNA[0] == "0":
        a=np.nan    
        c=np.nan
    if b <0:
        print ("pegRNA to target is less than 0!")
        exit()
    d = RTS_length - alt_length- b + 1 ## number of nucleotide to the RTT 5' end, from the target mutation (not including)
    if d <0:
        print ("pegRNA to target is less than 0!")
        exit()    
    index_list = ["nick_to_pegRNA","target_to_pegRNA","target_to_ngRNA","target_to_RTT5"]
    out = pd.DataFrame([a,b,c,d])
    out.index = index_list
    return out
    

    pass

def plot_ML_scatter(my_true,my_pred,output_file,max_value,min_value=-2):
	outer_df = pd.DataFrame()
	outer_df['true'] = my_true
	outer_df['pred'] = my_pred
	sns.set_style("ticks")
	plt.figure(figsize=(6,6))
	r,p = scipy.stats.pearsonr(my_true,my_pred)
	sr,p = scipy.stats.spearmanr(my_true,my_pred)
	x="true"
	y='pred'
	linewidth=3
	sns.regplot(data=outer_df,x=x,y=y,scatter=True,scatter_kws={'alpha':0.3,'s':4,"color":"grey"},
				line_kws={'linewidth':linewidth,'color':'black'},fit_reg=True)
	plt.legend(['R=%.2f\nr=%.2f'%(sr,r)])
	plt.xlim(min_value,max_value)
	plt.ylim(min_value,max_value)
	plt.xlabel("True efficiency")
	plt.ylabel("Predicted efficiency")
	plt.savefig("%s.corr.scatter.pdf"%(output_file),bbox_inches='tight')




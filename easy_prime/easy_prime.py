#!/usr/bin/env python

from utils import *

"""

"""

def my_args():
	mainParser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,description="easy_prime for pegRNA design")
	
	mainParser.add_argument('-f','--vcf_file',  help="input target mutations to look for pegRNAs",required=True)
	mainParser.add_argument('-c','--config',  help="A YAML file specifying parameters",default=None)
	mainParser.add_argument('-o','--output',  help="output name for csv file",default="easy_prime.output.csv")

	
	##------- add parameters above ---------------------
	args = mainParser.parse_args()	
	return args

def run_steps(t,**kwargs):
	t.init(**kwargs)
	t.search(**kwargs)
	t.featurize(**kwargs)
	t.predict(**kwargs)
	t.output(**kwargs)
	# print (t.topX.head())
	return t.topX

def main():

	args = my_args()
	parameters = get_parameters(args.config)
	
	## get a list of targets
	my_targets = []
	with open(args.vcf_file) as f:
		for line in f:
			line = line.strip()
			if "#" == line[0]:
				continue
			line = line.split()
			chr = line[0]
			pos = int(line[1])
			name = line[2]
			ref = line[3]
			alt = line[4]
			my_targets.append(target_mutation(chr,pos,name,ref,alt))
	
	## find best pegRNAs
	df_list = Parallel(n_jobs=-1)(delayed(run_steps)(t,**parameters) for t in my_targets)

	df = pd.concat(df_list)
	df.to_csv(args.output)
	
	
if __name__ == "__main__":
	main()



























# exec(open("test.py").read())
import utils
import _target_mutation
import imp
imp.reload(utils)
imp.reload(_target_mutation)

from utils import get_parameters, print_parameters
parameters = get_parameters("config.yaml")
print_parameters(parameters)

## get a list of targets
from _target_mutation import target_mutation
import pandas as pd
vcf = pd.read_csv("test.vcf",comment="#",sep="\t",header=None)
for i,r in vcf.iterrows():
	t = target_mutation(*r) 

t.init(**parameters)



import json
import math
from scipy.stats import poisson
import sys
import logging
from probabilistic import *
import numpy as np
from mpmath import *

from functools import reduce

logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")

mp.dps=512

#p_lambda(k) = ( (lambda^k)/(k!) ) * e ^ (- lambda)

calculatedPriors = []

with open(snakemake.input['probabilities'],'r') as infile, open(snakemake.output['checksum'],'w') as outfile:
	vals = [mpf(x.split('\t')[1]) for x in infile.read().splitlines()]
	vals = filter(lambda x : x != -mpf('inf'),vals)
	#print(list(vals))
	cs = reduce(logSpaceAdd,vals)
	#to actual probability
	outfile.write(nstr(cs,5))

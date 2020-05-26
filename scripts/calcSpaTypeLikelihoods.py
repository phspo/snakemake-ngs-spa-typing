import json
import math
from scipy.stats import poisson
import sys
import logging
import numpy as np
from mpmath import *
from probabilistic import *


mp.dps=snakemake.params['dps']
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")

#Read counts
expected_counts = json.load(open(snakemake.input['expected'],'r'))
actual_counts = json.load(open(snakemake.input['observed'],'r'))

#Create Set of valid kmers

validKmers = set(k  for spaTypeID in expected_counts for k in expected_counts[spaTypeID])

#Parse k
k = int(snakemake.params['k'])

#prior = 0

#with open(snakemake.input['prior'],'r') as priorFile:
#	prior = mpf(priorFile.read())

kmerError = -1
with open(snakemake.input['kmerError'],'r') as infile:
	kmerError = float(infile.read().split('\t')[1])
	
probabilities = {}
unexpectedLikelihoods = {}

for spaTypeID,expectedCounts in expected_counts.items():
	
	#print('Calculating Likelihood for SpaType: {}'.format(spaTypeID))
	
	if(len(expectedCounts) == 0):
		logging.warning("No Kmers found for spaType : {}, skipping ...".format(spaTypeID))
		continue
	
	tp,unexpectedLikelihood = calculateKmerLikelihood(actual_counts,expectedCounts,kmerError,k,validKmers)
	
	if tp >= 0:
		raise ValueError("Probability has reached a value >= 0! This doesn't make sense ... \n Probability {} Prior {}".format(probability,prior))
			
	probabilities[spaTypeID] = float(tp)
	unexpectedLikelihoods[spaTypeID] = float(unexpectedLikelihood)

with open(snakemake.output['likelihoods'],'w') as outfile:
	json.dump(probabilities,outfile)

with open(snakemake.output['unexpectedLikelihoods'],'w') as outfile:
	json.dump(unexpectedLikelihoods,outfile)


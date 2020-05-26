import json
import math
from scipy.stats import poisson
import sys
import logging
import numpy as np
from mpmath import *
from probabilistic import *
from multiprocessing import Pool, TimeoutError, Manager
from calcGenerativeTask import calculationTask

mp.dps=snakemake.params['dps']
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")

#Read counts
sequenceKmerProfiles = json.load(open(snakemake.input['counts'],'r'))
actual_counts = json.load(open(snakemake.input['observed'],'r'))

#Parse k
k = int(snakemake.params['k'])

baseErrorRate = 0.0
with open(snakemake.input['baseError'],'r') as infile:
	baseErrorRate = float(infile.read())

#Thread management
pool = Pool(processes=snakemake.threads)
manager = Manager()
print("Calculating Generative Probabilities using {} Threads".format(snakemake.threads))

probabilities = manager.dict()
#probabilities = {}

#HD-Likelihoods LN-Likelihoods, Caches
hdLikelihoods = manager.dict()
lnLikelihoods = manager.dict()

#hdLikelihoods = {}
#lnLikelihoods = {}


pool.map(calculationTask,[[spaTypeID,sequenceKmerProfiles,hdLikelihoods,lnLikelihoods,probabilities,actual_counts,baseErrorRate,k] for spaTypeID,sequenceKmerProfiles in sequenceKmerProfiles.items()]) 

pool.close() # Wait for all threads to finish ...

probabilities = dict(probabilities) #convert back to regular dict that can be written to json	

with open(snakemake.output['likelihoods'],'w') as outfile:
	json.dump(probabilities,outfile)


import json
import math
from scipy.stats import poisson
import sys
import logging

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import MaxNLocator

logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")

#Read counts
expected_counts = json.load(open(snakemake.input['expectedCounts'],'r'))
actual_counts = json.load(open(snakemake.input['observedCounts'],'r'))

predictedType = '???'
with open(snakemake.input['probabilities'],'r') as predictionFile:
	predictedType = predictionFile.read().splitlines()[0].split('\t')[0]

kmer_error = -1
with open(snakemake.input['kmerError'],'r') as errorFile:
	kmer_error = float(errorFile.read().splitlines()[0].split('\t')[1])

groundTruthType = snakemake.params['gt']

gtt_counts = expected_counts[groundTruthType]
prd_counts = expected_counts[predictedType]

gtt_correct = [actual_counts[x] for x in actual_counts if x in gtt_counts ]
gtt_error = [actual_counts[x] for x in actual_counts if not x in gtt_counts ]

prd_correct = [actual_counts[x]  for x in actual_counts if x in prd_counts]
prd_error = [actual_counts[x] for x in actual_counts if not x in prd_counts ]

fig, axs = plt.subplots(2, 2, figsize=(8, 8),sharex='col')

labels, counts = np.unique(gtt_correct, return_counts=True)
axs[0,0].bar(labels, counts, align='center')
axs[0,0].set_title("GTT Correct Kmers")

labels, counts = np.unique(gtt_error, return_counts=True)
axs[0,1].bar(labels, counts, align='center')
axs[0,1].set_title("GTT Error Kmers")

labels, counts = np.unique(prd_correct, return_counts=True)
axs[1,0].bar(labels, counts, align='center')
axs[1,0].set_title("PT Correct Kmers")

labels, counts = np.unique(prd_error, return_counts=True)
axs[1,1].bar(labels, counts, align='center')
axs[1,1].set_title("PT Error Kmers")

plt.savefig(snakemake.output["errors"])

plt.clf()

#Deviation histogram

deviationsGttError = []
deviationsGttActual = []
deviationsGttSum = []
deviationsPrdError = []
deviationsPrdActual = []
deviationsPrdSum = []

epsilonGtt = sum(actual_counts[x] for x in actual_counts)*kmer_error/sum(1 for kmer in actual_counts if not kmer in gtt_counts)
epsilonPrd = sum(actual_counts[x] for x in actual_counts)*kmer_error/sum(1 for kmer in actual_counts if not kmer in prd_counts)

for kmer in actual_counts:
	if kmer in gtt_counts:
		deviationsGttActual.append(actual_counts[kmer]-gtt_counts[kmer])
	else:
		deviationsGttError.append(actual_counts[kmer]-epsilonGtt)
	deviationsGttSum.append(actual_counts[kmer]-epsilonGtt if kmer not in gtt_counts else actual_counts[kmer]-gtt_counts[kmer])
	if kmer in prd_counts:
		deviationsPrdActual.append(actual_counts[kmer]-prd_counts[kmer])
	else:
		deviationsPrdError.append(actual_counts[kmer]-epsilonPrd)
	deviationsPrdSum.append(actual_counts[kmer]-epsilonPrd if kmer not in prd_counts else actual_counts[kmer]-prd_counts[kmer])


fig, axs = plt.subplots(2, 3, figsize=(12, 8),sharex='col')

labels, counts = np.unique(deviationsGttError, return_counts=True)
axs[0,0].bar(labels, counts, align='center')
axs[0,0].set_title("Errors GTT, Epsilon GTT={}".format(epsilonGtt))
labels, counts = np.unique(deviationsGttActual, return_counts=True)
axs[0,1].bar(labels, counts, align='center')
axs[0,1].set_title("Actual GTT")
labels, counts = np.unique(deviationsGttSum, return_counts=True)
axs[0,2].bar(labels, counts, align='center')
axs[0,2].set_title("Sum GTT")

labels, counts = np.unique(deviationsPrdError, return_counts=True)
axs[1,0].bar(labels, counts, align='center')
axs[1,0].set_title("Errors PRD, Epsilon GTT={}".format(epsilonGtt))
labels, counts = np.unique(deviationsPrdActual, return_counts=True)
axs[1,1].bar(labels, counts, align='center')
axs[1,1].set_title("Actual PRD")
labels, counts = np.unique(deviationsPrdSum, return_counts=True)
axs[1,2].bar(labels, counts, align='center')
axs[1,2].set_title("Sum PRD")

plt.savefig(snakemake.output['deviations'])

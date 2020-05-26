import json
import math
from scipy.stats import poisson
import sys
import logging


import numpy as np


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

with open(snakemake.output['likelihoodAnalysis'],'w') as outfile:

	countsPrd = expected_counts[predictedType]
	countsGtt = expected_counts[groundTruthType]

	kmersPrd = set(countsPrd.keys())
	kmersGtt = set(countsGtt.keys())
	
	kmersUnion = kmersPrd.union(kmersGtt)	
	kmersIntersection = kmersPrd.intersection(kmersGtt)

	kmersPrdExclusive = kmersUnion-kmersGtt
	kmersGttExclusive = kmersUnion-kmersPrd

	outfile.write("Exclusive Ground Truth Type Kmers:\n")
	for kmer in kmersGttExclusive:
		outfile.write(kmer+'\n')

	outfile.write("Exclusive Predicted Type Kmers:\n")
	for kmer in kmersPrdExclusive:
		outfile.write(kmer+'\n')

	outfile.write("Mutual kmers with count differences:\n")
	for kmer in kmersIntersection:
		if countsPrd[kmer] != countsGtt[kmer]:
			outfile.write('{}\tGTT:{}\tPRD:{}\tOBS:{}\n'.format(kmer,countsGtt[kmer],countsPrd[kmer],actual_counts[kmer]))


	epsilonGtt = sum(actual_counts[x] for x in actual_counts)*kmer_error/sum(1 for kmer in actual_counts if not kmer in countsGtt)
	epsilonPrd = sum(actual_counts[x] for x in actual_counts)*kmer_error/sum(1 for kmer in actual_counts if not kmer in countsPrd)


	outfile.write("Epsilon GTT: {}\n".format(epsilonGtt))
	outfile.write("Epsilon PRD: {}\n".format(epsilonPrd))


	exclusivePrdLikelihood = 0
	exclusiveGttLikelihood = 0
	intersectionPrdLikelihood = 0
	intersectionGttLikelihood = 0
	errorPrdLikelihood = 0
	errorGttLikelihood = 0

	for kmer in actual_counts:
		if not kmer in kmersUnion:
			pass
		else:
			if kmer in kmersPrdExclusive:
				expectedCount = countsPrd[kmer]
				observedCount = actual_counts[kmer]
				exclusivePrdLikelihood += poisson.logpmf(observedCount,expectedCount)
				expectedCount = epsilonGtt
				observedCount = actual_counts[kmer]
				errorGttLikelihood += poisson.logpmf(observedCount,expectedCount)
			elif kmer in kmersGttExclusive:
				expectedCount = countsGtt[kmer]
				observedCount = actual_counts[kmer]
				exclusiveGttLikelihood += poisson.logpmf(observedCount,expectedCount)
				expectedCount = epsilonPrd
				observedCount = actual_counts[kmer]
				errorPrdLikelihood += poisson.logpmf(observedCount,expectedCount)
			else:
				expectedCount = countsGtt[kmer]
				observedCount = actual_counts[kmer]
				intersectionGttLikelihood += poisson.logpmf(observedCount,expectedCount)
				expectedCount = countsPrd[kmer]
				observedCount = actual_counts[kmer]
				intersectionPrdLikelihood += poisson.logpmf(observedCount,expectedCount)


	outfile.write("Likelihood of exclusive PRD kmers: {}\n".format(exclusivePrdLikelihood))	
	outfile.write("Likelihood of exclusive GTT kmers: {}\n".format(exclusiveGttLikelihood))								
	outfile.write("Likelihood of exclusive GTT kmers in PRD: {}\n".format(errorPrdLikelihood))								
	outfile.write("Likelihood of exclusive PRD kmers in GTT: {}\n".format(errorGttLikelihood))								
	outfile.write("Likelihood of intersection kmers in PRD: {}\n".format(intersectionPrdLikelihood))								
	outfile.write("Likelihood of intersection kmers in GTT: {}\n".format(intersectionGttLikelihood))								
							
			

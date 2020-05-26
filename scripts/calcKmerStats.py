import json
import math
from scipy.stats import poisson
import sys
import logging
from decimal import *
from shared import *

expectedCounts = json.load(open(snakemake.input['expected'],'r'))
observedCounts = json.load(open(snakemake.input['observed'],'r'))
kmer_coverage_estimate = 0
kmerErrorRate = 0.0
k = int(snakemake.params['k'])

with open(snakemake.input['error'],'r') as infile:
	kmerErrorRate = float(infile.read())
	
with open(snakemake.input['coverage_estimate'],'r') as infile:
	lines = infile.read().splitlines()
	kmer_coverage_estimate = extractCoverageEstimate(lines,snakemake.config)

with open(snakemake.output['stats'],'w') as outfile:

	outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
		"SPA TYPE",
		"Expected Different Kmers",
		"Observed Different Kmers",
		"Ratio of Expected Different Kmers that were observed",
		"Total Expected Kmer Expected Count",
		"Total Observed Kmer Count",
		"Unexpected Different Kmers",
		"Total Unexpected Kmer Count",
		"True Kmer Coverage Estimate",
		"Expected Total Error Count"
		)
	)
	
	#
	size_of_observed_kmers = len(observedCounts)
	#
	total_kmer_count_observed = sum(observedCounts.values())
	
	for spaType in expectedCounts:
		#
		size_of_expected_kmers = len(expectedCounts[spaType])
		if (size_of_expected_kmers == 0):
			continue
		#
		of_which_we_observed = 0
		for kmer in expectedCounts[spaType]:
			if kmer in observedCounts:
				of_which_we_observed += 1
		ratio = of_which_we_observed / size_of_expected_kmers
		#
		total_kmer_count_expected = sum(expectedCounts[spaType].values())
		#
		unexpected = [kmer for kmer in observedCounts if not kmer in expectedCounts[spaType]]
		size_of_unexpected = len(unexpected)
		#
		total_kmer_count_observed_not_expected = sum(observedCounts[kmer] for kmer in unexpected)
		#
		expected_error = total_kmer_count_expected * (1/(1-kmerErrorRate)) -total_kmer_count_expected
		#
		outfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
			spaType,
			size_of_expected_kmers,
			size_of_observed_kmers,
			ratio,
			total_kmer_count_expected,
			total_kmer_count_observed,
			size_of_unexpected,
			total_kmer_count_observed_not_expected,
			kmer_coverage_estimate,
			expected_error)
		)


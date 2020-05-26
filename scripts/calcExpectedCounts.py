import json
from shared import *

counts = json.load(open(snakemake.input['counts'],'r'))

kmer_coverage_estimate = 0


with open(snakemake.input['kmerCoverageEstimate'],'r') as infile, open(snakemake.output[0],'w') as outfile:
	lines = infile.read().splitlines()
	
	kmer_coverage_estimate = extractCoverageEstimate(lines,snakemake.config)
	
	#TODO: Instead of list comprehension use some dict "update" function that might perform better
	expected_counts = {
		spaTypeID	:
		{
			#TODO: Rounding???
			kmer : count*kmer_coverage_estimate
			for kmer,count in spaType.items()
		} 
		for spaTypeID,spaType in counts.items()
	}
	
	json.dump(expected_counts,outfile)

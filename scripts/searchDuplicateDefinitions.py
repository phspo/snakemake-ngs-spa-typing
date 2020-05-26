import json
import re
import itertools
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


seqs = []
with open(snakemake.input['types'],'r') as infile:
	lines = infile.read().splitlines()
	for line in lines:
		split = line.split(',')
		name = split[0]
		value = split[1]
		sptRepeats = value.split('-')
		seqs.append((name,sptRepeats))
	
with open(snakemake.output['out'],'w') as outfile:
	#TODO: Progress Bar
	combs = itertools.combinations(seqs,2)
	for r1,r2 in combs:
		if len(r1[1]) != len(r2[1]): #Sequences are obviously not identical if they contain a different number of repeats
			continue
		duplicate = True
		for b1,b2 in zip(r1[1],r2[1]):
			if b1 != b2:
				duplicate = False
				break
		if duplicate:
			outfile.write('Duplicate: {} = {} \n'.format(r1[0],r2[0]))

import json
import re
import itertools
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import sys

repeats = list(SeqIO.parse(snakemake.input['repeats'],'fasta'))

with open(snakemake.output['out'],'w') as outfile:
	combs = itertools.combinations(repeats,2)
	total = len(repeats)*(len(repeats)-1)/2
	processed = 0
	for r1,r2 in combs:
			if str(r1.seq) == str(r2.seq):
				outfile.write('Duplicate: {} = {} \n'.format(r1.name,r2.name))
			processed += 1
			sys.stdout.write('\r Searching types ... {}%'.format(round(processed/total*100,2)))
			sys.stdout.flush()

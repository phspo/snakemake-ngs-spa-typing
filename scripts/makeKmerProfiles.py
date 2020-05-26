import json
import re
import math
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from shared import *

#Creates k-mer profile from a set of sequence, 

spaSequences = SeqIO.parse(snakemake.input['sequences'],'fasta')
k = int(snakemake.params['k'])

profiles = {}
counts = {}

for sequence in spaSequences:


	#print("creating kmer profile for spa type: " + str(sequence.id))

	kmers = {}

	parseKmers(kmers,sequence.seq,k)
	parseKmers(kmers,sequence.seq.reverse_complement(),k) #also note rev comp kmers

	counts[sequence.id] = kmers
	
	profiles[sequence.id] = normalizeKmers(kmers)

with open(snakemake.output['profiles'],'w') as outfile:
	json.dump(profiles,outfile)
	
with open(snakemake.output['counts'],'w') as outfile:
	json.dump(counts,outfile)

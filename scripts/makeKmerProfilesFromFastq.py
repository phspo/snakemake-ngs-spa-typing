from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from itertools import chain
import json
from shared import *

spaSeqs = json.load(open(snakemake.input['expectedCounts'],'r'))

reads1 = SeqIO.parse(snakemake.input['read1'],'fastq')
reads2 = SeqIO.parse(snakemake.input['read2'],'fastq')

k = int(snakemake.params['k'])


### Step 1: Identify the set of valid k-mers (k-mers that occur in any spa-type)

validKmers = {}

for x in spaSeqs:
    for kmer in spaSeqs[x]:
        validKmers[kmer] = 0


### Step 2: Count valid k-mers

stack = chain(reads1,reads2)

for read in stack:
    sequence = read.seq
    for i in range(len(sequence)+1-k):
        kmer = str(sequence[i:i+k])
        if kmer in validKmers:
            validKmers[kmer] += 1

### Step 3: Normalize

normalizedKmers = normalizeKmers(validKmers)

### Step 4: Generate output

with open(snakemake.output['counts'],'w') as outfile:
	json.dump(validKmers,outfile)

with open(snakemake.output['profile'],'w') as outfile:
	json.dump(normalizedKmers,outfile)

import json
import re
import math

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

import sys
from shared import *


import scipy
import scipy.signal

import pysam

import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")

from scipy.stats import poisson

from collections import Counter

#########################################################################################

k = int(snakemake.params['k'])

kmers = {}

f= open(snakemake.output['debug'],"w+")
f.write('RUN\n')

#Parse Meta Data

meta = {}

with open(snakemake.input['regions'],'r') as infile:
    for l in infile.read().splitlines():
        data = l.split('\t')
        meta[data[0]]=(int(data[1]),int(data[2]))

samfile = pysam.AlignmentFile(snakemake.input['alignment'])

readIDs = {}

for read in samfile.fetch():
	
	readID = read.query_name

	if readID not in readIDs:
		readIDs[readID] = 1
	else:
		logging.critical('{}: Read got processed twice!'.format(readID))

	#print("creating kmer profile for read: " + readID)
	
	#print(bitFlags)
	
	if read.is_secondary or read.is_unmapped:
		logging.info('{}: Secondary (or unmapped) Alignment ... skipping'.format(readID))
		continue
	

	reference_start = read.reference_start
	reference_end = read.reference_end
	#qstart = read.query_alignment_start
	#qend = read.query_alignment_end

	region_start,region_end = meta[read.reference_name]
	
	isReversed  = containsFlag(read.flag,4)

	startCutoff = (region_start - reference_start) if region_start > reference_start else 0
	endCutoff = (reference_end -region_end) if reference_end > region_end else 0
	#startCutoff = region_start - qstart if region_start > qstart else 0
	#endCutoff = qned -region_end if qend > region_end else 0


	#print(region_start,region_end,reference_start,reference_end)
	f.write('Reference: {} (RegionBorders: {}/{})\n'.format(read.reference_name,region_start,region_end))
	f.write('Cutaways [{}:{}]: (Alignment Reference Start/End {}/{})\n'.format(startCutoff,endCutoff,reference_start,reference_end))

	f.write('RawSequence: {}\n'.format(read.query_alignment_sequence))



	sequence = read.query_alignment_sequence

	if endCutoff > 0:
		sequence = sequence[:-endCutoff]
	if startCutoff > 0:
		sequence = sequence[startCutoff:]

	f.write('CutSequence: {}\n'.format(sequence[startCutoff:-endCutoff]))
	
	if isReversed: #Reverse
		sequence = Seq(sequence).reverse_complement()

	parseKmers(kmers,sequence,k) #see shared.py
			
#Normalize Vector

normalized_kmers = normalizeKmers(kmers)


	
	
with open(snakemake.output['profile'],'w') as outfile:
	json.dump(normalized_kmers,outfile)
	
with open(snakemake.output['counts'],'w') as outfile:
	json.dump(kmers,outfile)

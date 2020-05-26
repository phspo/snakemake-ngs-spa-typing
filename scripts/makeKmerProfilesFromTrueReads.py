from Bio import SeqIO
from itertools import chain
from shared import *
import json

k = int(snakemake.params['k'])

#Read the meta data (in order to cut reads appopriately)

regionXBorders = (-1,-1)

with open(snakemake.input['regionXMetaData'],'r') as metaFile:
	data = metaFile.read().split('\t')
	regionXBorders=(int(data[0]),int(data[1]))

regionXStart = regionXBorders[0]
regionXEnd = regionXBorders[1]

#Read the correct reads (did really originate from region x)

filteredReads = SeqIO.parse(snakemake.input['filteredReads'],'fastq')

kmers = {}
kmerOrigins = {}


for read in filteredReads:
	
	flags = {}

	startPosRead1,startPosRead2,r1strand,r2strand = parseDwgsimReadName(read.id)
	readLength = len(read)
	
	relevantSequence = read.seq


	readStart = -1
	readStrand = -1

	if (read.name.endswith('/1')):
		readStart = startPosRead1
		readStrand = r1strand
	elif(read.name.endswith('/2')):
		readStart = startPosRead2
		readStrand = r2strand

	#Detect Read Orientation

	readOrientation = -1

	if startPosRead2 < startPosRead1:
		readOrientation = 1 if read.name.endswith('/1') else 0 #reverse
	else:
		readOrientation = 0 if read.name.endswith('/1') else 1 

	flags['readStart'] = readStart
	flags['readStrand'] = readStrand

	flags['readOrientation'] = readOrientation
	
	readStart -= 1
	readEnd = readStart + readLength
	
	
	if regionXStart > readStart:
		clippingStart = regionXStart-readStart
		if readStrand == 0:
			relevantSequence = relevantSequence[clippingStart:]
		else: # read Strand 1
			relevantSequence = relevantSequence[:-clippingStart]
		flags['trimStart'] = clippingStart

	if readEnd > regionXEnd:
		clippingEnd = readEnd-regionXEnd
		if readStrand == 0:
			relevantSequence = relevantSequence[:-clippingEnd]
		else:
			relevantSequence = relevantSequence[clippingEnd:]
		flags['trimEnd'] = clippingEnd

	
	parseKmers(kmers,relevantSequence,k,kmerOrigins,read.name,flags)


kmerOriginsAnalysis = {}

#Post Processing: Analyze total fraction of reversals and complements
for kmer in kmerOrigins:
			
	total = 0
	rv = 0
	cp = 0
	ts = 0
	te = 0
	rs1 = 0
	ro1 = 0

	for read in kmerOrigins[kmer]: #Make more readable
		flags = read[2]
		total += 1
		if 'complement' in flags:
			cp += 1
		if 'reverse' in flags:
			rv += 1
		if 'trimStart' in flags:
			ts += 1
		if 'trimEnd' in flags:
			te += 1
		if flags['readStrand'] == 1:
			rs1 += 1
		if flags['readOrientation'] == 1:
			ro1 += 1
					
	kmerOriginsAnalysis[kmer] = {
		'total' : total ,
		'reversedKmers' : rv,
		'complementedKmers' : cp,
		'startTrimmed' : ts,
		'endTrimmed' : te,
		'readStrand1' : rs1,
		'readOrientation1' : ro1,
		'origins' : kmerOrigins[kmer]
	}


with open(snakemake.output['counts'],'w') as outfile:
	json.dump(kmers,outfile)

with open(snakemake.output['origins'],'w') as outfile:
	json.dump(kmerOriginsAnalysis,outfile)



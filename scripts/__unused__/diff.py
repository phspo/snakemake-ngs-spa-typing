from Bio import SeqIO
import sys
import json
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

record_dict = SeqIO.to_dict(SeqIO.parse(snakemake.input['spaSequences'], "fasta"),key_function = lambda x : int(x.id[1:]) )
types_dict = {}

#Read types into dict
with open(snakemake.input['types'],'r') as infile:
	spaTypes = infile.read().splitlines()
	
	for spaType in spaTypes:
		split = spaType.split(',')
		sid = int(split[0][1:])
		#if (snakemake.params['thinning']):   ->>> Thinning should be redundant here as this already works with the thinned types
		#	if sid > 999:
		#		continue
		value = split[1]
		sptRepeats =[ str(x) for x in value.split('-')]
		types_dict[sid] = sptRepeats

with open(snakemake.input['blastHits'],'r') as infile,open(snakemake.output[0],'w') as outfile:
	#Retrieve best Blast hit
	spaTypeID = int(infile.read().splitlines()[0].split('\t')[1][1:])
	predicted = int(snakemake.params['predictedType'])
	
	if (spaTypeID == predicted):
		outfile.write('Match')
		sys.exit(0)
	
	##Mismatch Analysis #####
	
	#Sequence Layer
	
	s1 = record_dict[spaTypeID]
	s2 = record_dict[predicted]

	outfile.write('Record for determined spa Type:\n')
	outfile.write(str(s1))
	outfile.write('\n\n')
	
	outfile.write('Record for predicted spa Type:\n')
	outfile.write(str(s2))
	outfile.write('\n\n')
	
	outfile.write('SpaType Sequence Difference (Possible Alignments) \n\n')

	for a in pairwise2.align.globalxx(s1.seq,s2.seq):
		outfile.write(format_alignment(*a))
		outfile.write('\n')	
	outfile.write('\n\n')
	
	##Repeat Layer	
	r1 = types_dict[spaTypeID]
	r2 = types_dict[predicted]
	outfile.write(str(r1))
	outfile.write('\n\n')
	outfile.write(str(r2))
	outfile.write('\n\n')
	
	'''

	outfile.write('SpaType Repeat Difference \n')
	for a in pairwise2.align.globalxx(r1,r2):
		outfile.write(format_alignment(*a))
	'''

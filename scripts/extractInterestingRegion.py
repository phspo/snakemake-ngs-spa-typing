import os.path
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import sys

#Extract best hits and determine region in between

contigID = 'UNDEFINED'
positions = []

#ARLP
with open(snakemake.input['arlpBlast'],'r') as infile:
	lines = infile.read().splitlines()
	if lines:
		bestHit = lines[0].split('\t') #Simply take the first Blast Hit (TODO: Search for fitting matches in a smarter way)
		contigID = bestHit[1]
		positions.append(int(bestHit[8]))
		positions.append(int(bestHit[9]))
	else:
		raise ValueError("BLAST File for ARLP is empty!")
		
#LLP
with open(snakemake.input['llpBlast'],'r') as infile:
	lines = infile.read().splitlines()
	if lines:
		bestHit = lines[0].split('\t') #Simply take the first Blast Hit (TODO: Search for fitting matches in a smarter way)
		if bestHit[1] != contigID:
			raise Exception('Best Matches don\'t share the same contig ID')
		positions.append(int(bestHit[8]))
		positions.append(int(bestHit[9]))
	else:
		raise ValueError("BLAST File for LLP is empty!")
		
positions = sorted(positions)

firstPosition = positions[1]
lastPosition  = positions[2]
		
print("Extracting from Scaffold {} [{}:{}]".format(contigID,firstPosition,lastPosition))		

contigs = SeqIO.parse(snakemake.input['scaffolds']+'/scaffolds.fasta','fasta')

for contig in contigs:
	if contig.id == contigID:
		with open(snakemake.output['out'],'w') as outfile:
			subsequence = contig[firstPosition:lastPosition].reverse_complement()
			outfile.write('> '+'Interesting Region (Hypothetical Protein A)'+'\n')
			outfile.write(str(subsequence.seq))
		sys.exit(0)
		
raise ValueError('Couldn\'t find scaffold ID [{}] in contig file'.format(contigID))
		

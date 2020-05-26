from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import bisect
from Bio.SubsMat import MatrixInfo as matlist
matrix = matlist.blosum62


##Parameters for Alignment

spatypes = SeqIO.parse(snakemake.input['seqs'], "fasta")

region = SeqIO.read(snakemake.input['hypProtA'],'fasta')

best = -float('inf')
bestAlignments = []

for spatype in spatypes:
	
	alignment = pairwise2.align.globaldx(region.seq,spatype.seq,matrix,one_alignment_only=True)
	
	score = alignment[0][2]
	
	if score > best:
		best = score
		bestAlignments = [] #Reset List
		bestAlignments.append((spatype.id,alignment[0])) #Tuple [0] = SpaType ID [1] = Actual Alignment
	elif score == best:
		bestAlignments.append((spatype.id,alignment[0]))
	print('Aligned {} to Region, Score : {}'.format(spatype.id,score))
		

with open(snakemake.output[0],'w') as outfile:
	for bestAl in bestAlignments:
		#outfile.write('Best Alignment found for SpaType with ID: {}\n'.format(bestAlignmentID))
		#outfile.write(bestAlignment[0]+'\n')
		#outfile.write(bestAlignment[1]+'\n')
		#outfile.write(str(bestAlignment[2])+'\n')
		#outfile.write(str(bestAlignment[3])+'\n')
		#outfile.write(str(bestAlignment[4]))
		outfile.write(str(bestAl[0])+'\t')
		outfile.write(str(bestAl[1][2])+'\t')
		outfile.write(str(bestAl[1][3])+'\t')
		outfile.write(str(bestAl[1][4])+'\t\n')

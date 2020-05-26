from Bio import SeqIO
import sys

#Filename of the data set
fid = snakemake.params['fid']

with open(snakemake.input['groundTruthFile'],'r') as gtfile,open(snakemake.output['groundTruthSequence'],'w') as outfile:

	groundTruthID = '?'

	for l in gtfile.read().splitlines():
		data = l.split('\t')
		if data[0] == fid:
			groundTruthID = data[1]
			break
	else:
		print("Can't find {} in ground truth file ... check that an entry exists!".format(fid))
		
		sys.exit(-1)

	proteinASequences = SeqIO.parse(snakemake.input['synProteinAs'],'fasta')

	for seq in proteinASequences:
		if seq.id.startswith(groundTruthID):
			SeqIO.write(seq,outfile,'fasta')
		break
	else:
		print("Can't find {} in syn protein a file ... check that an entry exists!".format(groundTruthID))
		
		sys.exit(-1)

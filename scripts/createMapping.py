from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


def cm(refGen,ctgAlgn,out):

	refGenome = SeqIO.read(refGen, "fasta")
	refGenomeRaw = str(refGenome.seq)
	mappedContigs = '?'*len(refGenomeRaw) #Initialize with Question Marks
	
	with open(ctgAlgn,'r') as alignmentFile:
		contigMatches = alignmentFile.read().splitlines()
		for contigMatch in contigMatches:
			content = contigMatch.split('\t')
			startingPos = int(content[3])


cm(
	snakemake.input['refGen'],
	snakemake.input['ctgAlgn'],
	snakemake.output[0]
)

import os.path
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

INDEX_ORGANISM_ID = snakemake.config['reference_genome_table_index_organism_id']
INDEX_START_POS = snakemake.config['reference_genome_table_index_start_pos']
INDEX_END_POS = snakemake.config['reference_genome_table_index_end_pos']
INDEX_STRAND = snakemake.config['reference_genome_table_index_strand']
INDEX_PROTEIN_NAME = snakemake.config['reference_genome_table_index_protein_id']

def main():
	extractSeq(
		snakemake.input['pt'],
		snakemake.params['pid'],
		snakemake.output['main'],
		snakemake.output['meta'],
		snakemake.input['refg']
	)


def extractSeq(pt,pid,out1,out2,fna):
	#Check if protein id exists and retrieve information from protein table
	with open(pt,'r') as tablefile:
		lines = tablefile.read().splitlines()[1:] #Skip first line -> Contains only headers
		for line in lines:
			#Retrieve the fields that are relevant
			entries = line.split('\t')
			organism = entries[INDEX_ORGANISM_ID]
			startPos = int(entries[INDEX_START_POS])
			endPos = int(entries[INDEX_END_POS])
			strand = entries[INDEX_STRAND]
			name = entries[INDEX_PROTEIN_NAME]
			if name == pid:
				record = SeqIO.read(fna, "fasta")
				subsequence = record.seq[startPos:endPos]
				#If the protein is located on the negative strand we need to reverse complement it
				if strand == '-':
					subsequence = subsequence.reverse_complement()
				with open(out1,'w') as outfile:
					outfile.write('> '+pid+' [Positive Strand Data]'+'\n') #Write ID and newline (as pseudo fasta)
					outfile.write(str(subsequence))
				with open(out2,'w') as outfile:
					outfile.write(str(startPos)+'\t')
					outfile.write(str(endPos)+'\t')
					outfile.write(strand+'\t')
					outfile.write(organism)
				#We found the relevant entry and can abort traversing the table
				return
		print("Did not find the entry: {} in the protein table ...".format(pid))
		
if __name__ == "__main__":
	main()

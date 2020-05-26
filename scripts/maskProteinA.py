import os.path
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

INDEX_ORGANISM_ID = snakemake.config['reference_genome_table_index_organism_id']
INDEX_START_POS = snakemake.config['reference_genome_table_index_start_pos']
INDEX_END_POS = snakemake.config['reference_genome_table_index_end_pos']
INDEX_PROTEIN_NAME = snakemake.config['reference_genome_table_index_protein_id']

def main():
	extractSeq(
		snakemake.input['pt'],
		snakemake.params['pid'],
		snakemake.output['main'],
		snakemake.input['refg']
	)


def extractSeq(pt,pid,out,fna):
	#Check if protein id exists and retrieve information from protein table
	with open(pt,'r') as tablefile:
		lines = tablefile.read().splitlines()[1:] #Skip first line -> Contains only headers
		for line in lines:
			#Retrieve the fields that are relevant
			entries = line.split('\t')
			organism = entries[INDEX_ORGANISM_ID]
			startPos = int(entries[INDEX_START_POS])
			endPos = int(entries[INDEX_END_POS])
			name = entries[INDEX_PROTEIN_NAME]
			if name == pid:
				record = SeqIO.read(fna, "fasta")
				fillerLength = endPos+1-startPos
				subsequencePrior = record.seq[0:startPos]
				subsequenceMasked = Seq('N'*fillerLength,generic_dna)
				subsequencePosterior = record.seq[endPos+1:]
				masked = SeqRecord(subsequencePrior+subsequenceMasked+subsequencePosterior,name='Masked Reference Genome',id='maskref',description='Masked Reference Genome (Protein A replaced with N)')
				SeqIO.write(masked,out,'fasta')
				
				return
		print("Did not find the entry: {} in the protein table ...".format(pid))
		
if __name__ == "__main__":
	main()

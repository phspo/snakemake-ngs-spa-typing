from Bio import SeqIO
import json

k = int(snakemake.params['k'])

spaSeqs = json.load(open(snakemake.input['kmerCounts'],'r'))

validKmers = {}

for x in spaSeqs:
    for kmer in spaSeqs[x]:
        validKmers[kmer] = 0



maskedRefSeq = SeqIO.read(snakemake.input['maskedReference'],'fasta').seq


with open(snakemake.output[0],'w') as outfile:
    for i in range(len(maskedRefSeq)+1-k):
        kmer = str(maskedRefSeq[i:i+k])
        if kmer in validKmers:
            outfile.write('kmer {} found in masked reference at position: {} (+ strand)\n'.format(kmer,i))
    else:
        outfile.write('Sense K-mers are unique to the protein A region and not found elsewhere in the genome')

    maskedRefSeq_rc = maskedRefSeq.reverse_complement()

    for i in range(len(maskedRefSeq)+1-k):
        kmer = str(maskedRefSeq_rc[i:i+k])
        if kmer in validKmers:
            outfile.write('kmer {} found in masked reference at position: {} (- strand)\n'.format(kmer,i))
    else:
        outfile.write('Complement K-mers are unique to the protein A region and not found elsewhere in the genome')


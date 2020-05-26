import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

insertionMarker = '[[[]]]'

outFile = snakemake.output['seqs']
metaFile = snakemake.output['metaInf']

frame = str(
	SeqIO.read(snakemake.input['proteinAFrame'],'fasta').seq
)

seqs = SeqIO.parse(snakemake.input['spaSequences'],'fasta')

#List of Synthetic Protein A Sequences
synthSeqs = []
metaInf = []

#Determine the location of the insertionMarker
globalLeft = frame.find(insertionMarker)

for seq in seqs:
	rawSeq = str(seq.seq)
	sequenceStr = frame.replace(insertionMarker,rawSeq)
	right = len(rawSeq)
	sequence = SeqRecord(Seq(sequenceStr,generic_dna),id=seq.id,name=seq.id+'_synth',description=seq.id+' Synthetic Protein A')
	synthSeqs.append(sequence)
	metaInf.append([seq.id,globalLeft,globalLeft+right])

#Write Output


SeqIO.write(synthSeqs,outFile,'fasta')

with open(metaFile,'w') as outFile:
	for inf in metaInf:
		outFile.write(inf[0]+'\t'+str(inf[1])+'\t'+str(inf[2])+'\n')

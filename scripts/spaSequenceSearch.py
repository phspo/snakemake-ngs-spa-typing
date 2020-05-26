import json
import re
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import regex

scaffolds = SeqIO.parse(snakemake.params['scaffolds'],'fasta')
types = SeqIO.parse(snakemake.params['spaSeqs'],'fasta')

with open(snakemake.output[0],'w') as outfile:
	for t in types:
		pattern = str(t.seq)
		pattern_r = regex.compile('(%s){e<=1}' % pattern)
		pattern_rev = str(t.seq.reverse_complement())
		pattern_rev_r = regex.compile('(%s){e<=1}' % pattern_rev)
		for scaffold in scaffolds:
			scaffoldSeqStr = str(scaffold.seq)
			plussearch = pattern_r.search(scaffoldSeqStr)
			if plussearch != None:
				for m in plussearch:
					outfile.write("(+) match in scaffold: {} for spa-sequence: {}, SUB/INS/DEL: {}/{}/{} \n".format(scaffold,t,m.fuzzy_counts[0],m.fuzzy_counts[1],m.fuzzy_counts[2]))
			minussearch = pattern_rev_r.search(scaffoldSeqStr)
			if minussearch != None:
				for m in minussearch:
					outfile.write("(-) match in scaffold: {} for spa-sequence: {}, SUB/INS/DEL: {}/{}/{} \n".format(scaffold,t,m.fuzzy_counts[0],m.fuzzy_counts[1],m.fuzzy_counts[2]))
	

import json
import re
import math
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq

read_profile = json.load(open(snakemake.input['readProfile'],'r'))
spa_profiles = json.load(open(snakemake.input['spaProfiles'],'r'))

distances = []

with open (snakemake.output[0],'w') as outfile:
	for profileID,profile in spa_profiles.items():
	
		keys_read = set(read_profile.keys())
		keys_spa = set(profile.keys())
		
		#print(keys_read)
		#print(keys_spa)
	
		all_kmers = keys_read.union(keys_spa)
		
		#shared = keys_read.intersection(keys_spa)
		
		#print("{} of {} kmers are shared".format(len(shared),len(all_kmers)))
		
		tmpsum = 0
		for kmer in all_kmers:
			pos1 = profile[kmer] if kmer in profile else 0
			pos2 = read_profile[kmer] if kmer in read_profile else 0
			tmpsum += (pos1-pos2)*(pos1-pos2)
		distance = math.sqrt(tmpsum)
		
		distances.append((profileID,distance))
		
	distances = sorted(distances,key=lambda x: x[1])
	
	for d in distances:
		outfile.write(str(d[0])+'\t'+str(d[1])+'\n')

import pysam


samfile = pysam.AlignmentFile(snakemake.input['alignment'])

counts = 0
plus = 0
minus = 0

for read in samfile.fetch():
	
	if not (read.is_unmapped or read.is_secondary):
		
		counts += 1
		
		if read.is_reverse:
			minus += 1
		else:
			plus += 1

with open(snakemake.output['strandbias'],'w') as outfile:
	outfile.write('{}\t{}\t{}'.format(counts,plus,minus))

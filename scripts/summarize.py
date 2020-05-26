import os

def getReadIDFromPath(path):
	return os.path.basename(os.path.dirname(path))

def summarize(input,output):
	with open(output,'w') as outfile:
		for f in input:
			with open(f,'r') as infile:
				data = infile.read().splitlines()
				bestHit = data[0].split('\t')
				outfile.write(getReadIDFromPath(f)+'\t'+bestHit[0]+'\t'+bestHit[1]+'\n')

summarize(snakemake.input['results'],snakemake.output['summary'])

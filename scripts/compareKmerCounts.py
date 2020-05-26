import json

trueCounts = json.load(open(snakemake.input['trueCounts'],'r'))
observedCounts = json.load(open(snakemake.input['observedCounts'],'r'))

keys = set().union(trueCounts.keys(),observedCounts.keys())


with open(snakemake.output['differences'],'w') as outfile:
	for x in keys:
		outfile.write(x+'\t')
		if x in observedCounts and x in trueCounts:
			outfile.write(str(observedCounts[x]-trueCounts[x]))
		elif x in observedCounts:
			outfile.write(str(observedCounts[x]))
		elif x in trueCounts:
			outfile.write(str(0-trueCounts[x]))
		outfile.write('\n')

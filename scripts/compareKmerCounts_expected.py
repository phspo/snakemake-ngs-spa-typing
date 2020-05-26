import json

spaType = 'unknown'

with open(snakemake.input['groundTruthFile'],'r') as infile:
	data = infile.read().splitlines()
	for line in data:
		lineData = line.split('\t')
		if lineData[0] == snakemake.params['inputFileID']:
			spaType = lineData[1]
			break

trueCounts = json.load(open(snakemake.input['trueCounts'],'r'))
expCounts = json.load(open(snakemake.input['expectedCounts'],'r'))[spaType]
diff = {}

keys = expCounts.keys()


for x in keys:
	if x in expCounts and x in trueCounts:
		diff[x]=(expCounts[x],trueCounts[x])
	elif x in expCounts:
		diff[x]=(expCounts[x],0)
	elif x in trueCounts:
		diff[x]=(0,trueCounts[x])

with open(snakemake.output['differences'],'w') as outfile:

	
	totalDiffShared = 0
	totalKmerCount = 0

	outfile.write('[ Shared Kmers ] \n')
	keys = set.intersection(set(expCounts.keys()),set(trueCounts.keys()))
	for x in sorted(keys):
		outfile.write('{} -> Expected: {} TrueReads: {} Ratio : {}\n'.format(x,expCounts[x],trueCounts[x],round(expCounts[x]/trueCounts[x],2)))
		totalDiffShared += (trueCounts[x]-expCounts[x])
		totalKmerCount += trueCounts[x]
	totalDiffUnobserved = 0

	outfile.write('[ Unobserved Kmers ] \n')
	keys = set(expCounts.keys()).difference(set(trueCounts.keys()))
	for x in sorted(keys):
		outfile.write('{} -> Expected: {}  \n'.format(x,expCounts[x]))
		totalDiffUnobserved -= expCounts[x]

	totalDiffUnexpected = 0
	outfile.write('[ Unexpected Kmers ] \n')
	keys = set(trueCounts.keys()).difference(set(expCounts.keys()))
	for x in sorted(keys):
		outfile.write('{} -> TrueReads: {}  \n'.format(x,trueCounts[x]))
		totalDiffUnexpected += trueCounts[x]
		totalKmerCount += trueCounts[x]

	outfile.write('[ Stats ] \n')
	outfile.write('Total Amount of Observed Kmers: {}\n'.format(totalKmerCount))
	outfile.write('Total Diff Shared Kmers: {}\n'.format(totalDiffShared))
	outfile.write('Total Diff Unobserved Kmers: {}\n'.format(totalDiffUnobserved))
	outfile.write('Total Diff Unexpected Kmers: {}\n'.format(totalDiffUnexpected))

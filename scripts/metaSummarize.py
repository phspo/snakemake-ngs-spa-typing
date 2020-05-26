import sys

import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")

groundTruth = {}

with open(snakemake.input['groundTruth'],'r') as groundTruthFile:
	lines = groundTruthFile.read().splitlines()
	for l in lines:
		d = l.split('\t')
		groundTruth[d[0]] = d[1]


def metaSummarize(input,output):
	with open(output,'w') as outFile:
		for f in input:
			hitCount = 0
			with open(f,'r') as inFile:
				data = inFile.read().splitlines()
				for line in data:
					readID = line.split('\t')[0]
					predictedTypeID = line.split('\t')[1]
					if not readID in groundTruth:
						logging.critical('ReadID: {} is not part of the ground truth provided! Terminating ...'.format(readID))
						print(groundTruth)
						sys.exit(-1)
					hitCount += 1 if  predictedTypeID == groundTruth[readID] else 0
			outFile.write(f+'\t'+str(hitCount)+'\n')


metaSummarize(snakemake.input['summary'],snakemake.output['meta'])

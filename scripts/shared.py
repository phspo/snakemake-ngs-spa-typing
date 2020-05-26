import math

def extractCoverageEstimateFile(file,config):
	with open(file,'r') as infile:
		lines = infile.read().splitlines()
		return extractCoverageEstimate(lines,config)


def extractCoverageEstimate(lines,config):
	if config['kmerCoverageEstimationMethod'] == 'alignment':
		return round(float(lines[0]))
	elif config['kmerCoverageEstimationMethod'] == 'countMean':
		return  round(float(lines[0].split('\t')[1]))
	else:
		return int(lines[1].split('\t')[1])


def parseDwgsimReadName(name):
	data = name.split('_')
	#print(data)
	startPosRead1 = int(data[1])
	startPosRead2 = int(data[2])

	strandRead1 = int(data[3])
	strandRead2 = int(data[4])
	return (startPosRead1,startPosRead2,strandRead1,strandRead2)


def parseWgsimReadName(name):
	data = name.split('_')
	#print(data)
	startPosRead1 = int(data[1])
	startPosRead2 = int(data[2])
	return (startPosRead1,startPosRead2)



def parseKmers(kmers,sequence,k,origins=None,readSource=None,flags=None):
	if len(sequence) < k:
		return
	for i in range(len(sequence)+1-k):
		kmer = str(sequence[i:i+k])
		if kmer in kmers:
			kmers[kmer] += 1
		else:
			kmers[kmer] = 1
		
		if readSource:
			if kmer in origins:
				origins[kmer].append((readSource,i,flags))
			else:
				origins[kmer] = [(readSource,i,flags)]

def containsFlag(integer,flagNumber):
	return (True if (integer & (1 << flagNumber)) else False)

#Normalizes a dictionary of k-mer counts (results in the length of the k-mer vector being 1)
def normalizeKmers(kmers):
	normalized_kmers = {}

	length = 0

	for kmer in kmers:
		length += kmers[kmer]*kmers[kmer]

	length = math.sqrt(length)

	for kmer in kmers:
		normalized_kmers[kmer] = kmers[kmer] / length

	return normalized_kmers


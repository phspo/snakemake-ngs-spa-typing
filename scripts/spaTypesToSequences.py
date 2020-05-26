from Bio import SeqIO

#Load repeat sequences in fasta format
repeats = SeqIO.parse(snakemake.input['repeats'],'fasta')

repeatsDict = {}

filterList = []

#Pre-Process Repeats

if 'filterList' in snakemake.input:
	with open(snakemake.input['filterList'],'r') as infile:
		filterList = infile.read().split(',')

for repeat in repeats:
	repeatsDict[repeat.id[1:]]=str(repeat.seq)

with open(snakemake.input['types'],'r') as infile, open(snakemake.output['out'],'w') as outfile:
	spaTypes = infile.read().splitlines()
	spaTypesExtended = []
	
	for spaType in spaTypes:
		split = spaType.split(',')
		name = split[0]
		if len(filterList) > 0:
			if not name in filterList:
				continue
		value = split[1]
		#print(name,value)
		sptRepeats = value.split('-')
		if not sptRepeats:
			continue
		sequence = ''
		for repeat in sptRepeats:
			sequence += repeatsDict[repeat]
		outfile.write('>'+name+'\n'+sequence+'\n')

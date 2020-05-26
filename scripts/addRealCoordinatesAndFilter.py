import sys

def filterAddRealCoords(
	repeatMatchFilePath,
	contigMatchFilePath,
	protAMetaData,
	outFilePath
	):
	#Create Dictionary of Contigs and Coordinates
	contigs = {}
	with open(contigMatchFilePath,'r') as infile:
		contigHits = infile.read().splitlines()[2:] #First two lines contain metadata
		for hit in contigHits:
			entries = hit.split('\t')
			contigs[entries[0]]=int(entries[3])
	#Extract Protein A start and end pos
	startPos = -1
	endPos = -1
	
	with open(protAMetaData,'r') as infile:
		data = infile.read().split('\t')
		startPos = int(data[0])
		endPos = int(data[1])
		
	#Write filtered matches		
	with open(repeatMatchFilePath,'r') as infile, open(outFilePath,'w') as outfile:
		blastHits = infile.read().splitlines()
		for hit in blastHits:
			entries = hit.split('\t')
			
			offset = contigs[entries[1]]
			print(entries[1]," Offset:",offset)
			hit_startPos  = int(entries[8])+offset
			hit_endPos = int(entries[9])+offset
			
			print(startPos,hit_startPos,hit_endPos,endPos)
			
			if (startPos <= hit_startPos) and (hit_endPos <= endPos): #Check if Match occurs in protein A region  
				outfile.write(
						entries[0]+'\t', #repeatID
						entries[2]+'\t', #percentageMatching
						entries[10]+'\t', #evalue
						####### Coordinates
						hit_startPos+'\t', #start
						hit_endPos+'\n' #end
				)
				
filterAddRealCoords(
	snakemake.input['repm'],
	snakemake.input['ctgm'],
	snakemake.input['pam'],
	snakemake.output[0]
)

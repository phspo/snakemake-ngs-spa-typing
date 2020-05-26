import sys


#Step 1: Read the protein table 

INDEX_START_POS = snakemake.params['idxstart']
INDEX_END_POS = snakemake.params['idxend']
INDEX_STRAND = snakemake.params['idxstrand']
INDEX_PROTEIN_ID = snakemake.params['idxprotein']

proteinAIdentifier = snakemake.params['idProtA']

spaType = snakemake.params['spaType']

startPos,endPos,strand = -1,-1,-1

#Check if protein id exists and retrieve information from protein table
with open(snakemake.input['proteintable'],'r') as tablefile:
	lines = tablefile.read().splitlines()[1:] #Skip first line -> Contains only headers
	for line in lines:
		#Retrieve the fields that are relevant
		entries = line.split('\t')
		startPos = int(entries[INDEX_START_POS])
		endPos = int(entries[INDEX_END_POS])
		strand = entries[INDEX_STRAND]
		name = entries[INDEX_PROTEIN_ID]
		if name == pid:
			print(name,startPos,endPos,strand)
			break
	else:
		sys.exit(-1)
		print("Did not find the entry: {} in the protein table ...".format(pid))

#Step 2: Read the reference genome
sequence = SeqIO.read(snakemake.input['referenceGenome'],'fasta')
seqPre = sequence[:startPos]
seqPost = sequence[endPos:]

#Step 3: Read synthetic protein A

proteinAs = SeqIO.to_dict(SeqIO.parse(snakemake.input['syntheticProteinAsSequences'],'fasta'))
synthProteinA = proteinAs[spaType]

#Read meta information
regionXStart = -1
regionXEnd = -1
with open(snakemake.input['syntheticProteinAsMetaFile'],'r') as metafile:
    lines = metafile.read().splitlines()
    lines = [line.split('\t') for line in lines]
    metaInf = {
        line[0] : (int(line[1]),int(line[2])) for line in lines 
    }
    regionXStart,regionXEnd = metaInf[spaType]

#print("Reine Sequenz")
#print(synthProteinA[regionXStart:regionXEnd].seq)

#print(synthProteinA.seq)
#Calculate offset and get total coordinates for region X

regionXLength = regionXEnd - regionXStart
proteinALength = len(synthProteinA)
#print(regionXLength,proteinALength)

paddingStart = proteinALength - regionXEnd
paddingEnd = regionXStart


artificalReference = seqPre+synthProteinA.reverse_complement()+seqPost
artificalReference.id = spaType
artificalReference.description = "Artifical MRSA Genome"
genomeLength = len(artificalReference)

#print("Eingebaute Sequenz")
#print(artificalReference[paddingStart+startPos:paddingStart+startPos+regionXLength].reverse_complement().seq)
#print(artificalReference[startPos:startPos+proteinALength].reverse_complement().seq)

with open(snakemake.output['syntheticReference'],'w') as outfile:
    SeqIO.write(artificalReference,outfile,'fasta')
    
with open(snakemake.output['syntheticReferenceMeta'],'w') as outfile:
    outfile.write(str(paddingStart+startPos)+'\t'+str(paddingStart+startPos+regionXLength))


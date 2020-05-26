k = int(snakemake.params['k'])

with open(snakemake.input['coverageEstimate'],'r') as covfile, open(snakemake.input['readLengthEstimate'],'r') as rlfile,open(snakemake.input['baseErrorEstimate'],'r') as errfile, open(snakemake.output['kmerCoverage'],'w') as outfile:
	coverageEstimate = round(float(covfile.read()))
	readLengthEstimate = round(float(rlfile.read()))
	errorEstimate = float(errfile.read())
	outfile.write(
		str(
			1/2*
			coverageEstimate * (readLengthEstimate-k+1) / readLengthEstimate *(
			1-(errorEstimate*k)	#error term
)
		
		)
	)

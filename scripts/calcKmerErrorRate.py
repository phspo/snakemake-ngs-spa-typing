
#########################################################################################

k = int(snakemake.params['k'])

expectedBaseErrorRate  = -1
with open(snakemake.input['baseError'],'r') as infile:
	expectedBaseErrorRate = float(infile.read())

kmerSurvivalP = (1-expectedBaseErrorRate)**k
error_rate_estimate = 1 - kmerSurvivalP

with open(snakemake.output['error'],'w') as outfile:
	outfile.write(str(error_rate_estimate))


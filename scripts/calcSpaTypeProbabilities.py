import json
from mpmath import *
from probabilistic import *

with open(snakemake.input['likelihoods'],'r') as likelihoods, open(snakemake.input['prior'],'r') as priorfile:
    likelihoods = json.load(likelihoods)
    prior = mpf(priorfile.read())
    probabilities = {}

    for spaType in likelihoods:
        if likelihoods[spaType] == None: #or likelihoods[spaType] < snakemake.config['likelihoodCutoff']:
            probabilities[spaType] = None
        else:
            probabilities[spaType] = likelihoods[spaType]-prior

    sorted_probs = sorted(
        [(k,v) for k,v in probabilities.items()],
        key = lambda x : -float('inf') if x[1] is None else x[1],
        reverse = True
    )

    with open(snakemake.output['probabilities'],'w') as outfile:
        for p in sorted_probs:
            outfile.write(str(p[0])+'\t'+nstr(p[1],20)+'\n')

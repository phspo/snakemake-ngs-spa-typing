import json
import math
import logging
from probabilistic import *
from mpmath import *

logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")

mp.dps=snakemake.params['dps']


#p_lambda(k) = ( (lambda^k)/(k!) ) * e ^ (- lambda)

with open(snakemake.input['likelihoods'],'r') as likelihoodFile:
    likelihoods = json.load(likelihoodFile)

    prior = None
    total = 0

    for spaTypeID in likelihoods:
        likelihood = likelihoods[spaTypeID]
        if likelihood is None: #Ignore ... well ignored types
            pass
        #elif likelihood < snakemake.config['likelihoodCutoff']:
        #    pass
        else:
            total += 1
            try:
                prior = logSpaceAdd(likelihood,prior)
                #print(prior)
            except Exception as e:
                logging.info(e)
                raise e


    prior -= log(1/total)

    logging.info('Finished calculating prior: {}'.format(prior))

    logging.info('Saving to file ...')

    #finally sort the list and save to file ...
    with open(snakemake.output['priorFilePath'],'w') as priorFile:
        priorFile.write(str(prior))

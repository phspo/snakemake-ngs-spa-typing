from probabilistic import *


def calculationTask(args):

	#this is some unelegant code but apparently multiprocessing library has limitations that force this
	
	spaTypeID = args[0]
	sequenceKmerProfiles = args[1]
	hdLikelihoods = args[2]
	lnLikelihoods = args[3]
	probabilities = args[4]
	actual_counts = args[5]
	baseErrorRate = args[6]
	k = args[7]

	if(len(sequenceKmerProfiles) == 0):
		logging.warning("No Kmers found for spaType : {}, skipping ...".format(spaTypeID))
		return -1

	print('Calculating Likelihood for SpaType: {}'.format(spaTypeID)) # Remove later for speedup

	
	tp = calculateKmerLikelihood_Generative(actual_counts,sequenceKmerProfiles,baseErrorRate,k,hdLikelihoods,lnLikelihoods)
	
	#Switch to logspace here!!!
	tp = memolog(tp)

	if tp >= 0:
		raise ValueError("Probability has reached a value >= 0! This doesn't make sense ... \n SpaTypeID {} Probability {}".format(spaTypeID,tp))
	
	#TODO: Conversion to float at this point?
	probabilities[spaTypeID] = float(tp)
	
	#print('Calculated a Log-Likelihood of : {}'.format(tp))

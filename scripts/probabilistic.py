import mpmath as mp

import logging

#define memoize function to improve performance TODO: Evaluate if this is indeed a speedup in terms of performance
memopower = mp.memoize(mp.power)
memolog1p = mp.memoize(mp.log1p)
memoexp = mp.memoize(mp.exp)
memolog = mp.memoize(mp.log)

'''

def calculateKmerLikelihood(observedCounts,expectedCounts,kmerError,k,validKmers):


	#We use a mpmath float to have a sufficiently precise numerical representation
	probability = mp.mpf(0)

	#Split probability for analysis purposes
	unexpectedProbability = mp.mpf(0)
	
	#Calculate default value for expected counts
	sharedKmers = set(observedCounts.keys()).difference(set(expectedCounts.keys()))
	expectedDefaultValue = kmerError*sum(observedCounts.values())/len(sharedKmers)
	#print("Expected Count for Zeros:"+str(expectedDefaultValue))

	for kmerID in validKmers:

		#If we didn't observe a kmer at all we assume we observed it 0 times
		observedCount = 0

		expectedCount = expectedDefaultValue

		#Otherwise we have counted before how many times we observed it
		if kmerID in observedCounts:
			observedCount = observedCounts[kmerID]
			
		if kmerID in expectedCounts:
			expectedCount = expectedCounts[kmerID]

		#Use log10 to calculate probabilities
		local_p = poisson.logpmf(observedCount,expectedCount)
		#local_p = mp.log(mp.exp(-expectedCount)*(mp.mpf(expectedCount)**mp.mpf(observedCount))/(mp.fac(observedCount)))

		#print(local_p)
		
		if (local_p == 0):
			logging.warning('Warning ...! local_p has reached {}'.format(local_p))
			logging.warning('kmerID:{} observed:{} expected: {}'.format(kmerID,observedCount,expectedCount))
			sys.exit(-1)
		else:
			try:
				probability += local_p
				if not (kmerID in expectedCounts):
					unexpectedProbability += local_p

			except ValueError as e:
				logging.warning('kmerID:{} observed:{} expected: {}'.format(kmerID,observedCount,expectedCount))
				logging.warning('calculated (raw) probability of {}\%'.format(local_p))
				logging.critical(str(e))
				sys.exit(-1)
			
		if (probability == 0):
			logging.critical('Warning: Total log probability has reached 0! float insufficient?')
			sys.exit(-1)

	return probability, unexpectedProbability
	
#TODO: Move to shared as this is prooooobably not a very probabilistic thing
def hamming_distance(seq1,seq2):
	#is not called very often due to memoization but maybe there are more efficient ways of doing this (quick bit operation)
	return sum(1 if x != y else 0 for x,y in zip(seq1,seq2))


USE_OWN_MEMOIZATION = True #Uses additional memoization besides mpmath funtionality
MEMOIZE_ON_PROCESS_LEVEL = True #Only uses a cache on a process level and not a shared cache between processes (avoids locks)

#Alternative Model Alex proposed ... let's see how well it does
def calculateKmerLikelihood_Generative(observedCounts,sequenceKmerProfiles,baseErrorRate,k,hdLikelihoods,likelihoods):
	
	#don't use global caches to avoid locks, just cache locally
	if MEMOIZE_ON_PROCESS_LEVEL:
		hdLikelihoods = {}
		likelihoods = {}
	
	#Count how often kmers occur
	nrOfSequenceKmers = sum(sequenceKmerProfiles[kmer] for kmer in sequenceKmerProfiles)

	#Assume uniform distribution for drawing generator kmers
	chanceOfDrawingKmer = 1/nrOfSequenceKmers

	#print("Chance of drawing kmer: {}".format(chanceOfDrawingKmer))
	#We start with a log probability of 0 -> Actual probability is 1 and then factor in

	probability = mp.mpf(1)

	#currentPos = 0
	#total = len(observedCounts)
	
	for kmerIDObs in observedCounts:
		
		#print("Processing kmer {}/{} ...".format(currentPos,total))
		
		#We start with a log probability of None -> Actual probability is 0
		#obsProbability = None
		obsProbability = 0
		
		factorObs = observedCounts[kmerIDObs]
		
		for kmerIDGen in sequenceKmerProfiles:
			#Let's see how many times the kmer is contained in the sequence
			factor = sequenceKmerProfiles[kmerIDGen]
			hd = hamming_distance(kmerIDObs,kmerIDGen)

			calculateHdLikelihood = True
			hdLikelihood = None

			if USE_OWN_MEMOIZATION:
				if hd in hdLikelihoods:
					#Maybe we already calculated this?
					calculateHdLikelihood = False
					hdLikelihood = hdLikelihoods[hd]


			
			if calculateHdLikelihood:
				#print("Calculating Likelihood for hd:{} for the first time ...".format(hd))
				#print((1-baseErrorRate),k,hd,(k-hd))
				intactBasesProbability = memopower((1-baseErrorRate),(k-hd))
				alteredBasesProbability = memopower(baseErrorRate,hd)
				hdLikelihood = intactBasesProbability * alteredBasesProbability
				if USE_OWN_MEMOIZATION:
					hdLikelihoods[hd] =  hdLikelihood
				#print("It is: {} ({}/{})".format(hdLikelihoods[hd],intactBasesProbability,alteredBasesProbability))


			obsProbability += factor*hdLikelihood


		obsProbability *= chanceOfDrawingKmer

		obsProbability = memopower(obsProbability,factorObs)


		if (mp.isnan(obsProbability)):
			raise ValueError("Probability is not a number: {}!".format(obsProbability))



		probability *= obsProbability


	if (mp.isnan(probability)):
		raise ValueError("LogProbability is not a number: {}!".format(probability))
	

	return probability
	

'''

def logSpaceAdd(update,before):
	'''
	Assume a >= b

	We want to calculate ln(exp(ln(a))+exp(ln(b))), thus
	ln(
	exp(ln(b))*(1+exp(ln(a)-ln(b)))
	)
	->
	ln(b) + ln(1+exp(ln(a)-ln(b)))
	->
	ln(b) + ln1p(exp(ln(a)-ln(b)))
	'''
	#As there is no neutral element of addition in log space we only start adding when two values are given
	if before == None:
		return update
	else:
		a = None
		b = None
		#Check which value is larger
		if update > before:
			b = update
			a = before
		else:
			b = before
			a = update

		x = mp.mpf(a)-mp.mpf(b)
		#print(update,before)
		xexp = memoexp(x)
		#print('Exp:',xexp)
		val = memolog1p(xexp)
		#print('Log:',val)
		val = val + b
		if val == 0:
			print('a:{} b:{} x:{} xexp:{} log1p:{}'.format(a,b,x,xexp,val))
			raise ValueError('LogSpace Addition has resulted in a value of 0: a = {} b = {} (possible underflow)'.format(a,b))
		#elif val > 0:
			#print('a:{} b:{} x:{} xexp:{} log1p:{}'.format(a,b,x,xexp,val))
			#raise ValueError('Prior Update has resulted in a value > 0: a = {} b = {} val = {}'.format(a,b,val))
		
		if before == val:
			#print('a:{} b:{} x:{} xexp:{} log1p:{}'.format(a,b,x,xexp,val))
			logging.warning('LogAddition had no effect: a = {} b = {} val = {}'.format(a,b,val))
			#raise ValueError('LogAddition had no effect: a = {} b = {} val = {}'.format(a,b,val))
			#raise ValueError('LogAddition had no effect!')
			
		if mp.isnan(val):
			raise ValueError('LogSpace Addition has resulted in a value of nan: a = {} b = {}'.format(a,b))

		if mp.fabs(val) > 1000: #At this point who cares let's round a bit
			val = mp.floor(val)
		
		return val

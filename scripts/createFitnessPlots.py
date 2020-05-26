import json
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from scipy.stats import poisson

observedKmers = {}
with open(snakemake.input['counts'],'r') as infile:
	observedKmers = json.load(infile)

rawCounts = [x for x in observedKmers.values()]


b=np.amax(rawCounts)

hist,_ = np.histogram(rawCounts,bins=b)

local_mins = argrelextrema(hist,np.less)[0] #choosing local min here, see DOI: 10.1093/bib/bbv029
threshold = local_mins[0]

vals = [x for x in observedKmers.values() if x > threshold]

b = np.arange(np.amax(vals)+1)
b_centers = np.arange(np.max(vals))+0.5
b_reduced = np.arange(np.max(vals))

plt.rcParams['figure.figsize'] = [10, 5]
hist,bin_edges,patches = plt.hist(vals,bins=b)

ratios = json.load(open(snakemake.input['ratios'],'r'))


filteredSpaTypes = []
with open(snakemake.input['probabilities'],'r') as infile:
	lines = infile.read().splitlines()[:3]
	for l in lines:
		filteredSpaTypes.append(l.split('\t')[0])



scores = {}


maxKmerCountX = np.argmax(hist)
maxKmerCountY = np.amax(hist)


for spaType in filteredSpaTypes:

	spaTypeRatios = ratios[spaType]

	def multiPeakPoisson(x,firstPeak,peakHeight):
		return sum( poisson.pmf(x,firstPeak*(i+1))*spaTypeRatios[i]*peakHeight for i in range(len(spaTypeRatios)))

	params,cm = curve_fit(multiPeakPoisson,b_reduced,hist,p0=[maxKmerCountX,maxKmerCountY],bounds=([0,0],[np.amax(vals)+1,1000*maxKmerCountY]),maxfev=5000) #TODO: maxfev exc handling

	x = b_centers
	y = multiPeakPoisson(b_reduced,*params)
	
	fitScore = 0

	for xv,yv in zip(x,y):
		fitScore += abs(yv - xv)

	#print(x,y)

	plt.plot(x,y,label='{}({}/{})'.format(spaType,fitScore,params))

plt.legend()
plt.savefig(snakemake.output[0])

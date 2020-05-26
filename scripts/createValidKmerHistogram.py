import json
import math
from scipy.stats import poisson
import sys
import logging

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import MaxNLocator

logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")

#Read counts
expected_counts = json.load(open(snakemake.input['expected'],'r'))
actual_counts = json.load(open(snakemake.input['observed'],'r'))

#Create Set of valid kmers

validButUnexpectedKmers = set(k  for spaTypeID in expected_counts for k in expected_counts[spaTypeID])-set(k for k in expected_counts[snakemake.params['gtt']])

withZeroCounts = [actual_counts[x] if x in actual_counts else 0 for x in validButUnexpectedKmers]
nonZeroCounts = [actual_counts[x] for x in validButUnexpectedKmers if x in actual_counts ]


fig, axs = plt.subplots(1, 2, figsize=(8, 4))

labels, counts = np.unique(withZeroCounts, return_counts=True)
axs[0].bar(labels, counts, align='center')
axs[0].set_xticks(labels)
axs[0].set_title("Mean={}".format(np.mean(withZeroCounts)))

labels, counts = np.unique(nonZeroCounts, return_counts=True)
axs[1].bar(labels, counts, align='center')
axs[1].set_xticks(labels)
axs[1].set_title("Mean={}".format(np.mean(nonZeroCounts)))

plt.savefig(snakemake.output['histogram'])

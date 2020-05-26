import json
import re
import math
import tarfile

from Bio import SeqIO

from collections import defaultdict

import gzip

import sys

import seaborn as sns, numpy as np
import matplotlib.pyplot as plt

import scipy
import scipy.signal

from itertools import chain

import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,format="%(asctime)s:%(levelname)s:%(message)s")


#########################################################################################



rf1 = snakemake.input['read1']
rf2 = snakemake.input['read2']

reads_1 = SeqIO.parse(rf1,'fastq')
reads_2 = SeqIO.parse(rf2,'fastq')

expectedBaseErrors = 0
totalBases = 0

for read in chain(reads_1,reads_2):

	sequence = read.seq
	phredScores = [10**(-x/10) for x in read.letter_annotations["phred_quality"]]
	
	totalBases += len(sequence)
	
	expectedBaseErrors += sum(phredScores)

expectedBaseErrorRate = expectedBaseErrors/totalBases

with open(snakemake.output['baseError'],'w') as outfile:
	outfile.write(str(expectedBaseErrorRate))

import json
import re
import tempfile
import networkx as nx
import sys

def k_mers(inputFilePath, k):
    with open(inputFilePath,'r') as inputFile:
        #Handle Negative Strand Reads?
        lines = inputFile.read().splitlines()
        rawReads = lines[1::2][::2] #Maybe this is a bit hacky?
        
        ret = []
        
        total = len(rawReads)
        processed = 0
        
        for line in rawReads:
            sys.stdout.write('\r') #carriage return to overwrite progress
            #print(line)
            ret+=[line[x:x+k] for x in range(len(line)-(k-1))]
            sys.stdout.write('Creating k-mers for set of reads ... {}%'.format(round(processed/total*100,2)))
            sys.stdout.flush()
            processed += 1
        return ret

def de_bruijn(kmers):
    g = nx.MultiDiGraph()
    k = len(kmers[0])
    print(k)
    for kmer in kmers:
        left = kmer[0:k-1]
        right = kmer[1:k]
        #print(left,right)
        g.add_edge(left,right)
    return g

km = k_mers(snakemake.input['reads_1'], snakemake.params['k'])+k_mers(snakemake.input['reads_2'], snakemake.params['k'])
graph = de_bruijn(km)
nx.nx_pydot.to_pydot(graph).write_png(snakemake.output['graph'])
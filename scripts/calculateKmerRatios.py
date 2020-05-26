from collections import Counter
import json

ratios = {}
spaSeqs = {}

with open(snakemake.input[0],'r') as infile:
    spaSeqs = json.load(infile)

k = int(snakemake.params['k'])

for spaType in spaSeqs:
    cnt = Counter(spaSeqs[spaType].values())
    total = sum(cnt.values())
    res = []
    for count in sorted(cnt.keys()):
        res.append(cnt[count]/total)
    ratios[spaType] = res

with open(snakemake.output[0],'w') as outfile:
    json.dump(ratios,outfile)

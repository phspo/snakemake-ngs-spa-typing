import seaborn as sns
import matplotlib.pyplot as plt

x = []
y = []

coverageData = {}
counts = {}

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none'
}
plt.rcParams.update(new_rc_params)

for f in snakemake.input:
	#print(f)
	fileid = f.split('/')[-2]
	metafile = snakemake.params['metaFiles'][f]
	data = open(metafile,'r').read().split('\t')
	start = int(data[0])
	end = int(data[1])
	#print(start,end)
	pileup = open(f,'r').read().splitlines()
	for l in pileup:
		data = l.split('\t')
		pos = int(data[1])
		coverage = int(data[3])
		if start <= pos < end:
			x.append(fileid)
			y.append(coverage)
			if not pos in coverageData:
				coverageData[pos] = 0
			if not pos in counts:
				counts[pos] = 0
			coverageData[pos] += coverage
			counts[pos] += 1
plt.figure(figsize=(16,8))
sns.boxplot(x,y,order=sorted(list(set(x))))
plt.ylabel('\\footnotesize{coverage}')
plt.savefig(snakemake.output[0])

x = []
y = []

for p in coverageData:
	x.append(p)
	y.append(coverageData[p]/counts[p])
plt.clf()
plt.figure(figsize=(10,6))

plt.plot(x,y)
plt.xlabel('\\footnotesize{position}')
plt.xticks(rotation=90)
plt.ylabel('\\footnotesize{coverage}')
plt.savefig(snakemake.output[1])

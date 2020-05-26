include: "scripts/shared.py"

from snakemake.utils import validate

#Validate configuration files

configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")


#Generate Input/Output Files from specified folder
inputIDs, = glob_wildcards('data/input/'+config['input_folder']+'/{id}'+config['input_read_1_ending'])
kmer_lengths = config['kmers']


#kmer_lengths = [24]

#Helper function that assembles required input files dependent on configuration settings
def get_input():
	input_list = []

	if config['generative_model']:
		input_list.append(expand('data/output/'+config['input_folder']+'/kmers/{kmer}/predictions.probabilistic_gen.tsv',kmer=kmer_lengths))
	if config['probabilistic_model']:
		input_list.append(expand('data/output/'+config['input_folder']+'/kmers/{kmer}/predictions.probabilistic_cov.tsv',kmer=kmer_lengths))
		if config['plot_top3_fit']:
			input_list.append(expand('data/output/'+config['input_folder']+'/kmers/{kmer}/{id}_top3fit.svg',kmer=kmer_lengths,id=inputIDs))
	if config['distance_model']:
		input_list.append(expand('data/output/'+config['input_folder']+'/kmers/{kmer}/predictions.euclidean.tsv',kmer=kmer_lengths))
	if config['assembly_model']:
		input_list.append(expand('data/output/'+config['input_folder']+'/{id}/exactMatches.tsv',id=inputIDs))
	if config['calc_strand_bias']:
		input_list.append(expand('data/output/'+config['input_folder']+'/{id}/strandbias.txt',id=inputIDs))
	if config['mapping_diff_analysis']:
		input_list.append(expand('data/output/'+config['input_folder']+'/methodAnalysis/{id}/mapping.comparison',id=inputIDs))
	if config['map_filtered_reads']:
		input_list.append(expand('data/output/'+config['input_folder']+'/methodAnalysis/{id}/alignmentToGroundTruthType.sorted.bam.bai',id=inputIDs))
	if config['verifyUniqueness']:
		input_list.append(expand('data/output/kmers/{kmer}/uniquenessTest.tsv',kmer=kmer_lengths))
	if config['kmer_stats_analysis']:
		input_list.append(expand('data/auxiliary/'+config['input_folder']+'/kmers/{kmer}/{id}/stats.tsv',kmer=kmer_lengths,id=inputIDs))

	return input_list

rule all:
	input:
		get_input()


##### load rules #####
include: "rules/assembly.snk"
include: "rules/shared.snk"
include: "rules/kmerApproach.snk"
include: "rules/coverageBased.snk"
include: "rules/probabilistic.snk"
include: "rules/euclidean.snk"

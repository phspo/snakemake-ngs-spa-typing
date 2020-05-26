Snakemake workflow for spa-typing of Next-Generation Sequencing data.

Execution of the workflow requires Snakemake, Conda and Singularity.

# Execution

Input reads (.fastq) should be moved to (or symlinked into) a folder in data/input.

Add a reference genome folder to the data/input folder containing a .fasta file as well as a gene table that contains an entry for protein A.

Add a spa repeat definitions file in .fasta format as well as a spa types definition file to the data/input folder. (Both files can be downloaded from https://spaserver.ridom.de/)

Create a copy of the config.example.yaml

​    mv config.example.yaml config.yaml

and modify the following entries:

spa_repeats_file (change to the name of your spa repeat definition file)

spa_types (change to the name of your spa type definition file)

input_folder (change to the name of your input reads folder)

input_read_1_ending,  input_read_2_ending (change to the file endings that mark your samples)

Execute with

​    snakemake --use-conda --use-singularity --cores 8

to run with 8 cores




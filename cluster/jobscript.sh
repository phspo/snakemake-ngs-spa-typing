#!/bin/bash

#Load Conda and Snakemake and Singularity for cluster execution if a module system exists on the HPC used for computation
#Depending on your HPC system's architecture you might have to adjust this file
module load Singularity/3.5.2
module load Miniconda/3_snakemake
module load Snakemake/5.10.0


#properties = {properties}

{exec_job}

#!/usr/bin/env python
from functions import *

singularity: "docker://vjbaskar/biotools:1.0"

##### Input values
# Genome raw files
GENOME_NAME = config["genome"]
ref_folder = "refs"
createdir(ref_folder)

#####
GENOME_FA = ref_folder + "/" + GENOME_NAME + ".fa"
GENOME_FAI = ref_folder + "/" + GENOME_NAME + ".fa.fai"


"""
rule all:
    input:
        fa = GENOME_FA, fai = GENOME_FAI
"""
# Download the genome if the ref is not provided
# Here the GENOME_NAME is supplied as params.genome into the shell command
rule download_fa:
    params:
        genome = GENOME_NAME
    output: 
        fa = GENOME_FA, fai = GENOME_FAI
    shell:
        """
        wget https://hgdownload.soe.ucsc.edu/goldenPath/{params.genome}/bigZips/{params.genome}.fa.gz -O {output.fa}.gz
        gunzip < {output.fa}.gz > {output.fa}
        samtools faidx {output.fa}
        """

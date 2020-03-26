#!/usr/bin/env python3
# Genome raw files
GENOME_NAME = config["genome"]

GENOME_FA= ref_folder + "/" + GENOME_NAME + ".fa"
GENOME_FAI=ref_folder + "/" + GENOME_NAME + ".fa.fai"

# BWA indexing file
BWA_INDEX_FILE =  ref_folder + "/" + GENOME_NAME + ".fa.sa"


# The main rule where the final output files are requested
rule all:
    input: BAMS
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
       # cp tmp_folder/{params.genome}.fa {output.fa}
       # rm -f {params.genome}.fa.gz
        samtools faidx {output.fa}
        #samtools faidx {output.fa} chr17 > temp
        #mv temp {output.fa}
        """


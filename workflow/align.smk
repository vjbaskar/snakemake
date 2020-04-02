#!/usr/bin/env python

singularity: "docker://vjbaskar/biotools:1.0"
from functions import *

# Input parameters from main snakemake file 
# SAMPLES = ["sample1","sample2" ...]
# GENOME_NAME="hg38" 

#sample_list = SAMPLES
gname = GENOME_NAME

# Create a reference folder
ref_folder = "refs"
createdir(ref_folder)

GENOME_FA= ref_folder + "/" + gname + ".fa"
GENOME_FAI=ref_folder + "/" + gname + ".fa.fai"
BWA_INDEX_FILE =  ref_folder + "/" + gname + ".fa.sa"

print(expand("bams/{sample}.bam", sample = SAMPLES))
BAMS=expand("bams/{sample}.bam", sample = SAMPLES)
# Download the genome if the ref is not provided
# Here the GENOME_NAME is supplied as params.genome into the shell command
rule download_fa:
    params:
        genome = gname
    output: 
        fa = GENOME_FA, fai = GENOME_FAI
    message:
        "Downloading genome files"
    shell:
        """
        wget https://hgdownload.soe.ucsc.edu/goldenPath/{params.genome}/bigZips/{params.genome}.fa.gz -O {output.fa}.gz
        gunzip < {output.fa}.gz > {output.fa}
        samtools faidx {output.fa}
        """

# Generate BWA index
# Note that the output is several files but we detect the last one created.
#       May need to change this later.
rule bwa_index:
    input:
        fa = GENOME_FA
    output:
        index_sa = BWA_INDEX_FILE
    message:
        "Generating bwa index"
    shell:
        """
        bwa index {input.fa}
        """


# Aligning using bwa mem
rule align:
    input:
        fa = GENOME_FA,
        f1 = "fastq/{sample}_R1.fq.gz",
        f2 = "fastq/{sample}_R2.fq.gz",
        index_sa = BWA_INDEX_FILE
    output:
        "bams/{sample}.raw.bam"
    shell:
        """
        bwa mem {input.fa} {input.f1} {input.f2} | samtools view -Sb - -o {output}
        """

rule sortbam:
    input:
        "bams/{sample}.raw.bam"
    output:
        "bams/{sample}.bam"
    shadow:
        "shallow"
    shell:
        """
        mkdir -p metrics
        samtools sort {input} -o {output}.tempfile
        java -jar /picard.jar MarkDuplicates I={output}.tempfile O={output} M=metrics/marked_dup_metrics.txt
        """

#!/usr/bin/env python
import pandas


"""
The needed fields in config file are:
    sample_filename = "something.tsv"
    genome = "hg38" etc
"""
sample_file = config["sample_file"]
data = pandas.read_csv(sample_file)
data.columns = ["sample_id","sample_name","condition"]

ref_folder = "refs"
SAMPLES = data['sample_id'].to_list()

# Get the names of the final output files
BAMS = expand("bams/{sample}.bam", sample = SAMPLES)
FASTQS = expand(["fastq/{sample}_R1.fq.gz", "fastq/{sample}_R2.fq.gz"], sample = SAMPLES)
#GENOME_NAME = "hg38"
GENOME_NAME = config["genome"]
print(f"Your genome = {GENOME_NAME}")

# Genome raw files
GENOME_FA= ref_folder + "/" + GENOME_NAME + ".fa"
GENOME_FAI=ref_folder + "/" + GENOME_NAME + ".fa.fai"

# BWA indexing file
BWA_INDEX_FILE =  ref_folder + "/" + GENOME_NAME + ".fa.sa"


# The main rule where the final output files are requested
rule all:
    input: FASTQS

include: "workflow/cramtofastq.smk"

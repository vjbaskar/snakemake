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
print(data)
SAMPLES = data['sample_id'].to_list()
print(SAMPLES)
# Get the names of the final output files
BAMS = expand("bams/{sample}.bam", sample = SAMPLES)
FASTQS = expand(["fastq/{sample}_R1.fq.gz", "fastq/{sample}_R2.fq.gz"], sample = SAMPLES)
GENOME_NAME = config["genome"] # Eg. hg38
print(f"Your genome = {GENOME_NAME}")


# The main rule where the final output files are requested
rule all:
    input: BAMS

include: "workflow/align.smk"





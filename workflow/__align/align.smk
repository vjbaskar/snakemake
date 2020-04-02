#!/usr/bin/env python

singularity: "docker://vjbaskar/biotools:1.0"


"""
# Get the names of the final output files
BAMS = expand("bams/{sample}.bam", sample = SAMPLES)
GENOME_NAME = "hg38"

# Genome raw files
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

# Generate BWA index
# Note that the output is several files but we detect the last one created.
#       May need to change this later.

rule bwa_index:
    input:
        fa = GENOME_FA
    output:
        index_sa = BWA_INDEX_FILE
    shell:
        """
        bwa index {input.fa}
        """


# Aligning using bwa mem
rule align:
    input:
        fa = GENOME_FA,
        f1 = "fastq/{sample}.pair1.fq.gz",
        f2 = "fastq/{sample}.pair2.fq.gz",
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
        #rm -f {output}.tempfile
        """

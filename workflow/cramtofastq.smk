#!/usr/bin/env python

singularity: "docker://vjbaskar/biotools:1.0"
import os
from functions import *


## directory definitions

ref_folder = "refs"
temp_folder = "temp_bams"

#createdir("cramlist")
createdir("bams")
createdir("fastq")
createdir(ref_folder)
createdir(temp_folder)
## Genome fa data
GENOME_FA= ref_folder + "/" + GENOME_NAME + ".fa"
GENOME_FAI=ref_folder + "/" + GENOME_NAME + ".fa.fai"

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
# For a given sample there may be one or more cram files.
# Convert these crams to bams
# Put all the bam file names of the sample into a file called sample.merge.txt
rule collatecrams:
    input:
        sample_folder = "data/{sample}",
        fa = GENOME_FA

    output:
        "bams/{sample}.merge.txt"
    shell:
        """
        crams=`ls {input.sample_folder}/finalfiles/*.cram`
        for file in ${{crams}}
        do
            f=`echo ${{file}} | sed -e 's/.cram/.bam/g'`
            samtools view -b -T {input.fa} ${{file}} -o - | samtools sort - -o ${{f}}
        done
        ls {input.sample_folder}/finalfiles/*.bam > {output}
        """
# Merge the bam files

rule mergebams:
    input:
        "bams/{sample}.merge.txt"
    output:
        "bams/{sample}.bam"
    shell:
        """
        samtools merge -b {input} {output}
        """

# Convert bam to fastq

rule bam2fq:
    input:
        bamfile = "bams/{sample}.bam"
    output:
        f1 = "fastq/{sample}_R1.fq.gz", f2 = "fastq/{sample}_R2.fq.gz"
    shell:
        """
        java -jar /picard.jar SamToFastq I={input.bamfile} F={output.f1} F2={output.f2}
        """


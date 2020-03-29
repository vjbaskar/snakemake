#!/usr/bin/env python
# Author: vijay <vjbaskar@gmail.com>
"""
Converts multiple cram files per sample to pe fastq.
Input vars: 
    genome: name of the genome used to create the cram file. 
            If you have non ucsc names, for eg: ensembl_GRCh37 etc. then you have to create a folder called "refs" and copy the genome.fa and genome.fai there.
Config file:
    profile/lsf/config.yaml: lsf config file to be used. 
    YOu can override it within the rule.
    bsub -M {resources.memory} -R \"select[mem>{resources.memory}] rusage[mem={resources.memory}] span[hosts={resources.cpus}]\" -o logs/{rule}.{wildcards}.cluster -e logs/{rule}.{wildcards}.cluster -J {rule}.{wildcards} -q normal -n {resources.cores}

"""

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

    resources:
        memory=lambda wildcards, attempt: attempt * 1000
    shell:
        """
        crams=`ls {input.sample_folder}/finalfiles/*.cram`
        rm -f {output}
        for file in ${{crams}}
        do
            f=`echo ${{file}} | sed -e 's/.cram/.bam/g'`
            samtools view -b -T {input.fa} ${{file}} -o - | samtools sort - -o ${{f}}
            echo "${{f}}" >> {output}
        done
        #ls {input.sample_folder}/finalfiles/*.bam > {output}
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


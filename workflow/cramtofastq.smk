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


rule cram2bam:
    input:
        "bams/{sample}.merge.txt"
    output:
        "bams/{sample}.bam"
    shell:
        """
        samtools merge -b {input} {output}
        """

rule bam2fq:
    input:
        bamfile = "bams/{sample}.bam"
    output:
        f1 = "fastq/{sample}_R1.fq.gz", f2 = "fastq/{sample}_R2.fq.gz", s = "fastq/{sample}_singletons.fq.gz", f = "fastq/{sample}.fq.gz"
    shell:
        """
        echo {input.bamfile} > {output.f1}
        echo {input.bamfile} > {output.f2}
        samtools fastq -1 {output.f1} -2 {output.f2} -o {output.f} -s {output.s} {input.bamfile} 
        """


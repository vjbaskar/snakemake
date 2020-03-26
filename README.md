# Snakemake file for align and marking duplicates

## Depends on:
   
   * Singularity
   * docker://vjbaskar/biotools:1.0 - has samtools, bwa and picard.
   ..* This docker can be generated from the github release https://github.com/vjbaskar/biotools/releases/tag/biotools-1.0


## Config files:
   
   * `config/smk.yaml`: main config file for snakemake
   * `config/lsf.json`: cluster config file for LSF

## Input:
    
    Set of sample names. 
    fastq/{sample_name}.pair1.fq.gz, fastq/{sample_name}.pair2.fq.gz

## Indexing.
    
    config/smk.yaml contains the genome name that you should change.
    This `config['genome']` is called in `Snakefile`
    The genome will be downloaded to refs from ucsc genome browser and indexed using bwa

## Align.

    

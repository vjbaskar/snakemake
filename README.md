# Snakemake file for align and marking duplicates

## Depends on:
   
   * Singularity
   * docker://vjbaskar/biotools:1.0 - has samtools, bwa and picard.
        - This docker can be generated from the github release https://github.com/vjbaskar/biotools/releases/tag/biotools-1.0


## Config files:
   
   * `config/smk.yaml`: main config file for snakemake
   * `config/lsf.json`: cluster config file for LSF

## Scripts:
   * The main script is `Snakemake`
   * The individual workflows are written in `workflow/*`
   * Python functions to be imported are in `scripts/*`
   * Example run is given in `snakemake.sh`

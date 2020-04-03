#!/bin/bash

# If you want use older clusterconfig use the following options
# --cluster " bsub -M {cluster.memory} -R {cluster.resources} -o {cluster.output} -e {cluster.error} -J {cluster.name} -q {cluster.queue} -n {cluster.nCPUs} " \
# --cluster-config config/lsf.json \



fullrun(){
    snakemake --use-singularity \
    --jobs 999 \
    --printshellcmds \
    --rerun-incomplete \
    --configfile config/smk.yaml \
    --profile profile/lsf/ \
    --restart-times 5 \
    --jobname "snakejob.{name}.{jobid}" \
    -p $@ #-n #--report report.html  #-n #-s align.snk

}

localrun(){
    snakemake --use-singularity \
    --printshellcmds \
    --rerun-incomplete \
    --configfile config/smk.yaml \
    -p $@ #-n #--report report.html  #-n #-s align.snk
}

fullrun -s bwa_align.smk

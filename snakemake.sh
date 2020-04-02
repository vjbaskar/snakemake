#mkdir -p bams logs/cluster metrics refs cramlist

fullrun(){
snakemake --use-singularity \
    --jobs 999 \
    --printshellcmds \
    --rerun-incomplete \
    --configfile config/smk.yaml \
    --cluster " bsub -M {cluster.memory} -R {cluster.resources} -o {cluster.output} -e {cluster.error} -J {cluster.name} -q {cluster.queue} -n {cluster.nCPUs} " \
    -s bwa_align.smk \
    -p $@ #-n #--report report.html  #-n #-s align.snk
}

localrun(){
    snakemake --use-singularity \
    --printshellcmds \
    --rerun-incomplete \
    --configfile config/smk.yaml \
    -s bwa_align.smk \
    -p $@ #-n #--report report.html  #-n #-s align.snk
}

fullrun

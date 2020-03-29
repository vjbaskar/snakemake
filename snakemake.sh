mkdir -p bams logs/cluster metrics refs cramlist
snakemake --use-singularity \
    --jobs 999 \
    --printshellcmds \
    --rerun-incomplete \
    --configfile config/smk.yaml \
    --profile profile/lsf/ \
    --restart-times 5 \
    --show-failed-logs \
    -p $@ #-n #--report report.html  #-n #-s align.snk
    #--cluster " bsub -M {cluster.memory} -R {cluster.resources} -o {cluster.output} -e {cluster.error} -J {cluster.name} -q {cluster.queue} -n {cluster.nCPUs} " \
    #--cluster-config config/lsf.json \
    #--configfile config/smk.yaml \


rule subsample_read:
    input:
        lambda wildcards:  os.path.join(raw_fastq_dir[wildcards.sample], "{cell_line}.fastq.gz")
    output:
        os.path.join(subsample_dir, 
                        f"subsample_n_{subsample_read_n}",
                        "{sample}/{cell_line}.fastq")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        n_reads = subsample_read_n,
        seed = config['random_seed']
    shell:
        """
        mkdir -p $(dirname {output})
        seqtk sample -s {params.seed} {input} {params.n_reads} > {output}
        """
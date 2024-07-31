configfile: "config/input_config.yaml"
results_dir = config["output_path"]

rule run_cellranger:
    input:
        #expand(".flag/{sample_id}_cellranger.done", sample_id=config["sample_name"])
        expand(".flag/{sample_id}_cellranger.done", sample_id=config["sample_name"])
    output:
        touch(".flag/cellranger.done")


rule _cellranger_mkref:
    input:
        gtf=config['reference']['gtf'],
        genome=config['reference']['genome']
    output:
        os.path.join(results_dir, 'cellranger/refdata_genecodev44')
    resources:
        mem_mb=256000,
        cpus_per_task=32
    params:
        genome='GRCh38'
    shell:
        """
        module load cellranger/8.0.1
        cellranger mkref --genome={params.genome} \
            --fasta={input.genome} \
            --genes={input.gtf} \
            --memgb=256 \
            --nthreads={resources.cpus_per_task} \
            --output-dir={output}
        """

rule _cellranger_count_sc:
    input:
        fastq = "/home/users/allstaff/you.yu/LongBench/sequencing_data/illumina_sc",
        ref = os.path.join(results_dir, 'cellranger/refdata_genecodev44')
    output:
        flag = touch('.flag/ill_sc_cellranger.done'),
        dir = report(directory(os.path.join(results_dir, 'cellranger/ill_sc')), htmlindex="outs/web_summary.html")
    resources:
        mem_mb=200000,
        cpus_per_task=16
    params:
        sample_id = "M000495_batch1_2_GEX"
    shell:
        """
        module load cellranger/8.0.1
        cellranger count --id=sc \
           --transcriptome={input.ref} \
           --fastqs={input.fastq} \
           --sample={params.sample_id} \
           --create-bam=true \
           --localcores={resources.cpus_per_task} \
           --output-dir={output.dir}
        """

rule _cellranger_count_sn:
    input:
        fastq = "/home/users/allstaff/you.yu/LongBench/sequencing_data/illumina_sn",
        ref = os.path.join(results_dir, 'cellranger/refdata_genecodev44')
    output:
        flag = touch('.flag/ill_sn_cellranger.done'),
        dir = report(directory(os.path.join(results_dir, 'cellranger/ill_sn')), htmlindex="outs/web_summary.html")
    resources:
        mem_mb=200000,
        cpus_per_task=16
    params:
        sample_id = "M000495_batch1_1_GEX"
    shell:
        """
        module load cellranger/8.0.1
        cellranger count --id=sn \
           --transcriptome={input.ref} \
           --fastqs={input.fastq} \
           --sample={params.sample_id} \
           --create-bam=true \
           --localcores={resources.cpus_per_task} \
           --output-dir={output.dir}
        """
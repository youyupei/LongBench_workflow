# Required variables (set in sr_sc_sn Snakefile):
## config["samples_fastq_dir"]  — dict: sample_name -> fastq directory
## config["sample_name"]        — list of sample names
## results_dir

######## Subsample ########
rule subsample_1M_sr_sc_sn:
    input:
        lambda wildcards: config["samples_fastq_dir"][wildcards.sample_name]
    output:
        temp(os.path.join(config["scratch_dir"], "sr_sc_sn/qc/subsample_fq/{sample_name}_subsampled1M.fastq"))
    resources:
        cpus_per_task=1,
        mem_mb=8000
    params:
        n_reads = 1000000,
        seed = config['random_seed']
    shell:
        """
        mkdir -p $(dirname {output})
        cat {input}/*_R2_001.fastq.gz | seqtk sample -s {params.seed} - {params.n_reads} > {output}
        """

######## NanoPlot ########
rule NanoPlot_sr_sc_sn:
    input:
        reads=rules.subsample_1M_sr_sc_sn.output
    output:
        os.path.join(results_dir, "qc/NanoPlot/{sample_name}/NanoPlot-data.tsv.gz")
    conda:
        config['conda']['NanoPlot']
    resources:
        cpus_per_task=16,
        mem_mb=16000
    shell:
        """
        output_dir=$(dirname {output})
        mkdir -p $output_dir
        NanoPlot --fastq {input.reads} --outdir $output_dir -t {resources.cpus_per_task} --raw --tsv_stats
        """

######## Subsample BAM for AlignQC ########
rule subsample_bam_for_alignqc:
    input:
        bam = os.path.join(results_dir, "cellranger/{sample_name}/outs/possorted_genome_bam.bam")
    output:
        bam = temp(os.path.join(config["scratch_dir"], "sr_sc_sn/qc/subsample_bam/{sample_name}_subsample10M.bam")),
        bai = temp(os.path.join(config["scratch_dir"], "sr_sc_sn/qc/subsample_bam/{sample_name}_subsample10M.bam.bai"))
    params:
        seed    = config['random_seed'],
        n_reads = "10M"
    resources:
        cpus_per_task = 1,
        mem_mb        = 16000
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view -F 4 {input.bam} | cut -f1 | sort -u > {output.bam}.read_ids
        shuf -n $(numfmt --from=si {params.n_reads}) \
            --random-source=<(yes {params.seed}) {output.bam}.read_ids \
            > {output.bam}.read_ids.subsampled
        rm {output.bam}.read_ids
        samtools view -b -N {output.bam}.read_ids.subsampled {input.bam} > {output.bam}
        samtools index {output.bam}
        rm {output.bam}.read_ids.subsampled
        """

######## AlignQC ########
rule alignQC_analysis_single:
    input:
        fa      = config['reference']['genome'],
        anno    = config['reference']['gtf_gz'],
        bam     = rules.subsample_bam_for_alignqc.output.bam,
        bai     = rules.subsample_bam_for_alignqc.output.bai
    output:
        out_dir = directory(os.path.join(results_dir, "qc/AlignQC/{sample_name}/")),
        tmp_dir = temp(directory(os.path.join(config["scratch_dir"], "sr_sc_sn/alignQC_tmp/{sample_name}")))
    resources:
        cpus_per_task = 16,
        mem_mb        = 200000,
        slurm_extra = "--partition=long --exclude=med-n02"
    priority: 10
    container: "docker://vacation/alignqc"
    shell:
        """
        mkdir -p {output.tmp_dir}
        alignqc analyze {input.bam} \
            -g {input.fa} \
            --gtf {input.anno} \
            --output_folder {output.out_dir} \
            --threads {resources.cpus_per_task} \
            --specific_tempdir {output.tmp_dir}
        """

rule alignQC_analysis:
    input:
        expand(
            rules.alignQC_analysis_single.output,
            sample_name=config["sample_name"]
        )
    output:
        touch(os.path.join(results_dir, ".flag/alignQC.done"))
    

######## Count bases ########
rule count_bases_in_fastq_sr_sc_sn:
    """Count total bases sequenced from raw FASTQ directory (all .fastq.gz files)."""
    input:
        lambda wildcards: config["samples_fastq_dir"][wildcards.sample_name]
    output:
        os.path.join(results_dir, "qc/base_counts/{sample_name}.total_bases")
    resources:
        cpus_per_task=1,
        mem_mb=4000
    shell:
        """
        mkdir -p $(dirname {output})
        zcat {input}/*.fastq.gz | awk 'NR%4==2{{sum+=length($0)}} END{{print sum}}' > {output}
        """


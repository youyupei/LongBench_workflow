configfile: "config/config.yaml"
results_dir = config["output_path"]

###################### Run QC pipelines ###############################
rule all_qc:
    input:
        expand(
            [
                os.path.join(results_dir, "qc/{sample}.{read_length_source}.read_length.txt"),
                os.path.join(results_dir, "qc/{sample}_coverage_picard.RNA_Metrics"),
                os.path.join(results_dir, "qc/{sample}_QualityScoreDistribution_picard.pdf"),
                os.path.join(results_dir, "qc/{sample}_align2genome.bam.picard_AlignmentSummaryMetrics.txt"),
                os.path.join(results_dir, "qc/{sample}_flame_coverage_plot.{flames_cov_plot_suffix}"),
                ".flag/{sample}_sqanti3.done"
            ],
            sample=config["sample_id"],
            read_length_source=["raw", "blaze_demux"],
            flames_cov_plot_suffix = ['png', 'CDS.png']
        )
    output:
        os.path.join(results_dir, "qc/.flag/all_qc.done")
    shell:
        """
        mkdir -p $(dirname {output})
        touch {output}
        """
###################### Rules ###############################

## Read length related
rule generate_demux_read_length_statistics:
    input:
        os.path.join(results_dir, "flames_out/{sample}/matched_reads.fastq")
    output:
        os.path.join(results_dir, "qc/{sample}.blaze_demux.read_length.txt")
    resources:
        cpus_per_task=1,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        mkdir -p $(dirname {output})

        # remove output file if it already exists
        if [ -f {output} ]; then
            rm {output}
        fi

        # if input is a directory, iterate through all the fastq 
        if [ -d {input} ]; then
            for f in {input}/*.fastq; do
                awk '{{if(NR%4==2){{print length($0)}}}}' $f >> {output}
            done
        else
            awk '{{if(NR%4==2){{print length($0)}}}}' {input} > {output}
        fi
        """

## Read length related
rule generate_raw_read_length_statistics:
    input:
        lambda wildcards: config["samples_fastq_dir"][wildcards.sample]
    output:
        os.path.join(results_dir, "qc/{sample}.raw.read_length.txt")
    resources:
        cpus_per_task=1,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        mkdir -p $(dirname {output})

        # remove output file if it already exists
        if [ -f {output} ]; then
            rm {output}
        fi

        # if input is a directory, iterate through all the fastq 
        if [ -d {input} ]; then
            for f in {input}/*.fastq; do
                awk '{{if(NR%4==2){{print length($0)}}}}' $f >> {output}
            done
        else
            awk '{{if(NR%4==2){{print length($0)}}}}' {input} > {output}
        fi
        """

# QualityScoreDistribution
rule picard_QualityScoreDistribution:
    input:
        bam = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam"),
        picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar"
    output:
        chart = os.path.join(results_dir, "qc/{sample}_QualityScoreDistribution_picard.pdf"),
        txt = os.path.join(results_dir, "qc/{sample}_QualityScoreDistribution_picard.txt")
    resources:
        cpus_per_task=1,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        module load picard-tools/2.26.11
        mkdir -p $(dirname {output.chart})
        java -jar {input.picard_dir} QualityScoreDistribution \
            I={input.bam} \
            O={output.txt} \
            CHART={output.chart}
        """

# Coverage
rule picard_coverage_data:
    input:
        bam = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam"),
        refFlat = '/home/users/allstaff/you.yu/LongBench/reference_files/GRCh38/refFlat.txt',
        picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar",
        flag = ".flag/flames_{sample}.done"
    output:
        os.path.join(results_dir, "qc/{sample}_coverage_picard.RNA_Metrics")
    resources:
        cpus_per_task=16,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        module load picard-tools/2.26.11
        mkdir -p $(dirname {output})
        java -jar {input.picard_dir} \
            CollectRnaSeqMetrics  \
            I={input.bam} O={output} \
            REF_FLAT={input.refFlat} \
            STRAND_SPECIFICITY=NONE
        """

# BamIndexStats
rule picard_CollectAlignmentSummaryMetrics:
    input:
        bam = os.path.join(results_dir,"flames_out/{sample}/{bam}"),
        picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar",
        ref = config['reference']['genome']
    output:
        os.path.join(results_dir, "qc/{sample}_{bam}.picard_AlignmentSummaryMetrics.txt")
    resources:
        cpus_per_task=16,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        module load picard-tools/2.26.11
        mkdir -p $(dirname {output})
        java -jar {input.picard_dir} CollectAlignmentSummaryMetrics \
          R={input.ref} \
          I={input.bam} \
          O={output}
        """

rule flame_coverage_plot:
    input:
        bam = os.path.join(results_dir,"flames_out/{sample}/realign2transcript.bam"),
        gtf = config['reference']['gtf']
    output:
        os.path.join(results_dir, "qc/{sample}_flame_coverage_plot.{suffix}")
    resources:
        cpus_per_task=1,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        mkdir -p $(dirname {output})
        Rscript scripts/plot_coverage.R {input.bam} {input.gtf} {output}
        """


## SQANTI3

rule sqanti3:
    input:
        reads=os.path.join(results_dir, "flames_out/{sample}/matched_reads_subsampled.fastq"),
        gtf=config['reference']['gtf'],
        genome=config['reference']['genome'],
    
    output:
        touch(".flag/{sample}_sqanti3.done")
    conda: 
        config['software']['sqanti3_env']
    resources:
        cpus_per_task=16,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params: 
        sqanti3_qc_script = os.path.join(config['software']['sqanti3_dir'], "sqanti3_qc.py"),
        outidr = lambda w:  os.path.join(results_dir, f"qc/sqanti3/{w.sample}"),
        additional_arg = "--skipORF"
    shell:
        """
        {params.sqanti3_qc_script} {input.reads} {input.gtf} {input.genome} {params.additional_arg} --fasta  --cpus {resources.cpus_per_task} --force_id_ignore -d {params.outidr} --report pdf
        """
    
rule _subsample_reads_for_sqanti3:
    input:
        os.path.join(results_dir, "flames_out/{sample}/matched_reads.fastq")
    output:
        os.path.join(results_dir, "flames_out/{sample}/matched_reads_subsampled.fastq")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        n_reads = 200
    shell:
        """
        seqtk sample -s100 {input} {params.n_reads} > {output}
        """


# rule get_intronic_pct:
#     input:
#         expand(os.path.join(results_dir, "qc/{sample}_coverage_picard.RNA_Metrics"), sample=config['sample_id'])
#     output:
#         os.path.join(results_dir, "qc/coverage_intronic_pct.txt")
#     localrule: True
#     run:
#         import pandas as pd
#         import os

#         def read_metrics_class(file):
#             with open(file, 'r') as f:
#                 lines = f.readlines()
#                 for i, line in enumerate(lines):
#                     if line.startswith("## METRICS CLASS"):
#                         matrix_df = pd.read_csv(file, skiprows=i+1, nrows=1, sep="\t")
#                     if line.startswith("## HISTOGRAM"):
#                         coverage_df = pd.read_csv(file, skiprows=i+1, sep="\t")
#                         break
#             return matrix_df, coverage_df 
        
#         intronic_pcts = []
#         # Create the directory if it does not exist
#         if not os.path.exists(os.path.dirname(output)):
#             os.makedirs(os.path.dirname(output))

#         with open(output, 'w') as f:
#             for i, file in enumerate(input):
#                 df, _ = read_metrics_class(file)
#                 intronic_pct = df['PCT_INTRONIC_BASES'].values[0]
#                 f.write(f"{config['sample_id'][i]}\t{intronic_pct}\n")
#                 intronic_pcts.append(intronic_pct)


# rule seqkit_qc_state:
#     input:
#         config["samples_fastq_dir"]["ont_sn"],
#         config["samples_fastq_dir"]["ont_sc"],
#         config["output_path"] + "/flames_out/ont_sn/matched_reads.fastq",
#         config["output_path"] + "/flames_out/ont_sc/matched_reads.fastq"
#     output:
#         os.path.join(results_dir, "qc/seqkit_qc_state.tsv")
#     resources:
#         cpus_per_task=8,
#         mem_mb=32000,
#         slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
#     shell:
#         """
#         mkdir -p $(dirname {output})
#         /home/users/allstaff/you.yu/software/bin/seqkit stats  -Ta {input} > {output}
#         """
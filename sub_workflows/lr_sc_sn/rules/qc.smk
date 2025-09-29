# Required variables
## config["sample_id"]
## config['conda']['sqanti3']
## config["samples_fastq_dir"]
## config["reference"]["picard_reference"]
## config["reference"]["genome"]
## config["reference"]["gtf"]
## config["random_seed"]
## config["software"]["sqanti3_dir"]
## results_dir

######## Demultiplexing ########
# NOTE: this is the only QC rule that uses all reads, the following rules use reads successfully assigned to a cell

rule Demultiplexing_plot:
    input:
        expand(
            os.path.join(results_dir, "flames_out/{sample}/summary.txt"),
            sample=config["sample_id"]
        )
    output:
        report(os.path.join(results_dir, "plots/qc/demultiplexing_state.pdf"), category="QC", subcategory="SC_demultiplex")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        sample_names = config["sample_id"]
    localrule: True
    script:
        os.path.join(config['sub_wf_dir'],'scripts/plot_demultiplexing_state.R')



######## Subsample ########
rule subsample_2M_reads:
    input:
        lambda wildcards: config["samples_fastq_dir"][wildcards.sample]
    output:
        os.path.join(scratch_dir, "qc/subsample_fq/{sample}_matched_reads_subsampled.fastq")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        n_reads = 2000000,
        seed = config['random_seed']
    shell:
        """
        mkdir -p $(dirname {output})
        seqtk sample -s {params.n_reads} {input} {params.n_reads} > {output}
        """


rule Filter_and_Subsample_bam_for_qc:
    input:
        bam = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam")
    output:
        bam = temp(os.path.join(scratch_dir,"subsample_data/{sample}/genome_map_subsample_rate_{subsample_rate}.bam")),
        bai = temp(os.path.join(scratch_dir,"subsample_data/{sample}/genome_map_subsample_rate_{subsample_rate}.bam.bai"))
    resources:
        cpus_per_task=1,
        mem_mb=8000
    params:
        seed = config['random_seed'],
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view --subsample {wildcards.subsample_rate} --subsample-seed {params.seed} -F2034 {input.bam} -h | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

######## NanoPlot ########
rule NanoPlot:
    input:
        reads=os.path.join(scratch_dir, "qc/subsample_fq/{sample}_matched_reads_subsampled.fastq")
    output:
        os.path.join(results_dir, "qc/NanoPlot/{sample}/NanoPlot-data.tsv.gz")
    conda:
        config['conda']['NanoPlot']
    resources:
        cpus_per_task=16,
        mem_mb=16000
    shell:
        """
        output_dir=$(dirname {output})
        mkdir -p $output_dir
        NanoPlot --fastq {input.reads} --outdir $output_dir -t {resources.cpus_per_task} --raw  --tsv_stats
        """

# Coverage
rule picard_coverage_data:
    input:
        bam = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam"),
        refFlat = config['reference']['picard_reference'],
        picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar",
        flag = os.path.join(results_dir, ".flag/flames_{sample}.done")
    output:
        os.path.join(results_dir, "qc/coverage/{sample}.picard.RNA_Metrics")
    resources:
        cpus_per_task=16,
        mem_mb=16000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        module load picard-tools
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
        mem_mb=16000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        module load picard-tools
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
        report(
            os.path.join(results_dir, "qc/coverage/{sample}.flame.coverage_plot.{suffix}"), 
            category="QC", subcategory="Coverage",  labels=lambda wildcards: {
              "Sample": "{sample}",
              "Region": "CDS" if "CDS" in wildcards.suffix else "Whole transcript",
              "figure": "flame_coverage_plot"
          })
    resources:
        cpus_per_task=1,
        mem_mb=32000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        script = os.path.join(config['main_wf_dir'],'scripts/flames_coverage_plot.R')
    shell:
        """
        mkdir -p $(dirname {output})
        Rscript  {params.script} {input.bam} {input.gtf} {output}
        """

## SQANTI3
rule sqanti3:
    input:
        reads=os.path.join(scratch_dir, "qc/subsample_fq/{sample}_matched_reads_subsampled.fastq"),
        gtf=config['reference']['gtf'],
        genome=config['reference']['genome']
    output:
        dir = directory(os.path.join(results_dir, "qc/sqanti3/{sample}/"))
        # html = report(os.path.join(results_dir, "qc/sqanti3/{sample}/matched_reads_subsampled_SQANTI3_report.html"),
        #             category="QC", subcategory="SQANTI3") # this is not included as it is too large
    conda: 
        config['conda']['sqanti3']
    resources:
        cpus_per_task=16,
        mem_mb=32000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params: 
        sqanti3_qc_script = os.path.join(config['software']['sqanti3_dir'], "sqanti3_qc.py"),
        additional_arg = "--skipORF"
    shell:
        """
        {params.sqanti3_qc_script} {input.reads} {input.gtf} {input.genome} {params.additional_arg} --fasta  --cpus {resources.cpus_per_task} --force_id_ignore -d {output} --report html
        """



rule RSeQC_gene_body_coverage:
    input:
        bed=config['reference']['bed_housekeeping_genes'], 
        genome_bam = os.path.join(scratch_dir,"subsample_data/{sample}/genome_map_subsample_rate_0.1.bam"),
        genome_bai = os.path.join(scratch_dir,"subsample_data/{sample}/genome_map_subsample_rate_0.1.bam.bai")
    output:
        report(os.path.join(results_dir, "qc/RSeQC/{sample}.geneBodyCoverage.curves.pdf"))
    resources:
        cpus_per_task=1,
        mem_mb=16000
    conda:
        config['conda']['RSeQC']
    shell:
        """
        mkdir -p $(dirname {output})
        geneBody_coverage.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}
        """


rule RSeQC_junction_saturation:
    input:
        bed=config['reference']['gtf'] + '.bed', # I have rules to convert gtf to bed
        genome_bam = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam")
    output:
        report(os.path.join(results_dir, "qc/RSeQC/{sample}.junctionSaturation_plot.pdf"))
    resources:
        cpus_per_task=8,
        mem_mb=128000
    conda:
        config['conda']['RSeQC']
    shell:
        """
        mkdir -p $(dirname {output})
        junction_saturation.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}
        """


rule RSeQC_junction_annotation:
    input:
        bed=config['reference']['gtf'] + '.bed', # I have rules to convert gtf to bed
        genome_bam = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam")
    output:
        report(os.path.join(results_dir, "qc/RSeQC/{sample}.splice_events.pdf")),
        report(os.path.join(results_dir, "qc/RSeQC/{sample}.splice_junction.pdf"))
    resources:
        cpus_per_task=1,
        mem_mb=8000
    conda:
        config['conda']['RSeQC']
    shell:
        """
        mkdir -p $(dirname {output[0]})
        junction_annotation.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}
        """



## AlignQC
rule alignQC_analysis:
    input: 
        #fa=config['reference']['genome'], 
        fa="/vast/scratch/users/you.yu/LongBench/reference/GRCh38.primary_assembly.genome.fa",
        anno=config['reference']['gtf_gz'],
        genome_bam = os.path.join(scratch_dir,"subsample_data/{sample}/genome_map_subsample_rate_0.01.bam"),
        bai = os.path.join(scratch_dir,"subsample_data/{sample}/genome_map_subsample_rate_0.01.bam.bai")
    # conda:
    #     config['conda']['AlignQC']
    resources:
        cpus_per_task=16,
        mem_mb=800000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au" # tmp use bonus qos
   
    output:
        output = directory(os.path.join(results_dir, "qc/AlignQC/{sample}/")),
        tmp_dir = temp(directory(os.path.join(scratch_dir, "alignQC_tmp","{sample}")))
    priority: 10
    container: "docker://vacation/alignqc"
    retries: 3
    shell:
        """
        mkdir -p $(dirname {output.output})
        mkdir -p {output.tmp_dir}
        alignqc analyze {input.genome_bam} \
            -g {input.fa} \
            --gtf {input.anno} \
            --output_folder {output.output} \
            --threads {resources.cpus_per_task} \
            --specific_tempdir {output.tmp_dir}
        """

# Pseudobulk QC
rule _pseudobulk_read_count_single_run:
    input:
        fastq = results_dir + "/PseudoBulkAlignment/{sample}_{cell_line}_pseudo_bulk.fastq"
    output:
        txt = temp(results_dir + "/PseudoBulkQC/{sample}_{cell_line}_pseudo_bulk_read_count.txt")
    resources:
        cpus_per_task=1
    shell:
        """
        mkdir -p $(dirname {output.txt})
        expr $(( $(wc -l < {input.fastq}) / 4 )) > {output.txt}
        """

rule pseudobulk_read_count:
    input:
        expand(rules._pseudobulk_read_count_single_run.output[0], 
            cell_line = config['cell_line_list'],
            sample = config['sample_id']) 
    output:
        results_dir + "/PseudoBulkQC/pseudo_bulk_read_count.csv"
    resources:
        cpus_per_task=1
    script:
        "../scripts/combine_read_count.R"




# # Commented out as sqanti3 provides similar information
#     # Read length related
#     rule generate_demux_read_length_statistics:
#         input:
#             os.path.join(results_dir, "flames_out/{sample}/matched_reads.fastq")
#         output:
#             os.path.join(results_dir, "qc/{sample}.blaze_demux.read_length.txt")
#         resources:
#             cpus_per_task=1,
#             mem_mb=32000,
#             slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
#         shell:
#             """
#             mkdir -p $(dirname {output})
#     
#             # remove output file if it already exists
#             if [ -f {output} ]; then
#                 rm {output}
#             fi
#     
#             # if input is a directory, iterate through all the fastq 
#             if [ -d {input} ]; then
#                 for f in {input}/*.fastq; do
#                     awk '{{if(NR%4==2){{print length($0)}}}}' $f >> {output}
#                 done
#             else
#                 awk '{{if(NR%4==2){{print length($0)}}}}' {input} > {output}
#             fi
#             """

    ## Read length related
    # rule generate_raw_read_length_statistics:
    #     input:
    #         lambda wildcards: config["samples_fastq_dir"][wildcards.sample]
    #     output:
    #         os.path.join(results_dir, "qc/{sample}.raw.read_length.txt")
    #     resources:
    #         cpus_per_task=1,
    #         mem_mb=32000,
    #         slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    #     shell:
    #         """
    #         mkdir -p $(dirname {output})
    # 
    #         # remove output file if it already exists
    #         if [ -f {output} ]; then
    #             rm {output}
    #         fi
    # 
    #         # if input is a directory, iterate through all the fastq 
    #         if [ -d {input} ]; then
    #             for f in {input}/*.fastq; do
    #                 awk '{{if(NR%4==2){{print length($0)}}}}' $f >> {output}
    #             done
    #         else
    #             awk '{{if(NR%4==2){{print length($0)}}}}' {input} > {output}
    #         fi
    #         """

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
    #         slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    #     shell:
    #         """
    #         mkdir -p $(dirname {output})
    #         /home/users/allstaff/you.yu/software/bin/seqkit stats  -Ta {input} > {output}
    #         """


# ######## QualityScoreDistribution ######## use nanoPlot instead
# rule picard_QualityScoreDistribution:
#     input:
#         bam = os.path.join(results_dir,"flames_out/{sample}/align2genome.bam"),
#         picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar"
#     output:
#         chart = report(
#             os.path.join(results_dir, "qc/{sample}_QualityScoreDistribution_picard.pdf"), 
#             category="QC", subcategory="QualityScoreDistribution",  
#             labels={
#               "Sample": "{sample}",
#               "figure": "QualityScoreDistribution_picard"
#             }),
#         txt = os.path.join(results_dir, "qc/{sample}_QualityScoreDistribution_picard.txt")
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000
#         #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
#     shell:
#         """
#         module load picard-tools
#         mkdir -p $(dirname {output.chart})
#         java -jar {input.picard_dir} QualityScoreDistribution \
#             I={input.bam} \
#             O={output.txt} \
#             CHART={output.chart}
#         """

###################### Run QC pipelines ###############################
rule qc:
    input:
        expand(
            [
                #os.path.join(results_dir, "qc/coverage/{sample}.flame.coverage_plot.{flames_cov_plot_suffix}"),
                #os.path.join(results_dir, "qc/sqanti3/{sample}"),
                os.path.join(results_dir, "qc/NanoPlot/{sample}/NanoPlot-data.tsv.gz"),
                #os.path.join(results_dir, "qc/RSeQC/{sample}.geneBodyCoverage.curves.pdf"),
                #os.path.join(results_dir, "qc/coverage/{sample}.picard.RNA_Metrics"),
                os.path.join(results_dir, "qc/RSeQC/{sample}.junctionSaturation_plot.pdf"),
                #os.path.join(results_dir, "qc/RSeQC/{sample}.splice_events.pdf"),
                #os.path.join(results_dir, "qc/RSeQC/{sample}.splice_junction.pdf"),
                # rules.alignQC_analysis.output[0],
            ],
            sample=config["sample_id"],
            #read_length_source=["raw"],
            flames_cov_plot_suffix = ['png', 'CDS.png']
        ),
        
        Demultiplexing_plot = os.path.join(results_dir, "plots/qc/demultiplexing_state.pdf")
    output:
        os.path.join(results_dir, "qc/.flag/qc.done")
    shell:
        """
        mkdir -p $(dirname {output})
        touch {output}
        """
###################################################################################

rule all_alignQC_analysis:
    input:
        expand(
            rules.alignQC_analysis.output[0], 
            sample=config["sample_id"])
    output:
        touch(os.path.join(config['flag_dir'], "lr_sc_sn_bulk_alginQC.done"))

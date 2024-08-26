results_dir = config["output_path"]
barcode_list = config['barcode_list']
cell_line_to_barcode = {cl: bc for d in barcode_list for bc, cl in d.items()}



###################### Run QC pipelines ###############################
rule qc:
    input:
        expand(
            [
                os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.flame.coverage_plot.{flames_cov_plot_suffix}"),
                os.path.join(results_dir, "qc/sqanti3/{sample}_{cell_line}/"),
                os.path.join(results_dir, "qc/NanoPlot/{sample}_{cell_line}/"),
                #os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.RSeQC.geneBodyCoverage.pdf"),
                os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.picard.RNA_Metrics")
            ],
            sample=config["sample_id"],
            flames_cov_plot_suffix = ['png'],
            cell_line = cell_line_to_barcode.keys()
        )
    output:
        os.path.join(results_dir, "qc/.flag/qc.done")
    shell:
        """
        mkdir -p $(dirname {output})
        touch {output}
        """
###################################################################################

######## Subsample ########
rule _subsample_1M_reads:
    input:
        lambda wildcards:  os.path.join(config["samples_fastq_dir"][wildcards.sample], "{cell_line}.fastq")
    output:
        os.path.join(results_dir, "qc/subsample_fq/{sample}_{cell_line}_subsampled1M.fastq")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        n_reads = 1000000,
        seed = config['random_seed']
    shell:
        """
        mkdir -p $(dirname {output})
        seqtk sample -s {params.n_reads} {input} {params.n_reads} > {output}
        """


######## Read length (subsampled fastq) ########
# rule generate_read_length_statistics:
#     input:
#         os.path.join(results_dir, "qc/subsample_fq/{sample}.{cell_line}_matched_reads_subsampled.fastq")
#     output:
#         os.path.join(results_dir, "qc/{sample}.{cell_line}.read_length.txt")
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000
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

######## NanoPlot ########
rule NanoPlot:
    input:
        reads=os.path.join(results_dir, "qc/subsample_fq/{sample}_{cell_line}_subsampled1M.fastq")
    output:
        directory(os.path.join(results_dir, "qc/NanoPlot/{sample}_{cell_line}/"))
    conda:
        config['conda']['NanoPlot']
    resources:
        cpus_per_task=16,
        mem_mb=32000
    shell:
        """
        output_dir={output}
        mkdir -p $output_dir
        NanoPlot --fastq {input.reads} --outdir $output_dir -t {resources.cpus_per_task} --raw  --tsv_stats
        """

# Coverage
rule picard_coverage_data:
    input:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
        refFlat = '/home/users/allstaff/you.yu/LongBench/reference_files/GRCh38/refFlat.txt',
        picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar"
    output:
        os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.picard.RNA_Metrics")
    resources:
        cpus_per_task=16,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
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
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
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
        bam = os.path.join(results_dir,"TranscriptAlignment/{sample}_{cell_line}.sorted.bam"),
        gtf = config['reference']['gtf']
    output:
        report(
            os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.flame.coverage_plot.{suffix}"), 
            category="QC", subcategory="Coverage",  labels=lambda wildcards: {
              "Sample": "{sample}",
              "Region": "CDS" if "CDS" in wildcards.suffix else "Whole transcript",
              "figure": "flame_coverage_plot"
          })
    resources:
        cpus_per_task=1,
        mem_mb=32000
        #slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
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
        reads=os.path.join(results_dir, "qc/subsample_fq/{sample}_{cell_line}_subsampled1M.fastq"),
        gtf=config['reference']['gtf'],
        genome=config['reference']['genome']
    output:
        dir = directory(os.path.join(results_dir, "qc/sqanti3/{sample}_{cell_line}/"))
        # html = report(os.path.join(results_dir, "qc/sqanti3/{sample}/matched_reads_subsampled_SQANTI3_report.html"),
        #             category="QC", subcategory="SQANTI3") # this is not included as it is too large
    conda: 
        config['conda']['sqanti3']
    resources:
        cpus_per_task=16,
        mem_mb=32000
        #slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params: 
        sqanti3_qc_script = os.path.join(config['software']['sqanti3_dir'], "sqanti3_qc.py"),
        additional_arg = "--skipORF"
    shell:
        """
        module unload R
        {params.sqanti3_qc_script} {input.reads} {input.gtf} {input.genome} {params.additional_arg} --fasta  --cpus {resources.cpus_per_task} --force_id_ignore -d {output} --report html
        """



rule RSeQC_gene_body_coverage:
    input:
        bed=config['reference']['bed_human'], # I have rules to convert gtf to bed
        genome_bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
    output:
        os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.RSeQC.geneBodyCoverage.pdf")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    conda:
        config['conda']['RSeQC']
    shell:
        """
        mkdir -p $(dirname {output})
        geneBody_coverage.py -i {input.genome_bam} -r {input.bed} -o {output}
        """



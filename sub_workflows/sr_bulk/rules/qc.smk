from os.path import join

# Step 1: Trim reads using fastp
rule fastp:
    input:
        R1 = join(input_fastq_dirs['ill_bulk'], "{cell_line}_R1.fastq.gz"),
        R2 = join(input_fastq_dirs['ill_bulk'], "{cell_line}_R2.fastq.gz")
    output:
        R1=os.path.join(scratch_dir, "trimmed_fq/{cell_line}_R1.fastq.gz"),
        R2=os.path.join(scratch_dir, "trimmed_fq/{cell_line}_R2.fastq.gz"),
        html=join(results_dir, "qc/fastp/{cell_line}.html"),
        json=join(results_dir, "qc/fastp/{cell_line}.json")
    conda:
        config['conda']['fastp']
    resources:
        cpus_per_task=16,
        mem_mb=32000
    shell:
        """
        mkdir -p $(dirname {output[0]})
        fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2}  --thread {resources.cpus_per_task} -j {output.json} -h {output.html}
        """


# Entire qc worflow head
rule qc:
    input:
        expand(
            [
                rules.fastp.output.R1,
                rules.fastp.output.R2
            ],
            cell_line = config['cell_lines']
        )
    output:
        touch(join(results_dir, "qc/.flag/qc.done"))

###################### Run QC pipelines ###############################
# rule qc:
#     input:
#         expand(
#             [
#                 os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.flame.coverage_plot.{flames_cov_plot_suffix}"),
#                 # os.path.join(results_dir, "qc/sqanti3/{sample}_{cell_line}/"),
#                 os.path.join(results_dir, "qc/NanoPlot/{sample}_{cell_line}/NanoPlot-data.tsv.gz"),
#                 os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.picard.RNA_Metrics"),
#                 os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.geneBodyCoverage.curves.pdf"),
#                 os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.junctionSaturation_plot.pdf"),
#                 os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.splice_events.pdf"),
#                 os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.splice_junction.pdf"),
#                 os.path.join(results_dir, "qc/aligment_summary/{alignment_type}/{sample}_{cell_line}.picard_AlignmentSummaryMetrics.txt"),
#             ],
#             sample=config["sample_id"],
#             flames_cov_plot_suffix = ['png'],
#             cell_line = cell_line_to_barcode.keys(),
#             alignment_type = ['GenomeAlignment', 'TranscriptAlignment']
#         )
#     output:
#         os.path.join(results_dir, "qc/.flag/qc.done")
#     shell:
#         """
#         mkdir -p $(dirname {output})
#         touch {output}
#         """
###################################################################################

######## Subsample ########
# rule subsample_2M_reads:
#     input:
#         lambda wildcards:  glob.glob(os.path.join(input_fastq_dirs[wildcards.sample], f"{wildcards.cell_line}.fastq*"))[0]
#     output:
#         os.path.join(scratch_dir, "subsample_fq/raw/{sample}_{cell_line}_subsampled2M.fastq")
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000
#     params:
#         n_reads = 2000000,
#         seed = config['random_seed']
#     shell:
#         """
#         mkdir -p $(dirname {output})
#         seqtk sample -s {params.seed} {input} {params.n_reads} > {output}
#         """

# ######## Count reads ########
# rule count_reads_in_fastq:
#     input:
#         lambda wildcards:  glob.glob(os.path.join(input_fastq_dirs[wildcards.sample], f"{wildcards.cell_line}.fastq*"))[0]
#     output:
#         os.path.join(results_dir, "qc/read_counts/{sample}_{cell_line}.count")
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000
#     shell:
#         """
#         mkdir -p $(dirname {output})
#         wc -l {input} | awk '{{print $1/4}}' > {output}
#         """
# 
# ######## NanoPlot ########
# rule NanoPlot:
#     input:
#         reads=rules.subsample_2M_reads.output
#     output:
#         os.path.join(results_dir, "qc/NanoPlot/{sample}_{cell_line}/NanoPlot-data.tsv.gz")
#     conda:
#         config['conda']['NanoPlot']
#     resources:
#         cpus_per_task=16,
#         mem_mb=32000
#     shell:
#         """
#         output_dir=$(dirname {output})
#         mkdir -p $output_dir
#         NanoPlot --fastq {input.reads} --outdir $output_dir -t {resources.cpus_per_task} --raw  --tsv_stats
#         """
# 
# # Coverage
# rule picard_coverage_data:
#     input:
#         bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
#         refFlat = '/home/users/allstaff/you.yu/LongBench/reference_files/GRCh38/refFlat.txt',
#         picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar"
#     output:
#         os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.picard.RNA_Metrics")
#     resources:
#         cpus_per_task=16,
#         mem_mb=32000,
#         slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
#     shell:
#         """
#         module load picard-tools
#         mkdir -p $(dirname {output})
#         java -jar {input.picard_dir} \
#             CollectRnaSeqMetrics  \
#             I={input.bam} O={output} \
#             REF_FLAT={input.refFlat} \
#             STRAND_SPECIFICITY=NONE
#         """
# 
# # BamIndexStats
# rule picard_CollectAlignmentSummaryMetrics:
#     input:
#         bam = os.path.join(results_dir,"{alignment_type}/{sample}_{cell_line}.sorted.bam"),
#         picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar",
#         ref = lambda w: config['reference']['genome' if w.alignment_type == 'GenomeAlignment' else 'transcript']
#     output:
#         os.path.join(results_dir, "qc/aligment_summary/{alignment_type}/{sample}_{cell_line}.picard_AlignmentSummaryMetrics.txt")
#     resources:
#         cpus_per_task=16,
#         mem_mb=32000,
#         slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
#     shell:
#         """
#         module load picard-tools
#         mkdir -p $(dirname {output})
#         java -jar {input.picard_dir} CollectAlignmentSummaryMetrics \
#           R={input.ref} \
#           I={input.bam} \
#           O={output}
#         """
# 
# rule flame_coverage_plot:
#     input:
#         bam = os.path.join(results_dir,"TranscriptAlignment/{sample}_{cell_line}.sorted.bam"),
#         gtf = config['reference']['gtf']
#     output:
#         report(
#             os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.flame.coverage_plot.{suffix}"), 
#             category="QC", subcategory="Coverage",  labels=lambda wildcards: {
#               "Sample": "{sample}",
#               "Region": "CDS" if "CDS" in wildcards.suffix else "Whole transcript",
#               "figure": "flame_coverage_plot"
#           })
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000
#         #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
#     params:
#         script = os.path.join(config['main_wf_dir'],'scripts/flames_coverage_plot.R')
#     shell:
#         """
#         mkdir -p $(dirname {output})
#         Rscript  {params.script} {input.bam} {input.gtf} {output}
#         """
# 
# 
# ## SQANTI3
# rule sqanti3:
#     input:
#         reads=os.path.join(results_dir, "qc/subsample_fq/{sample}_{cell_line}_subsampled1M.fastq"),
#         gtf=config['reference']['gtf'],
#         genome=config['reference']['genome']
#     output:
#         dir = directory(os.path.join(results_dir, "qc/sqanti3/{sample}_{cell_line}/"))
#         # html = report(os.path.join(results_dir, "qc/sqanti3/{sample}/matched_reads_subsampled_SQANTI3_report.html"),
#         #             category="QC", subcategory="SQANTI3") # this is not included as it is too large
#     conda: 
#         config['conda']['sqanti3']
#     resources:
#         cpus_per_task=16,
#         mem_mb=32000
#         #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
#     params: 
#         sqanti3_qc_script = os.path.join(config['software']['sqanti3_dir'], "sqanti3_qc.py"),
#         additional_arg = "--skipORF"
#     shell:
#         """
#         {params.sqanti3_qc_script} {input.reads} {input.gtf} {input.genome} {params.additional_arg} --fasta  --cpus {resources.cpus_per_task} --force_id_ignore -d {output} --report html
#         """
# 
# 
# 
# rule Subsample_bam_for_RSeQC:
#     input:
#         bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
#     output:
#         bam = temp(os.path.join(results_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_rate_{subsample_rate}.bam")),
#         bai = temp(os.path.join(results_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_rate_{subsample_rate}.bam.bai"))
#     resources:
#         cpus_per_task=4,
#         mem_mb=64000
#     params:
#         seed = config['random_seed'],
#     shell:
#         """
#         mkdir -p $(dirname {output.bam})
#         samtools view --subsample {wildcards.subsample_rate} --subsample-seed {params.seed} {input.bam} -h | samtools sort -o {output.bam}
#         samtools index {output.bam}
#         """
# 
# 
# rule RSeQC_gene_body_coverage:
#     input:
#         bed=config['reference']['bed_housekeeping_genes'], 
#         genome_bam = os.path.join(results_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_rate_0.01.bam")
#     output:
#         report(os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.geneBodyCoverage.curves.pdf"))
#     resources:
#         cpus_per_task=4,
#         mem_mb=64000
#     conda:
#         config['conda']['RSeQC']
#     shell:
#         """
#         mkdir -p $(dirname {output})
#         geneBody_coverage.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}_{wildcards.cell_line}
#         """
# 
# 
# rule RSeQC_junction_saturation:
#     input:
#         bed=config['reference']['bed_human'], 
#         genome_bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
#     output:
#         report(os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.junctionSaturation_plot.pdf"))
#     resources:
#         cpus_per_task=1,
#         mem_mb=64000
#     conda:
#         config['conda']['RSeQC']
#     shell:
#         """
#         mkdir -p $(dirname {output})
#         junction_saturation.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}_{wildcards.cell_line}
#         """
# 
# 
# rule RSeQC_junction_annotation:
#     input:
#         bed=config['reference']['bed_human'], # I have rules to convert gtf to bed
#         genome_bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
#     output:
#         report(os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.splice_events.pdf")),
#         report(os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.splice_junction.pdf"))
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000
#     conda:
#         config['conda']['RSeQC']
#     shell:
#         """
#         mkdir -p $(dirname {output[0]})
#         junction_annotation.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}_{wildcards.cell_line}
#         """
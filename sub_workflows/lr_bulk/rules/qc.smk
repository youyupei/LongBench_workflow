import glob

cell_line_to_barcode = {cl: bc for d in barcode_list for bc, cl in d.items()}
###################################################################################

# ######## Qscore filtering (chopper) ########
# rule chopper_qscore_filtering:
#     input:
#         lambda wildcards:  os.path.join(input_fastq_dirs[wildcards.sample], "{cell_line}.fastq")
#     output:
#         os.path.join(scratch_dir, "Q{q_threshold}_filter_fastqs/{sample}_{cell_line}.fastq")
#     resources:
#         cpus_per_task=8,
#         mem_mb=32000
#     conda:
#         config['conda']['NanoPlot']
#     shell:
#         """
#         mkdir -p $(dirname {output})
# 
#         # if a is 0, make a soft link from input to output, else run chopper
#         if [ {wildcards.q_threshold} -eq 0 ]; then
#             ln -s {input} {output}
#         else
#             chopper -q {wildcards.q_threshold} -i {input} -o {output} --threads {resources.cpus_per_task}
#         fi
#         """


######## Subsample ########
rule subsample_2M_reads:
    input:
        lambda wildcards:  glob.glob(os.path.join(input_fastq_dirs[wildcards.sample], f"{wildcards.cell_line}.fastq*"))[0]
    output:
        os.path.join(scratch_dir, "subsample_fq/raw/{sample}_{cell_line}_subsampled2M.fastq")
    resources:
        cpus_per_task=1,
        mem_mb=24000
    params:
        n_reads = 2000000,
        seed = config['random_seed']
    shell:
        """
        mkdir -p $(dirname {output})
        seqtk sample -s {params.seed} {input} {params.n_reads} > {output}
        """


rule Subsample_bam_for_qc:
    input:
        results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
    output:
        bam = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_{n_reads}.bam"),
        bai = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_{n_reads}.bam.bai")
    resources:
        cpus_per_task=1,
        mem_mb=8000
    params:
        seed = config['random_seed'],
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view {input} | cut -f1 | sort | uniq > {output.bam}.read_ids
        shuf -n $(numfmt --from=si  {wildcards.n_reads}) --random-source=<(yes {params.seed}) {output.bam}.read_ids  > {output.bam}.read_ids.subsampled
        rm {output.bam}.read_ids
        samtools view -b -N {output.bam}.read_ids.subsampled {input} > {output.bam}
        samtools index {output.bam}
        """

######## Count reads ########
rule count_reads_in_fastq:
    input:
        lambda wildcards:  glob.glob(os.path.join(input_fastq_dirs[wildcards.sample], f"{wildcards.cell_line}.fastq*"))[0]
    output:
        os.path.join(results_dir, "qc/read_counts/{sample}_{cell_line}.count")
    resources:
        cpus_per_task=1,
        mem_mb=8000
    shell:
        """
        mkdir -p $(dirname {output})
        wc -l {input} | awk '{{print $1/4}}' > {output}
        """

######## NanoPlot ########
rule NanoPlot:
    input:
        reads=rules.subsample_2M_reads.output
    output:
        os.path.join(results_dir, "qc/NanoPlot/{sample}_{cell_line}/NanoPlot-data.tsv.gz")
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

######## Coverage ########
rule picard_coverage_data:
    input:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
        refFlat = '/home/users/allstaff/you.yu/LongBench/reference_files/GRCh38/refFlat.txt',
        picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar"
    output:
        os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.picard.RNA_Metrics")
    resources:
        cpus_per_task=16,
        mem_mb=16000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
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
        bam = os.path.join(results_dir,"{alignment_type}/{sample}_{cell_line}.sorted.bam"),
        picard_dir = "/home/users/allstaff/you.yu/project/software/picard.jar",
        ref = lambda w: config['reference']['genome' if w.alignment_type == 'GenomeAlignment' else 'transcript']
    output:
        os.path.join(results_dir, "qc/aligment_summary/{alignment_type}/{sample}_{cell_line}.picard_AlignmentSummaryMetrics.txt")
    resources:
        cpus_per_task=16,
        mem_mb=16000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
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
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        script = os.path.join(config['main_wf_dir'],'scripts/flames_coverage_plot.R')
    shell:
        """
        mkdir -p $(dirname {output})
        Rscript  {params.script} {input.bam} {input.gtf} {output}
        """


# # SQANTI3
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


## RSeQC
rule RSeQC_gene_body_coverage:
    input:
        bed=config['reference']['bed_human'], 
        genome_bam = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_3M.bam"),
        genome_bai = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_3M.bam.bai")
    output:
        report(os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.geneBodyCoverage.curves.pdf")),
        os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.geneBodyCoverage.r")
    resources:
        cpus_per_task=4,
        mem_mb=16000
    conda:
        config['conda']['RSeQC']
    shell:
        """
        mkdir -p $(dirname {output})
        geneBody_coverage.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}_{wildcards.cell_line}
        """


rule RSeQC_junction_saturation:
    input:
        bed=config['reference']['bed_human'], 
        genome_bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
    output:
        report(os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.junctionSaturation_plot.pdf")),
        os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.junctionSaturation_plot.r")
    resources:
        cpus_per_task=8,
        mem_mb=64000
    conda:
        config['conda']['RSeQC']
    shell:
        """
        mkdir -p $(dirname {output})
        junction_saturation.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}_{wildcards.cell_line}
        """


rule RSeQC_junction_annotation:
    input:
        bed=config['reference']['bed_human'], # I have rules to convert gtf to bed
        genome_bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
    output:
        report(os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.splice_events.pdf")),
        report(os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.splice_junction.pdf"))
    resources:
        cpus_per_task=1,
        mem_mb=8000
    conda:
        config['conda']['RSeQC']
    shell:
        """
        mkdir -p $(dirname {output[0]})
        junction_annotation.py -i {input.genome_bam} -r {input.bed} -o $(dirname {output[0]})/{wildcards.sample}_{wildcards.cell_line}
        """



## AlignQC
rule alignQC_samflag_filter:
    input:
        genome_bam = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_3M.bam"),
        genome_bai = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_3M.bam.bai")
    output:
        os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_3M.F2304.bam")
    resources:
        cpus_per_task=1,
        mem_mb=8000
    shell:
        """
        mkdir -p $(dirname {output})
        samtools view -b -F 2304 {input.genome_bam} > {output}
        samtools index {output}

        """

rule alignQC_analysis_subsample:
    input: 
        fa=config['reference']['genome'], 
        anno=config['reference']['gtf_gz'],
        genome_bam = rules.alignQC_samflag_filter.output
    resources:
        cpus_per_task=12,
        mem_mb=500000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au" #
    output:
        output = directory(os.path.join(results_dir, "qc/AlignQC/{sample}_{cell_line}/")),
        tmp_dir = temp(directory(os.path.join(scratch_dir, "alignQC_tmp","{sample}_{cell_line}")))
    priority: 10
    container: "docker://vacation/alignqc"
    #retries: 3
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
# TOP aligned length
rule find_top_aligned_length:
    input:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam"
        #bai = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.sorted.bam.bai"
    output:
        os.path.join(results_dir, "qc/TopAlignedRead/{sample}_{cell_line}.tsv")
    resources:
        cpus_per_task=1,
        mem_mb=8000
    params:
        script = os.path.join(config['main_wf_dir'],'scripts/find_longest_reads_in_bam.py')
    shell:
        """
        mkdir -p $(dirname {output})
        python3 {params.script} {input.bam} {output} -n 100
        """


# use rule alignQC_analysis_subsample as alignQC_analysis_full with:
#     input: 
#         fa=config['reference']['genome'], 
#         anno=config['reference']['gtf_gz'],
#         genome_bam = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_1.bam"),
#         genome_bai = os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_subsample_1.bam.bai")
#     output:
#         output = directory(os.path.join(results_dir, "qc/AlignQC_full/{sample}_{cell_line}/")),
#         tmp_dir = temp(directory(os.path.join(scratch_dir, "alignQC_full_tmp","{sample}_{cell_line}")))
#     resources:
#         cpus_per_task=8,
#         mem_mb=400000,
#         slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au --qos=bonus" # tmp use bonus qos

### Target QC rules
###################### Run QC pipelines ###############################
rule qc:
    input:
        expand(
            [
                # os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.flame.coverage_plot.{flames_cov_plot_suffix}"), # flames coverage plot
                # os.path.join(results_dir, "qc/sqanti3/{sample}_{cell_line}/"), # SQANTI3
                os.path.join(results_dir, "qc/AlignQC/{sample}_{cell_line}/"), # AlignQC
                os.path.join(results_dir, "qc/NanoPlot/{sample}_{cell_line}/NanoPlot-data.tsv.gz"), # NanoPlot
                # os.path.join(results_dir, "qc/coverage/{sample}_{cell_line}.picard.RNA_Metrics"), # picard coverage
                # os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.geneBodyCoverage.curves.pdf"), # RSeQC gene body coverage
                os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.junctionSaturation_plot.pdf"), # RSeQC junction saturation
                os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.splice_events.pdf"), # RSeQC_junction_annotation
                os.path.join(results_dir, "qc/RSeQC/{sample}_{cell_line}.splice_junction.pdf"), # RSeQC_junction_annotation
                # os.path.join(results_dir, "qc/aligment_summary/{alignment_type}/{sample}_{cell_line}.picard_AlignmentSummaryMetrics.txt"), # picard alignment summary
                os.path.join(results_dir, "qc/TopAlignedRead/{sample}_{cell_line}.tsv")
            ],
            sample=config["sample_id"],
            flames_cov_plot_suffix = ['png'],
            cell_line = cell_line_to_barcode.keys(),
            alignment_type = ['GenomeAlignment', 'TranscriptAlignment']
        )
    output:
        os.path.join(results_dir, "qc/.flag/qc.done")
    shell:
        """
        mkdir -p $(dirname {output})
        touch {output}
        """


rule all_alignQC_analysis_subsample:
    input:
        expand(
            os.path.join(results_dir, "qc/AlignQC/{sample}_{cell_line}/"), 
            sample=config["sample_id"],
            cell_line = cell_line_to_barcode.keys()
        )
    output:
        touch(os.path.join(config['flag_dir'], "lr_bulk_alginQC.done"))

# rule all_alignQC_analysis_full:
#     input:
#         expand(
#             rules.alignQC_analysis_full.output.output, 
#             sample=config["sample_id"],
#             cell_line = cell_line_to_barcode.keys()
#         )
#     output:
#         touch(os.path.join(config['flag_dir'], "lr_bulk_alginQC_full.done"))

# rule all_sqanti3:
#     input:
#         expand(
#             rules.sqanti3.output.dir, 
#             sample=config["sample_id"],
#             cell_line = cell_line_to_barcode.keys()
#         )
#     output:
#         touch(os.path.join(config['flag_dir'], "lr_bulk_sqanti3.done"))
figures_output_path = config["figures_output_path"]

lr_sc_sn_result_dir = sub_wf_config['lr_sc_sn']["output_path"]
lr_bulk_result_dir = sub_wf_config['lr_bulk']["output_path"]
sr_bulk_result_dir = sub_wf_config['sr_bulk']["output_path"]

from os.path import join

# Total number
rule read_number_plot:
    input:
        bulk_read_count = expand(
            join(lr_bulk_result_dir, "qc/read_counts/{sample}_{cell_line}.count"),
            sample = sub_wf_config['lr_bulk']["sample_id"],
            cell_line = sub_wf_config['lr_bulk']["cell_lines"]
        ),
        sc_blaze_summary = expand(
            join(lr_sc_sn_result_dir, "flames_out/{sample}/summary.txt"),
            sample = sub_wf_config['lr_sc_sn']["sample_id"]
        ),
        sr_bulk_fastp_json = expand(
            join(sr_bulk_result_dir, "qc/fastp/{cell_line}.json"),
            cell_line = sub_wf_config['lr_bulk']["cell_lines"]
        )
    output:
        report(join(figures_output_path, "qc/read_number_plot.pdf"), 
                category = "QC", subcategory = "Read number")
    params:
        bulk_sample_name = expand("{sample}_{cell_line}", 
                            sample = sub_wf_config['lr_bulk']["sample_id"],
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"]),
        sc_sample_name = sub_wf_config['lr_sc_sn']["sample_id"],
        sr_sample_name = expand("Illumina_{cell_line}", 
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"])
    script:
        join(config['main_wf_dir'], "scripts/read_count_plot.R")

    


rule lr_read_length_plot:
    # get the read length from the read length file and plot
    input: 
        expand(
            join(lr_bulk_result_dir, "qc/NanoPlot/{sample}_{cell_line}/NanoPlot-data.tsv.gz"),
            sample = sub_wf_config['lr_bulk']["sample_id"],
            cell_line = sub_wf_config['lr_bulk']["cell_lines"]
        ),
        expand(
            join(lr_sc_sn_result_dir, "qc/NanoPlot/{sample}/NanoPlot-data.tsv.gz"),
            sample = sub_wf_config['lr_sc_sn']["sample_id"]
        ),
        # Not looking good when including the sr_bulk data
        # expand(
        #     join(sr_bulk_result_dir, "qc/NanoPlot/{cell_line}/NanoPlot-data.tsv.gz"),
        #     sample = sub_wf_config['lr_bulk']["sample_id"],
        #     cell_line = sub_wf_config['lr_bulk']["cell_lines"]
        # )
    output:
        report(join(figures_output_path, "qc/read_length_and_quality_plot.pdf"), 
                category = "QC", subcategory = "Read length and quality")
    params:
        sample_id = expand("{sample}_{cell_line}", 
                            sample = sub_wf_config['lr_bulk']["sample_id"],
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"]) + \
                         sub_wf_config['lr_sc_sn']["sample_id"]
    resources:
        cpus_per_task=1,
        mem_mb=1000

    script:
        join(config['main_wf_dir'], "scripts/read_length_and_quality_plot.R")

rule RSeQC_gene_body_coverage_plot:
    input:
        lr_bulk = \
            expand(
                join(lr_bulk_result_dir, "qc/RSeQC/{sample}_{cell_line}.geneBodyCoverage.r"),
                sample = sub_wf_config['lr_bulk']["sample_id"],
                cell_line = sub_wf_config['lr_bulk']["cell_lines"]
            ),
        sr_bulk = \
            expand(
                join(sr_bulk_result_dir, "qc/RSeQC/{cell_line}.geneBodyCoverage.r"),
                cell_line = sub_wf_config['lr_bulk']["cell_lines"]
            )
    output:
        report(join(figures_output_path, "qc/gene_body_coverage_plot.pdf"), 
                category = "QC", subcategory = "Gene body coverage")
    params:
        lr_sample_name = expand("{sample}_{cell_line}", 
                            sample = sub_wf_config['lr_bulk']["sample_id"],
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"]),
        sr_sample_name = expand("Illumina_{cell_line}", 
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"])
    resources:
        cpus_per_task=1,
        mem_mb=1000
    script:
        join(config['main_wf_dir'], "scripts/genebody_coverage_plot.R")


rule RSeQC_junction_saturation_plot_known:
    input:
        lr_bulk = \
            expand(
                join(lr_bulk_result_dir, "qc/RSeQC/{sample}_{cell_line}.junctionSaturation_plot.r"),
                sample = sub_wf_config['lr_bulk']["sample_id"],
                cell_line = sub_wf_config['lr_bulk']["cell_lines"]
            ),
        sr_bulk = \
            expand(
                join(sr_bulk_result_dir, "qc/RSeQC/{cell_line}.junctionSaturation_plot.r"),
                cell_line = sub_wf_config['sr_bulk']["cell_lines"]
            )
    output:
        report(join(figures_output_path, "qc/KnownJunctionSaturation_plot.pdf"), 
                category = "QC", subcategory = "Gene body coverage")
    params:
        lr_sample_name = expand("{sample}_{cell_line}", 
                            sample = sub_wf_config['lr_bulk']["sample_id"],
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"]),
        sr_sample_name = expand("Illumina_{cell_line}", 
                            cell_line = sub_wf_config['sr_bulk']["cell_lines"])
    resources:
        cpus_per_task=1,
        mem_mb=1000
    script:
        join(config['main_wf_dir'], "scripts/junctionSaturation_plot_known.R")

rule RSeQC_junction_saturation_plot_novel:
    input:
        lr_bulk = \
            expand(
                join(lr_bulk_result_dir, "qc/RSeQC/{sample}_{cell_line}.junctionSaturation_plot.r"),
                sample = sub_wf_config['lr_bulk']["sample_id"],
                cell_line = sub_wf_config['lr_bulk']["cell_lines"]
            ),
        sr_bulk = \
            expand(
                join(sr_bulk_result_dir, "qc/RSeQC/{cell_line}.junctionSaturation_plot.r"),
                cell_line = sub_wf_config['sr_bulk']["cell_lines"]
            )
    output:
        report(join(figures_output_path, "qc/NovelJunctionSaturation_plot.pdf"), 
                category = "QC", subcategory = "Gene body coverage")
    params:
        lr_sample_name = expand("{sample}_{cell_line}", 
                            sample = sub_wf_config['lr_bulk']["sample_id"],
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"]),
        sr_sample_name = expand("Illumina_{cell_line}", 
                            cell_line = sub_wf_config['sr_bulk']["cell_lines"])
    resources:
        cpus_per_task=1,
        mem_mb=1000
    script:
        join(config['main_wf_dir'], "scripts/junctionSaturation_plot_novel.R")
# 
# rule summarise_sqanti:
#     input: 
#         expand(join(lr_bulk_result_dir, "qc/sqanti3/{sample}_{cell_line}/"),
#             sample = sub_wf_config['lr_bulk']["sample_id"],
#             cell_line = sub_wf_config['lr_bulk']["cell_lines"]),
#         expand(join(lr_sc_sn_result_dir, "qc/sqanti3/{sample}/"),
#                 sample = sub_wf_config['lr_sc_sn']["sample_id"])
#     output:
#         report(join(figures_output_path, "qc/sqanti_summary.pdf"), 
#                 category = "QC", subcategory = "SQANTI3")
#     params:
#         sample_names = expand("{sample}_{cell_line}", 
#                             sample = sub_wf_config['lr_bulk']["sample_id"],
#                             cell_line = sub_wf_config['lr_bulk']["cell_lines"]) + \
#                          sub_wf_config['lr_sc_sn']["sample_id"],
#         utilities_path = join(config['software']['sqanti3_dir'],
#                                         "utilities")
#     script:
#         join(config['main_wf_dir'], "scripts/summarise_sqanti.R")



rule main_qc_plot:
    input:
        rules.lr_read_length_plot.output,
        #rules.summarise_sqanti.output,
        rules.read_number_plot.output,
        rules.RSeQC_gene_body_coverage_plot.output,
        rules.RSeQC_junction_saturation_plot_known.output,
        rules.RSeQC_junction_saturation_plot_novel.output
    output:
        touch(join(figures_output_path, "qc/.flag.qc_plot.done"))
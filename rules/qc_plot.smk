combined_result_dir = config["combined_output_path"]

lr_sc_sn_result_dir = sub_wf_config['lr_sc_sn']["output_path"]
lr_bulk_result_dir = sub_wf_config['lr_bulk']["output_path"]

# Total number
rule read_number_plot:
    input:
        bulk_read_count = expand(
            os.path.join(lr_bulk_result_dir, "qc/read_counts/{sample}_{cell_line}.count"),
            sample = sub_wf_config['lr_bulk']["sample_id"],
            cell_line = sub_wf_config['lr_bulk']["cell_lines"]
        ),
        sc_blaze_summary = expand(
            os.path.join(lr_sc_sn_result_dir, "flames_out/{sample}/summary.txt"),
            sample = sub_wf_config['lr_sc_sn']["sample_id"]
        )
    output:
        report(os.path.join(combined_result_dir, "qc/read_number_plot.pdf"), 
                category = "QC", subcategory = "Read number")
    params:
        bulk_sample_name = expand("{sample}_{cell_line}", 
                            sample = sub_wf_config['lr_bulk']["sample_id"],
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"]),
        sc_sample_name = sub_wf_config['lr_sc_sn']["sample_id"]
    script:
        os.path.join(config['main_wf_dir'], "scripts/read_count_plot.R")

    


rule lr_read_length_plot:
    # get the read length from the read length file and plot
    input: 
        expand(
            os.path.join(lr_bulk_result_dir, "qc/NanoPlot/{sample}_{cell_line}/NanoPlot-data.tsv.gz"),
            sample = sub_wf_config['lr_bulk']["sample_id"],
            cell_line = sub_wf_config['lr_bulk']["cell_lines"]
        ),
        expand(
            os.path.join(lr_sc_sn_result_dir, "qc/NanoPlot/{sample}/NanoPlot-data.tsv.gz"),
            sample = sub_wf_config['lr_sc_sn']["sample_id"]
        )
    output:
        report(os.path.join(combined_result_dir, "qc/read_length_and_quality_plot.pdf"), 
                category = "QC", subcategory = "Read length and quality")
    params:
        sample_id = expand("{sample}_{cell_line}", 
                            sample = sub_wf_config['lr_bulk']["sample_id"],
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"]) + \
                         sub_wf_config['lr_sc_sn']["sample_id"]
    resources:
        cpus_per_task=2,
        mem_mb=64000

    script:
        os.path.join(config['main_wf_dir'], "scripts/read_length_and_quality_plot.R")


rule summarise_sqanti:
    input: 
        expand(os.path.join(lr_bulk_result_dir, "qc/sqanti3/{sample}_{cell_line}/"),
            sample = sub_wf_config['lr_bulk']["sample_id"],
            cell_line = sub_wf_config['lr_bulk']["cell_lines"]),
        expand(os.path.join(lr_sc_sn_result_dir, "qc/sqanti3/{sample}/"),
                sample = sub_wf_config['lr_sc_sn']["sample_id"])
    output:
        report(os.path.join(combined_result_dir, "qc/sqanti_summary.pdf"), 
                category = "QC", subcategory = "SQANTI3")
    params:
        sample_names = expand("{sample}_{cell_line}", 
                            sample = sub_wf_config['lr_bulk']["sample_id"],
                            cell_line = sub_wf_config['lr_bulk']["cell_lines"]) + \
                         sub_wf_config['lr_sc_sn']["sample_id"],
        utilities_path = os.path.join(config['software']['sqanti3_dir'],
                                        "utilities")
    script:
        os.path.join(config['main_wf_dir'], "scripts/summarise_sqanti.R")

rule combined_qc_plot:
    input:
        rules.lr_read_length_plot.output,
        #rules.summarise_sqanti.output,
        rules.read_number_plot.output
    output:
        touch(os.path.join(combined_result_dir, "qc/combined_qc_plot.done"))
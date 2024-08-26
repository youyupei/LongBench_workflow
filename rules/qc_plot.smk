combined_result_dir = config["combined_output_path"]

lr_sc_sn_result_dir = sub_wf_config['lr_sc_sn']["output_path"]
lr_bulk_result_dir = sub_wf_config['lr_bulk']["output_path"]


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
    script:
        os.path.join(config['main_wf_dir'], "scripts/read_length_and_quality_plot.R")

# 
# rule summarise_sqanti:
#     input: 
#         expand(os.path.join(results_dir, "qc/sqanti3/{sample}_{cell_line}/"),
#             sample = sub_wf_config['lr_bulk']["sample_id"],
#             cell_line = sub_wf_config['lr_bulk']["cell_lines"]),
#         expand(os.path.join(results_dir, "qc/sqanti3/{sample}/"),
#                 sample = sub_wf_config['lr_sc_sn']["sample_id"])
#     output:
#         report(os.path.join(combined_result_dir, "qc/sqanti_summary.pdf"), 
#                 category = "QC", subcategory = "SQANTI3")
#     params:
#         sample_id = expand("{sample}_{cell_line}", 
#                             sample = sub_wf_config['lr_bulk']["sample_id"],
#                             cell_line = sub_wf_config['lr_bulk']["cell_lines"]) + \
#                          sub_wf_config['lr_sc_sn']["sample_id"],
#         utilities_path = os.path.join(config['software']['sqanti3_dir'],
#                                         "utilities")
#     script:
#         os.path.join(config['main_wf_dir'], "scripts/summarise_sqanti.R")
# # working directory: config['main_wf_dir']
# from os.path import join
# 
# rule rmd_bulk_de_human:
#     input:
#         rmd = join(config['main_wf_dir'], "rmarkdown/Bulk_DE_analysis_human.Rmd"),
#         sirv_ercc_gtf= '/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_ERCC_longSIRV_multi-fasta_20210507.gtf',
#         sequins_gtf= '/vast/projects/LongBench/reference_files/Sequin_resources/v2.4/rnasequin_annotation_2.4.gtf',
#         sequins_tsv= '/vast/projects/LongBench/reference_files/Sequin_resources/v2.4/rnasequin_isoforms_2.4.tsv',
#         human_gtf= '/vast/projects/LongBench/reference_files/GRCh38/gencode.v44.annotation.gtf',
#         sirv_csv= '/home/users/allstaff/you.yu/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_set4_concentration.csv',
#         ont_bulk_oarfish_dir= '/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/ont_bulk',
#         pb_bulk_oarfish_dir= '/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/pb_bulk',
#         ont_drnd_oarfish_dir= '/vast/projects/LongBench/analysis/lr_bulk/result/oarfish_cov_output/dRNA_bulk',
#         ill_bulk_salmon_dir= '/vast/projects/LongBench/analysis/sr_bulk/result/salmon/salmon_quant',
#         bulk_meta= '/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt'
#     output:
#         figures = directory(join(config['figures_output_path'], "bulk_de_human/")),
#         html_report = join(config['reports_output_path'], "bulk_de_human.html"),
#         rst_rds = join(config['reports_output_path'], "RDS/bulk_de_human_summary.rds"),
#         cache_dir = directory(join(scratch_dir, ".rmd_Cache/bulk_de_human/"))
#     params:
#         random_seed = 2024
#     resources:
#         mem_mb = 16000,
#         cpus_per_task = 1
#     shell:
#         """
#         mkdir -p $(dirname {output.html_report})
#         mkdir -p {output.cache_dir}
#         mkdir -p {output.figures}
#         mkdir -p $(dirname {output.rst_rds}/)
#         Rscript -e 'rmarkdown::render("{input.rmd}", 
#             output_file="{output.html_report}", 
#             params=list(
#                 sirv_ercc_gtf="{input.sirv_ercc_gtf}",
#                 sequins_gtf="{input.sequins_gtf}",
#                 sequins_tsv="{input.sequins_tsv}",
#                 human_gtf="{input.human_gtf}",
#                 sirv_csv="{input.sirv_csv}",
#                 ont_bulk_oarfish_dir="{input.ont_bulk_oarfish_dir}",
#                 pb_bulk_oarfish_dir="{input.pb_bulk_oarfish_dir}",
#                 ont_drnd_oarfish_dir="{input.ont_drnd_oarfish_dir}",
#                 ill_bulk_salmon_dir="{input.ill_bulk_salmon_dir}",
#                 bulk_meta="{input.bulk_meta}",
#                 cache_dir="{output.cache_dir}/",
#                 random_seed={params.random_seed},
#                 fig.path="{output.figures}/",
#                 out.rds="{output.rst_rds}"
#             )
#         )'
#         """
# 
# rule rmd_sc_sn_de:
#     input:
#         rmd = join(config['main_wf_dir'], "rmarkdown/SC_DE_analysis.Rmd"),
#         bulk_rds = rules.rmd_bulk_de_human.output.rst_rds
#     output:
#         figures = directory(join(config['figures_output_path'], "sc_sn_de/")),
#         html_report = join(config['reports_output_path'], "sc_sn_de.html"),
#         cache_dir = directory(join(scratch_dir, ".rmd_Cache/sc_sn_de/"))
#     params:
#         random_seed = 2024
#     resources:
#         mem_mb = 16000,
#         cpus_per_task = 1
#     shell:
#         """
#         mkdir -p $(dirname {output.html_report})
#         mkdir -p {output.cache_dir}
#         mkdir -p {output.figures}
#         Rscript -e 'rmarkdown::render("{input.rmd}", 
#             output_file="{output.html_report}", 
#             params=list(
#                 random_seed={params.random_seed},
#                 cache_dir="{output.cache_dir}/",
#                 ont_sc_dir= "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sc",
#                 ont_sn_dir= "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sn",
#                 pb_sc_dir= "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sc",
#                 pb_sn_dir= "/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sn",
#                 human_gtf= "/vast/projects/LongBench/reference_files/GRCh38/gencode.v44.annotation.gtf",
#                 bulk_meta= "/vast/projects/LongBench/sequencing_data/illumina_bulk/metadata.txt",
#                 bulk_rds= "{input.bulk_rds}",
#                 fig.path="{output.figures}/"
#             )
#         )'
#         """
# 
# rule run_de_rmarkdown: 
#     input:
#         rules.rmd_bulk_de_human.output,
#         rules.rmd_sc_sn_de.output
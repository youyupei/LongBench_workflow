results_dir = config["output_path"]

# Analysis with demuxlet
rule rmd_clustering_and_annotation_with_demuxlet:
    input:
        in_10x = os.path.join(results_dir, 'cellranger/{sc_sample}/outs/filtered_feature_bc_matrix'),
        demuxlet_out = os.path.join(results_dir, 'demuxlet/{sc_sample}.demuxlet.best'),
        #fn_empty_gene_count = "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/emptydrop/{ont_sc_sample}_flames/gene_count.csv",
        rmd = os.path.join(config['sub_wf_dir'],  "rmarkdown/ill_sc_with_demuxlet.Rmd"),
        preprocessing_script = os.path.join(config['sub_wf_dir'],  "rmarkdown/sc_processing.R")
    output:
        html = report(os.path.join(results_dir,  "reports/{sc_sample}_clustering_and_annotation.html"),
                        category="SC analysis", subcategory="Clustering and annotation",  labels=lambda wildcards: {
                        "Sample": "{sc_sample}"
                    }),
        out_rds = os.path.join(results_dir,  "reports/RDS/{sc_sample}_annotated.rds")
    params:
        report_title = "SC analysis of {sc_sample}",
        cache_dir = os.path.join(results_dir,  "reports/.cache/clustering_and_annotation_{sc_sample}/"),
        seed = config['random_seed']
    resources:
        mem_mb = 30000,
        cpus_per_task=8
    shell:
        """
        Rscript -e '
            rmarkdown::render("{input.rmd}", 
                                output_file="{output.html}", 
                                params=list(
                                    report_title="{params.report_title}",
                                    random_seed={params.seed},
                                    cache_dir="{params.cache_dir}",
                                    preprocessing_script="{input.preprocessing_script}",
                                    out_rds="{output.out_rds}",
                                    in_10x="{input.in_10x}",
                                    demuxlet_out="{input.demuxlet_out}"
                                )
            )
        '
        """


rule run_all_rmarkdown:
    input:
        expand(os.path.join(results_dir,  "reports/RDS/{sc_sample}_annotated.rds"),
               sc_sample=config['sample_id'])
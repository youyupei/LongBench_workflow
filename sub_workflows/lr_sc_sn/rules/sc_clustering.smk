
results_dir = config["output_path"]

# SC analysis
rule rmd_lr_sc_clustering_and_annotation:
    input:
        fn_flames_gene_quant = os.path.join(results_dir,"flames_out/{sample}/gene_count.csv"),
        demuxlet_out = os.path.join(results_dir,"demuxlet/{sample}.demuxlet.best"),
        fn_empty_gene_count = os.path.join(results_dir, "emptydrop/{sample}_flames/gene_count.csv"),
        rmd = os.path.join(config['sub_wf_dir'],  "rmarkdown/lr_sc_clustering_annotation.Rmd"),
        preprocessing_script = os.path.join(config['sub_wf_dir'],  "rmarkdown/sc_long_preprocessing.R")
    output:
        html = os.path.join(results_dir,  "reports/{sample,.*_sc.*}_clustering_and_annotation.html"),
        annotated_so_rds = os.path.join(results_dir,  "reports/RDS/{sample,.*_sc.*}_annotated.rds")
    params:
        cache_dir = os.path.join(results_dir,  "reports/.cache/clustering_and_annotation_{sample,.*_sc.*}/"),
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
                                    random_seed={params.seed},
                                    cache_dir="{params.cache_dir}",
                                    sc_long_preprocessing_script="{input.preprocessing_script}",
                                    out_annotated_so_rds="{output.annotated_so_rds}",
                                    fn_flames_gene_quant="{input.fn_flames_gene_quant}",
                                    demuxlet_out="{input.demuxlet_out}"
                                )
            )
        '
        """


# SC analysis
rule rmd_lr_sn_clustering_and_annotation:
    input:
        fn_flames_gene_quant = os.path.join(results_dir, "flames_out/{sample}/gene_count.csv"),
        fn_empty_gene_count = os.path.join(results_dir, "emptydrop/{sample}_flames/gene_count.csv"),
        rmd = os.path.join(config['sub_wf_dir'],  "rmarkdown/lr_sn_clustering_annotation.Rmd"),
        preprocessing_script = os.path.join(config['sub_wf_dir'],  "rmarkdown/sc_long_preprocessing.R"),
        ref_sc_rds = os.path.join(results_dir,  "reports/RDS/ont_sc_annotated.rds")
    output:
        html = os.path.join(results_dir,  "reports/{sample,.*_sn.*}_clustering_and_annotation.html"),
        annotated_so_rds = os.path.join(results_dir,  "reports/RDS/{sample,.*_sn.*}_annotated.rds")
    params:
        cache_dir = os.path.join(results_dir,  "reports/.cache/clustering_and_annotation_{sample,.*_sn.*}/"),
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
                                    random_seed={params.seed},
                                    cache_dir="{params.cache_dir}",
                                    sc_long_preprocessing_script="{input.preprocessing_script}",
                                    out_annotated_so_rds="{output.annotated_so_rds}",
                                    fn_flames_gene_quant="{input.fn_flames_gene_quant}",
                                    ref_sc_rds="{input.ref_sc_rds}"
                                )
            )
        '
        """


rule run_all_rmarkdown:
    input:
        expand(os.path.join(results_dir,  "reports/RDS/{sample}_annotated.rds"),
               sample=config["sample_id"])
    output:
        touch(os.path.join(results_dir,  ".flag/clustering_and_annotation_done"))


rule plot_silhouette:
    input:
        rules.run_all_rmarkdown.output
    output:
        report(os.path.join(results_dir,  'plots/silhouette_score.pdf'))
    script:
        "../scripts/plot_silhouette_score.R"
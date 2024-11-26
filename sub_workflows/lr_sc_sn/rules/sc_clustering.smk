
results_dir = config["output_path"]

# SC analysis
rule rmd_clustering_and_annotation:
    input:
        fn_flames_gene_quant = os.path.join(results_dir,"flames_out/{sample}/gene_count.csv"),
        vireo_out = os.path.join(results_dir,"vireo/{sample}_vireo"),
        #fn_empty_gene_count = os.path.join(results_dir, "emptydrop/{sample}_flames/gene_count.csv"),
        rmd = os.path.join(config['sub_wf_dir'],  "rmarkdown/lr_sc_sn_clustering_annotation.Rmd"),
        preprocessing_script = os.path.join(config['sub_wf_dir'],  "rmarkdown/sc_long_preprocessing.R")
    output:
        html = os.path.join(results_dir,  "reports/{sample}_clustering_and_annotation.html"),
        #fig_pdf = os.path.join(results_dir,  "plots/{sample}_clustering_and_annotation.pdf"),
        annotated_so_rds = os.path.join(results_dir,  "reports/RDS/{sample}_annotated.rds")
    params:
        cache_dir = os.path.join(results_dir,  "reports/.cache/clustering_and_annotation_{sample}/"),
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
                                    vireo_out="{input.vireo_out}"
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
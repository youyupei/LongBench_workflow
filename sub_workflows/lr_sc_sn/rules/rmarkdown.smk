
# SC analysis
rule rmd_ont_sc_clustering_and_annotation:
    input:
        fn_flames_gene_quant = "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/flames_out/{ont_sc_sample}/gene_count.csv",
        demuxlet_out = "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/demuxlet/{ont_sc_sample}.demuxlet.best",
        fn_empty_gene_count = "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/emptydrop/{ont_sc_sample}_flames/gene_count.csv",
        rmd = os.path.join(config['sub_wf_dir'],  "rmarkdown/ont_sc_clustering_annotation.Rmd"),
        preprocessing_script = os.path.join(config['sub_wf_dir'],  "rmarkdown/sc_long_preprocessing.R")
    output:
        html = os.path.join(results_dir,  "reports/{ont_sc_sample}_clustering_and_annotation.html"),
        annotated_so_rds = os.path.join(results_dir,  "reports/RDS/sc/{ont_sc_sample}_annotated.rds")
    params:
        cache_dir = os.path.join(results_dir,  "reports/.cache/clustering_and_annotation_{ont_sc_sample}/"),
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
rule rmd_ont_sn_clustering_and_annotation:
    input:
        fn_flames_gene_quant = "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/flames_out/{ont_sn_sample}/gene_count.csv",
        fn_empty_gene_count = "/home/users/allstaff/you.yu/LongBench/analysis/lr_sc_sn/result/emptydrop/{ont_sn_sample}_flames/gene_count.csv",
        rmd = os.path.join(config['sub_wf_dir'],  "rmarkdown/ont_sn_clustering_annotation.Rmd"),
        preprocessing_script = os.path.join(config['sub_wf_dir'],  "rmarkdown/sc_long_preprocessing.R"),
        ref_sc_rds = os.path.join(results_dir,  "reports/RDS/sc/ont_sc_clean_annotated.rds")
    output:
        html = os.path.join(results_dir,  "reports/{ont_sn_sample}_clustering_and_annotation.html"),
        annotated_so_rds = os.path.join(results_dir,  "reports/RDS/sn/{ont_sn_sample}_annotated.rds")
    params:
        cache_dir = os.path.join(results_dir,  "reports/.cache/clustering_and_annotation_{ont_sn_sample}/"),
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
        expand(os.path.join(results_dir,  "reports/RDS/sc/{ont_sc_sample}_annotated.rds"),
               ont_sc_sample=[x for x in config['sample_id'] if "ont_sc" in x]),
        expand(os.path.join(results_dir,  "reports/RDS/sn/{ont_sn_sample}_annotated.rds"),
               ont_sn_sample=[x for x in config['sample_id'] if "ont_sn" in x])

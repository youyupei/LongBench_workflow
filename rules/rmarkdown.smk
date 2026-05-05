# This file manages all the rmarkdown files in the rmarkdown directory.
#
# Note: snakemake input paths are listed for dependency tracking only.
# Each Rmd file uses hard-coded paths internally so it can be run independently
# of the snakemake workflow.
#
# Caching: versioned_knit.R passes cache_dir/fig.path to each Rmd.
# knitr chunk caches are reused when the Rmd `version:` field is unchanged,
# and invalidated (fresh run) when the version is bumped.
#
# ─── Dependency graph ────────────────────────────────────────────────────────
#
#  [GTF]                    [quantification outputs / raw data dirs]
#    │                                      │
#    ▼                                      │
#  rmd_intronic_gene_and_exon_count         │          (Stage 1 – independent)
#    │  intronic_gene_and_exon_count.rds    │         ┌─ rmd_QC_plots
#    │                                      │         └─ rmd_mutation_analysis
#  r_tx2gene_map                            │
#    │  Tx2Gene.map.rds                     │
#    ├──────────────────────────────────────┤
#    │                                      │
#    ▼                                      ▼
#  r_get_bulk_DGE_objects        rmd_bulk_identification_analysis_20M
#    │  bulk_DGE.obj.rds         rmd_bulk_GC_analysis_20M
#    │
#    ├─ rmd_bulk_identification_analysis ──► bulk_identification.rds ──┐
#    │                                                                  │
#    ├─ rmd_bulk_DE_analysis ─────────────► bulk_DE.rds ───────────────┼──► rmd_sc_pseudo_bulk_analysis
#    │                                          │                       │
#    ├─ rmd_bulk_quantification_analysis        ├─ rmd_bulk_DE_analysis_continue
#    ├─ rmd_bulk_quantification_analysis_isoquant                       │
#    ├─ rmd_bulk_quantification_cross_tool_comparison                   │
#    ├─ rmd_bulk_DE_spikeins                    └─ rmd_bulk_DE_analysis_subsample_20M
#    └─ rmd_bulk_DE_spikeins
#
#  r_get_sc_DGE_objects ──► sc_DGE.obj.rds ──── ──────────────────────► rmd_sc_pseudo_bulk_analysis
#
#  rmd_sc_sn_umap ──► sc_sn_filtered_so.rds
#  rmd_sc_clustering_annotation   (reads FLAMES / Vireo / CellRanger directly)
#  rmd_sn_clustering_annotation   (reads FLAMES / Vireo / CellRanger directly)
#
# ─────────────────────────────────────────────────────────────────────────────

rmd_dir = join(config['main_wf_dir'], 'rmarkdown')
rmd_output_dir = join(config['main_wf_dir'], 'rmarkdown')
flag_dir = config['flag_dir']
knit_script = join(rmd_dir, 'versioned_knit.R')
Rscript_dir = join(rmd_dir, 'Rscript')

lr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'], config["sub_wf_config"]["lr_sc_sn"]))
main_wf_config = config_parser.load_configfile(join(config['main_wf_dir'], config["sub_wf_config"]["main_wf"]))


# ─── Stage 0: Preprocessing R scripts ────────────────────────────────────────

rule r_tx2gene_map:
    output:
        rds = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/Tx2Gene.map.rds'
    resources:
        cpus_per_task = 4,
        mem_mb = 32000
    script:
        join(Rscript_dir, 'Tx2Gene.map.R')


rule r_get_bulk_DGE_objects:
    input:
        rules.r_tx2gene_map.output,
    output:
        bulk_DGE_object = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DGE.obj.rds'
    resources:
        cpus_per_task = 1,
        mem_mb = 8000
    script:
        join(Rscript_dir, 'Bulk.DGElist.preprocessing.R')


rule r_get_bulk_DGE_objects_IP_filtered:
    input:
        rules.r_tx2gene_map.output,
    output:
        bulk_DGE_object = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DGE.IP_filtered.obj.rds'
    resources:
        cpus_per_task = 1,
        mem_mb = 8000
    script:
        join(Rscript_dir, 'Bulk.DGElist.preprocessing.IP_filtered.R')


rule r_get_sc_DGE_objects:
    input:
        rules.r_tx2gene_map.output,
    output:
        sc_DGE_object = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/sc_DGE.obj.rds'
    resources:
        cpus_per_task = 1,
        mem_mb = 8000
    script:
        join(Rscript_dir, 'Sc.DGElist.preprocessing.R')


rule rmd_intronic_gene_and_exon_count:
    input:
        gtf = '/vast/projects/LongBench/reference_files/GRCh38/gencode.v44.annotation.gtf'
    output:
        rds = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/intronic_gene_and_exon_count.rds'
    resources:
        cpus_per_task = 16,
        mem_mb = 32000
    script:
        join(rmd_dir, 'get_intronic_gene.R')


# ─── Stage 1: Independent analyses ───────────────────────────────────────────

rule rmd_QC_plots:
    input:
        rmd = join(rmd_dir, 'QC_plot.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.rmd_QC_plots.done'))
    resources:
        cpus_per_task = 1,
        mem_mb = 8000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


rule rmd_mutation_analysis:
    input:
        rmd = join(rmd_dir, 'Mutation_analysis.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.mutation_analysis.done'))
    resources:
        cpus_per_task = 1,
        mem_mb = 8000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "module load pandoc && Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


# ─── Stage 2a: Bulk identification (→ bulk_identification.rds) ───────────────

rule rmd_bulk_identification_analysis:
    input:
        rmd = join(rmd_dir, 'Bulk_identification.Rmd'),
        knit_script = knit_script,
        bulk_dge = rules.r_get_bulk_DGE_objects.output.bulk_DGE_object,
        intronic_rds = rules.rmd_intronic_gene_and_exon_count.output.rds
    output:
        rds = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_identification.rds'
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


rule rmd_bulk_identification_analysis_IP_filtered:
    input:
        rmd = join(rmd_dir, 'Bulk_identification_IP_filtered.Rmd'),
        knit_script = knit_script,
        bulk_dge = rules.r_get_bulk_DGE_objects_IP_filtered.output.bulk_DGE_object,
        intronic_rds = rules.rmd_intronic_gene_and_exon_count.output.rds
    output:
        rds = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_identification_IP_filtered.rds'
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    retries: 3
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        module load pandoc
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """


# ─── Stage 2b: Bulk DE summary (→ bulk_DE.rds) ───────────────────────────────

rule rmd_bulk_DE_analysis:
    input:
        rmd = join(rmd_dir, 'Bulk_DE_Summary.Rmd'),
        knit_script = knit_script,
        bulk_dge = rules.r_get_bulk_DGE_objects.output.bulk_DGE_object
    output:
        rds = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DE.rds'
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    retries: 3
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        module load curl
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """


# ─── Stage 2c: Quantification comparisons ────────────────────────────────────

rule rmd_bulk_quantification_analysis:
    input:
        rmd = join(rmd_dir, 'Bulk_quantification_analysis.Rmd'),
        knit_script = knit_script,
        bulk_dge = rules.r_get_bulk_DGE_objects.output.bulk_DGE_object,
        tx2gene = rules.r_tx2gene_map.output.rds
    output:
        touch(join(flag_dir, 'rmd.Bulk_quantification_analysis.done'))
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


rule rmd_bulk_quantification_analysis_isoquant:
    input:
        rmd = join(rmd_dir, 'Bulk_quantification_analysis_IsoQuant.Rmd'),
        knit_script = knit_script,
        bulk_dge = rules.r_get_bulk_DGE_objects.output.bulk_DGE_object,
        tx2gene = rules.r_tx2gene_map.output.rds
    output:
        touch(join(flag_dir, 'rmd.Bulk_quantification_analysis_IsoQuant.done'))
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        module load pandoc
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """


rule rmd_bulk_quantification_cross_tool_comparison:
    input:
        rmd = join(rmd_dir, 'Bulk_quantification_cross_tool_comparison.Rmd'),
        knit_script = knit_script,
        bulk_dge = rules.r_get_bulk_DGE_objects.output.bulk_DGE_object
    output:
        touch(join(flag_dir, 'rmd.Bulk_quantification_cross_tool_comparison.done'))
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


# ─── Stage 2d: Spike-in DE ────────────────────────────────────────────────────

rule rmd_bulk_DE_spikeins:
    input:
        rmd = join(rmd_dir, 'Bulk_DE.spikeins.Rmd'),
        knit_script = knit_script,
        bulk_dge = rules.r_get_bulk_DGE_objects.output.bulk_DGE_object,
        tx2gene = rules.r_tx2gene_map.output.rds
    output:
        touch(join(flag_dir, 'rmd.bulk_DE_spikeins.done'))
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    retries: 3
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        module load curl
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """


# ─── Stage 2e: 20M subsample analyses ────────────────────────────────────────

rule rmd_bulk_identification_analysis_20M:
    input:
        rmd = join(rmd_dir, 'Bulk_identification_20M.Rmd'),
        knit_script = knit_script,
        tx2gene = rules.r_tx2gene_map.output.rds
    output:
        touch(join(flag_dir, 'rmd.bulk_identification_20M.done'))
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


rule rmd_bulk_GC_analysis_20M:
    input:
        rmd = join(rmd_dir, 'Bulk_GC_analysis_20M.Rmd'),
        knit_script = knit_script,
        tx2gene = rules.r_tx2gene_map.output.rds,
        intronic_rds = rules.rmd_intronic_gene_and_exon_count.output.rds
    output:
        touch(join(flag_dir, 'rmd.Bulk_GC_analysis_20M.done'))
    resources:
        cpus_per_task = 4,
        mem_mb = 16000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


# ─── Stage 3: Downstream bulk (require bulk_DE.rds) ──────────────────────────

rule rmd_bulk_DE_analysis_continue:
    input:
        rmd = join(rmd_dir, 'Bulk_DE_analysis_continue.Rmd'),
        knit_script = knit_script,
        bulk_de = rules.rmd_bulk_DE_analysis.output.rds
    output:
        touch(join(flag_dir, 'rmd.bulk_DE_analysis_continue.done'))
    resources:
        cpus_per_task = 4,
        mem_mb = 16000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        module load pandoc
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """


rule rmd_bulk_DE_analysis_subsample_20M:
    input:
        rmd = join(rmd_dir, 'Bulk_DE_Summary_20M.Rmd'),
        knit_script = knit_script,
        bulk_de = rules.rmd_bulk_DE_analysis.output.rds
    output:
        touch(join(flag_dir, 'rmd.rmd_bulk_DE_analysis_20M.done'))
    resources:
        cpus_per_task = 8,
        mem_mb = 32000
    retries: 3
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        module load curl
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """


# ─── Stage 4: Single-cell preprocessing (→ sc_sn_filtered_so.rds) ────────────

rule rmd_sc_sn_umap:
    input:
        rmd = join(rmd_dir, 'sc_sn_umap.Rmd'),
        knit_script = knit_script
    output:
        rds = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/sc_sn_filtered_so.rds'
    resources:
        cpus_per_task = 2,
        mem_mb = 64000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


rule rmd_sc_clustering_annotation:
    input:
        rmd = join(rmd_dir, 'sc_clustering_annotation.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.sc_clustering_annotation.done'))
    resources:
        cpus_per_task = 4,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


rule rmd_sn_clustering_annotation:
    input:
        rmd = join(rmd_dir, 'sn_clustering_annotation.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.sn_clustering_annotation.done'))
    resources:
        cpus_per_task = 4,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


# ─── Stage 5: Single-cell pseudobulk ─────────────────────────────────────────

rule rmd_sc_pseudo_bulk_analysis:
    input:
        rmd = join(rmd_dir, 'SC_identification_DE_analysis.Rmd'),
        knit_script = knit_script,
        pre_steps = [
            rules.lr_sc_sn_pseudo_bulk_oarfish_map_n_quant.output,
            '/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkQC/pseudo_bulk_read_count.csv',
        ],
        bulk_de = rules.rmd_bulk_DE_analysis.output.rds,
        bulk_ident = rules.rmd_bulk_identification_analysis.output.rds,
        sc_dge = rules.r_get_sc_DGE_objects.output.sc_DGE_object
    output:
        touch(join(flag_dir, 'rmd.sc_pseudo_bulk_analysis.done'))
    resources:
        cpus_per_task = 1,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        "Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}"


rule _rmd_sc_pseudo_bulk_analysis_subsample_single_run:
    input:
        rmd = join(rmd_dir, 'SC_identification_DE_analysis.Rmd'),
        knit_script = knit_script,
        pre_steps = [
            '/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkQC/pseudo_bulk_read_count.csv',
        ],
        bulk_de = rules.rmd_bulk_DE_analysis.output.rds,
        bulk_ident = rules.rmd_bulk_identification_analysis.output.rds,
        out_dir_cov = [
            os.path.join(main_wf_config['output_path'], "rarefraction_analysis/oarfish/ont_sc/{subsample_size,.*M}"),
            os.path.join(main_wf_config['output_path'], "rarefraction_analysis/oarfish/ont_sn/{subsample_size,.*M}"),
            os.path.join(main_wf_config['output_path'], "rarefraction_analysis/oarfish/pb_sc/{subsample_size,.*M}"),
            os.path.join(main_wf_config['output_path'], "rarefraction_analysis/oarfish/pb_sn/{subsample_size,.*M}")
        ]
    output:
        touch(join(flag_dir, 'rmd.sc_pseudo_bulk_analysis_subsample_{subsample_size,.*M}.done'))
    resources:
        cpus_per_task = 1,
        mem_mb = 32000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        tmp_rmd=$(dirname {input.rmd})/$(basename {input.rmd} .Rmd)_subsample_{wildcards.subsample_size}.rmd
        cp {input.rmd} $tmp_rmd
        sed -i "s|/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sc|{input.out_dir_cov[0]}|g" $tmp_rmd
        sed -i "s|/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sn|{input.out_dir_cov[1]}|g" $tmp_rmd
        sed -i "s|/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sc|{input.out_dir_cov[2]}|g" $tmp_rmd
        sed -i "s|/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sn|{input.out_dir_cov[3]}|g" $tmp_rmd
        Rscript {input.knit_script} $tmp_rmd {params.rmd_output_dir}
        rm $tmp_rmd
        """


rule rmd_sc_pseudo_bulk_analysis_subsample:
    input:
        expand(rules._rmd_sc_pseudo_bulk_analysis_subsample_single_run.output, subsample_size=['30M'])
    output:
        touch(join(flag_dir, 'rmd.sc_pseudo_bulk_analysis_subsample.done'))
    localrule: True
    shell:
        "touch {output}"


# ─── Aggregator ───────────────────────────────────────────────────────────────

rule rmd_dRNA_run_comparison:
    input:
        rmd = join(rmd_dir, 'dRNA_run_comparison.Rmd'),
        knit_script = knit_script,
        oarfish = expand(
            "/vast/projects/LongBench/analysis/dRNA_comparison/result/oarfish_cov_output/{sample}/{cell_line}",
            sample=["dRNA_bulk_run1", "dRNA_bulk_run2"],
            cell_line=sub_wf_config['lr_bulk']['cell_lines']),
        read_counts = expand(
            "/vast/projects/LongBench/analysis/dRNA_comparison/result/qc/read_counts/{sample}_{cell_line}.count",
            sample=["dRNA_bulk_run1", "dRNA_bulk_run2"],
            cell_line=sub_wf_config['lr_bulk']['cell_lines'])
    output:
        touch(join(flag_dir, 'rmd.dRNA_run_comparison.done'))
    resources:
        cpus_per_task=4,
        mem_mb=16000
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        export RSTUDIO_PANDOC=/stornext/System/data/software/rhel/9/base/tools/pandoc/3.2/bin
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """


rule knit_rmarkdown:
    input:
        # Stage 1
        rules.rmd_QC_plots.output,
        rules.rmd_mutation_analysis.output,
        # Stage 2
        rules.rmd_bulk_identification_analysis.output,
        rules.rmd_bulk_DE_analysis.output,
        rules.rmd_bulk_quantification_analysis.output,
        rules.rmd_bulk_quantification_analysis_isoquant.output,
        rules.rmd_bulk_quantification_cross_tool_comparison.output,
        rules.rmd_bulk_DE_spikeins.output,
        rules.rmd_bulk_identification_analysis_20M.output,
        rules.rmd_bulk_GC_analysis_20M.output,
        # Stage 3
        rules.rmd_bulk_DE_analysis_continue.output,
        rules.rmd_bulk_DE_analysis_subsample_20M.output,
        # Stage 4
        rules.rmd_sc_sn_umap.output,
        rules.rmd_sc_clustering_annotation.output,
        rules.rmd_sn_clustering_annotation.output,
        # Stage 5
        rules.rmd_sc_pseudo_bulk_analysis.output,
        rules.rmd_sc_pseudo_bulk_analysis_subsample.output,
        # dRNA comparison
        rules.rmd_dRNA_run_comparison.output
    default_target: True
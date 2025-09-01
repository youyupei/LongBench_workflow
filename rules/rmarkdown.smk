# This file manages all the rmarkdown files in the rmarkdown directory
rmd_dir = join(config['main_wf_dir'], 'rmarkdown')
rmd_output_dir = join(config['main_wf_dir'], 'rmarkdown')
flag_dir = config['flag_dir']
knit_script = join(rmd_dir, 'versioned_knit.R')
Rscript_dir = join(rmd_dir, 'Rscript')


lr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_sc_sn"]))
main_wf_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["main_wf"]))

# Part 1 preprocessing R scripts
## 1. Get the transcript to gene mapping for human, sequins and SIRV
rule r_tx2gene_map:
    output:
        rds = '/home/users/allstaff/you.yu/LongBench/analysis/workflow/rmarkdown/RDS/Tx2Gene.map.rds'
    resources:
        cpus_per_task=4,
        mem_mb=32000
    script:
        join(Rscript_dir, 'Tx2Gene.map.R')


## 1. Unfiltered bulk DGE in both Gene and Transcript level
## 2. Unfiltered scRNA-seq pseudobulk DGE in both Gene and Transcript level
rule r_get_bulk_DGE_objects:
    input:
        rules.r_tx2gene_map.output,
    output:
        bulk_DGE_object = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DGE.obj.rds'
    resources:
        cpus_per_task=1,
        mem_mb=8000
    script:
        join(Rscript_dir, 'Bulk.DGElist.preprocessing.R')

rule r_get_sc_DGE_objects:
    input:
        rules.r_tx2gene_map.output,
    output:
        sc_DGE_object = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/sc_DGE.obj.rds'
    resources:
        cpus_per_task=1,
        mem_mb=8000
    script:
        join(Rscript_dir, 'Sc.DGElist.preprocessing.R')

## 3. Getting intronic gene and exon count
rule rmd_intronic_gene_and_exon_count:
    input:
        gtf = '/vast/projects/LongBench/reference_files/GRCh38/gencode.v44.annotation.gtf'
    output:
        rds = '/home/users/allstaff/you.yu/LongBench/analysis/workflow/rmarkdown/RDS/intronic_gene_and_exon_count.rds'
    resources:
        cpus_per_task=16,
        mem_mb=32000
    script:
        join(rmd_dir, 'get_intronic_gene.R')


# QC analysis
## BUlK
rule rmd_QC_plots: 
    # The input path is hard-coded here
    input:
        rmd = join(rmd_dir, 'QC_plot.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.rmd_QC_plots.done'))
    resources:
        cpus_per_task=1
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """


# Note: In the smk file, the snakemake input path is not directly used in the Rmd file.
# Instead, the Rmd file uses the hard-coded path to the required file.
# This is to make each Rmd file independent of the snakemake workflow so can be run independently.
# BUlK analysis
rule rmd_BCV_plot:
    input:
        rmd = join(rmd_dir, 'Bulk_BCV.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.BCV_plot.done'))
    resources:
        cpus_per_task=1
    params:
        rmd_output_dir = rmd_output_dir
    shell:
        """
        Rscript {input.knit_script} {input.rmd}  {params.rmd_output_dir}
        """


rule rmd_bulk_identification_analysis:
    input:
        rmd = join(rmd_dir, 'Bulk_identification.Rmd'),
        knit_script = knit_script,
        rds = '/home/users/allstaff/you.yu/LongBench/analysis/workflow/rmarkdown/RDS/intronic_gene_and_exon_count.rds'
    output:
        touch(join(flag_dir, 'rmd.rmd_bulk_identification_analysis.done'))
        # generates "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_identification.rds"
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=8,
        mem_mb=32000
    shell:
        """
        Rscript {input.knit_script} {input.rmd}  {params.rmd_output_dir}
        """


rule rmd_bulk_quantification_analysis:
    input:
        rmd = join(rmd_dir, 'Bulk_quantification_analysis.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.Bulk_quantification_analysis.done'))
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=8,
        mem_mb=32000
    shell:
        """
        Rscript {input.knit_script} {input.rmd}  {params.rmd_output_dir}
        """

rule rmd_bulk_DE_analysis:
    input:
        rds = rules.r_load_DGE_objects.output.bulk_DGE_object,
        rmd = join(rmd_dir, 'Bulk_DE_Summary.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.rmd_bulk_DE_analysis.done'))
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=8,
        mem_mb=32000
    retries: 3
    shell:
        """
        module load curl
        Rscript {input.knit_script} {input.rmd}  {params.rmd_output_dir}
        """

rule rmd_bulk_DE_analysis_subsample_20M:
    input:
        rmd = join(rmd_dir, 'Bulk_DE_Summary_20M.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.rmd_bulk_DE_analysis_20M.done'))
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=8,
        mem_mb=32000
    retries: 3
    shell:
        """
        module load curl
        Rscript {input.knit_script} {input.rmd}  {params.rmd_output_dir}
        """

rule rmd_bulk_DTU_anlaysis:
    input:
        rmd = join(rmd_dir, 'Bulk_DTU.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.bulk_DTU_anlaysis.done'))
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=16,
        mem_mb=128000
    shell:
        """
        module load gcc
        Rscript {input.knit_script} {input.rmd}  {params.rmd_output_dir}
        """

# SC
## PseudoBulk
rule rmd_sc_pseudo_bulk_analysis:
    input:
        rmd = join(rmd_dir, 'SC_identification_DE_analysis.Rmd'),
        knit_script = knit_script,
        pre_steps = [
            rules.lr_sc_sn_pseudo_bulk_oarfish_map_n_quant.output,
            '/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkQC/pseudo_bulk_read_count.csv',
            # Needs '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DE.rds', 
            rules.rmd_bulk_DE_analysis.output,
            # needs "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_identification.rds" # Bulk identification result
            rules.rmd_bulk_identification_analysis.output
        ]
    output:
        touch(join(flag_dir, 'rmd.sc_pseudo_bulk_analysis.done'))
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=1,
        mem_mb=32000
    shell:
        """
        Rscript {input.knit_script} {input.rmd}  {params.rmd_output_dir}
        """

rule _rmd_sc_pseudo_bulk_analysis_subsample_single_run:
    input:
        rmd = join(rmd_dir, 'SC_identification_DE_analysis.Rmd'),
        knit_script = knit_script,
        pre_steps = [
            '/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkQC/pseudo_bulk_read_count.csv',
            # Needs '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DE.rds', 
            rules.rmd_bulk_DE_analysis.output,
            # needs "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_identification.rds" # Bulk identification result
            rules.rmd_bulk_identification_analysis.output
        ],
        out_dir_cov = [
            os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/ont_sc/{subsample_size,.*M}"),
            os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/ont_sn/{subsample_size,.*M}"),
            os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/pb_sc/{subsample_size,.*M}"),
            os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/pb_sn/{subsample_size,.*M}")
        ]
    output:
        touch(join(flag_dir, 'rmd.sc_pseudo_bulk_analysis_subsample_{subsample_size,.*M}.done'))
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=1,
        mem_mb=32000
    shell:
        """
        tmp_rmd=$(dirname {input.rmd})/$(basename {input.rmd} .Rmd)_subsample_{wildcards.subsample_size}.rmd
        cp {input.rmd} $tmp_rmd
        # replace the xxxx with yyyyy from the tmp_rmd in place
        sed -i "s|/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sc|{input.out_dir_cov[0]}|g" $tmp_rmd
        sed -i "s|/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/ont_sn|{input.out_dir_cov[1]}|g" $tmp_rmd
        sed -i "s|/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sc|{input.out_dir_cov[2]}|g" $tmp_rmd
        sed -i "s|/vast/projects/LongBench/analysis/lr_sc_sn/result/PseudoBulkOarfishCov/pb_sn|{input.out_dir_cov[3]}|g" $tmp_rmd
        Rscript {input.knit_script} $tmp_rmd  {params.rmd_output_dir}
        rm $tmp_rmd
        """

rule rmd_sc_pseudo_bulk_analysis_subsample:
    input:
        expand(rules._rmd_sc_pseudo_bulk_analysis_subsample_single_run.output, subsample_size = ['30M'])
    output:
        touch(join(flag_dir, 'rmd.sc_pseudo_bulk_analysis_subsample.done'))
    resources:
        cpus_per_task=1,
        mem_mb=32000
    localrule: True
    shell:
        """
        touch {output}
        
        """

# rule rmd_sc_preprocessing:
# 
# rule rmd_sn_preprocessing:

rule rmd_sc_sn_plot:
    input:
        rmd = join(rmd_dir, 'sc_sn_plot.Rmd'),
        knit_script = knit_script,
    output:
        rds = '/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/sc_sn_filtered_so.rds'
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=2,
        mem_mb=64000
    shell:
        """
        Rscript {input.knit_script} {input.rmd} {params.rmd_output_dir}
        """

# 
# rule rmd_sc_sn_marker_analysis:
#     input:
#         rds = rules.rmd_sc_sn_plot.output.rds
# 
# rule rmd_sc_sn_pseudo_bulk_analysis:


# rare faction analysis
rule rmd_rarefaction_analysis:
    input:
        rmd = join(rmd_dir, 'rarefaction_analysis.Rmd'),
        knit_script = knit_script,
        pre_steps = [
            rules.rarefraction_analysis.output
        ]
    output:
        touch(join(flag_dir, 'rmd.rarefaction_analysis.done'))
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=1,
        mem_mb=32000
    shell:
        """
        Rscript {input.knit_script} {input.rmd}  {params.rmd_output_dir}
        """ 


rule knit_rmarkdown:
    input:
        rules.rmd_BCV_plot.output,
        rules.rmd_QC_plots.output,
        rules.rmd_bulk_DE_analysis.output,
        rules.rmd_bulk_DE_analysis_subsample_20M.output,
        rules.rmd_bulk_DTU_anlaysis.output,
        rules.rmd_sc_pseudo_bulk_analysis.output,
        rules.rmd_bulk_quantification_analysis.output,
        rules.rmd_sc_pseudo_bulk_analysis_subsample.output,
        rules.rmd_sc_sn_plot.output
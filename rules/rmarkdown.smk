# This file manages all the rmarkdown files in the rmarkdown directory
rmd_dir = join(config['main_wf_dir'], 'rmarkdown')
rmd_output_dir = join(config['main_wf_dir'], 'rmarkdown')
flag_dir = config['flag_dir']
knit_script = join(rmd_dir, 'versioned_knit.R')

lr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_sc_sn"]))
main_wf_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["main_wf"]))

# Note: In the smk file, the snakemake input path is not directly used in the Rmd file.
# Instead, the Rmd file uses the hard-coded path to the required file.
# This is to make each Rmd file independent of the snakemake workflow so can be run independently.

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
        knit_script = knit_script
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
        rmd = join(rmd_dir, 'Bulk_DE_Summary.Rmd'),
        knit_script = knit_script
    output:
        touch(join(flag_dir, 'rmd.rmd_bulk_DE_analysis.done'))
        # generates "/vast/projects/LongBench/analysis/workflow/rmarkdown/RDS/bulk_DE.rds"
    params:
        rmd_output_dir = rmd_output_dir
    resources:
        cpus_per_task=8,
        mem_mb=32000
    shell:
        """
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
        cpus_per_task=10,
        mem_mb=32000
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
            rules.lr_sc_sn_pseudo_bulk_oarfish_quant.output,
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

rule knit_rmarkdown:
    input:
        rules.rmd_BCV_plot.output,
        rules.rmd_QC_plots.output,
        rules.rmd_bulk_DE_analysis.output,
        rules.rmd_bulk_DTU_anlaysis.output,
        rules.rmd_sc_pseudo_bulk_analysis.output,
        rules.rmd_bulk_quantification_analysis.output,
        rules.rmd_sc_pseudo_bulk_analysis_subsample.output,
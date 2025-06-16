from os.path import join
# lr_bulk_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_bulk"]))
# sr_bulk_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_bulk"]))
lr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_sc_sn"]))
# sr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_sc_sn"]))
main_wf_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["main_wf"]))
# main_conda = join(config['main_wf_dir'], config['conda_config']['main'])

rule lr_sc_bam_downsample:
    input:
        bam=lr_sc_sn_config['output_path'] + "/PseudoBulkAlignment/Transcriptome/{sample}_{cell_line}.bam",
        csv=lr_sc_sn_config['output_path'] + "/PseudoBulkQC/pseudo_bulk_read_count.csv" 
    output:
        join(main_wf_config['tmp_dir'] , "PseudoBulkAlignment/{sample}_{cell_line}.read_id_subsampled.{subsample_rate}.bam")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        tmp_dir = main_wf_config['tmp_dir'],
        seed = config['random_seed']
    shell:
        """
        mkdir -p $(dirname {output})
        samtools view {input.bam} | cut -f1 | sort | uniq > {output}.read_ids
        awk 'BEGIN {{srand({params.seed})}} {{if (rand() < {wildcards.subsample_rate}) print $0}}' {output}.read_ids > {output}.read_ids.subsampled
        rm {output}.read_ids
        samtools view -b -N {output}.read_ids.subsampled {input.bam} > {output}
        rm {output}.read_ids.subsampled
        """


rule lr_sc_oarfish_cov_rare_fraction_analysis:
    input:
        bam = rules.lr_sc_bam_downsample.output,
        ref = lr_sc_sn_config["reference"]["transcript"]
    output:
        out_dir_cov = directory(os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/sc_sn/{sample}/{subsample_rate,0.*}/{cell_line}"))
    conda:
        join(config['main_wf_dir'], config["conda_config"]["oarfish"])
    resources:
        cpus_per_task=16,
        mem_mb=64000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A"
    priority: 100
    shell:
        """
        mkdir -p {output.out_dir_cov}
        oarfish --alignments {input.bam} --threads {resources.cpus_per_task} --output {output.out_dir_cov}/ --model-coverage  -d . --filter-group no-filters --num-bootstraps 50
        """

rule lr_sc_link_oarfish_cov_full:
    input:
        rules.lr_sc_sn_run_oarfish_cov.output
    output:
        directory(os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/sc_sn/{sample}/full/{cell_line}"))
    localrule: True
    shell:
        """
        mkdir -p $(dirname {output})
        ln -s {input} {output}
        sleep 1
        touch -h {output}
        """

# Entire quantification worflow head
rule lr_sc_rarefraction_analysis:
    input:
        expand(
            [
                rules.lr_sc_oarfish_cov_rare_fraction_analysis.output[0],
                rules.lr_sc_bam_downsample.output[0],
                rules.lr_sc_link_oarfish_cov_full.output[0]
            ],
            cell_line = lr_sc_sn_config['cell_line_list'],
            sample = [x for x in lr_sc_sn_config['sample_id'] if "sc" in x],
            subsample_rate =  main_wf_config['rarefaction_rate_sc']
        ),
        expand(
            [
                rules.lr_sc_oarfish_cov_rare_fraction_analysis.output[0],
                rules.lr_sc_bam_downsample.output[0],
                rules.lr_sc_link_oarfish_cov_full.output[0]
            ],
            cell_line = lr_sc_sn_config['cell_line_list'],
            sample = [x for x in lr_sc_sn_config['sample_id'] if "sn" in x],
            subsample_rate =  main_wf_config['rarefaction_rate_sn']
        )
    output:
        touch(join(config['flag_dir'], "lr_sc_rarefraction_analysis.done"))
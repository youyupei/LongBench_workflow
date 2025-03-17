from os.path import join
lr_bulk_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_bulk"]))
sr_bulk_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_bulk"]))
lr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_sc_sn"]))
sr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_sc_sn"]))
main_wf_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["main_wf"]))
# main_conda = join(config['main_wf_dir'], config['conda_config']['main'])

rule lr_bam_downsample:
    input:
        lr_bulk_config['output_path'] + "/TranscriptAlignment/{sample}_{cell_line}.bam"
    output:
        join(main_wf_config['tmp_dir'] , "TranscriptAlignment/{sample}_{cell_line}.read_id_subsampled.{subsample_size}.bam")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        tmp_dir = main_wf_config['tmp_dir'],
        seed = config['random_seed']
    shell:
        """
        mkdir -p $(dirname {output})
        samtools view {input} | cut -f1 | sort | uniq > {output}.read_ids
        shuf -n $(numfmt --from=si  {wildcards.subsample_size}) --random-source=<(yes {params.seed}) {output}.read_ids  > {output}.read_ids.subsampled
        rm {output}.read_ids
        samtools view -b -N {output}.read_ids.subsampled {input} > {output}
        """


rule oarfish_cov_rare_fraction_analysis:
    input:
        bam = rules.lr_bam_downsample.output,
        ref = lr_bulk_config["reference"]["transcript"]
    output:
        out_dir_cov = directory(os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/{sample}/{subsample_size,.*M}/{cell_line}"))
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

rule link_lr_bulk_oarfish_cov_full:
    input:
        rules.lr_bulk_oarfish_cov.output # directory(os.path.join(results_dir,"oarfish_cov_output/{sample}/{cell_line}"))
    output:
        directory(os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/{sample}/full/{cell_line}"))
    localrule: True
    shell:
        """
        sleep 5
        ln -s {input} {output}
        """

# SR bulk
## Step 2.1: Index the transcriptome
rule salmon_quant_downsample:
    input:
        R1 = rules.sr_bulk_fastp.output.R1,
        R2 = rules.sr_bulk_fastp.output.R2,
        index = rules.sr_bulk_salmon_index.output
    output:
         directory(join(main_wf_config['output_path'], "rarefraction_analysis/salmon/{subsample_size,.*M}/{cell_line}"))
    conda:
        join(config['main_wf_dir'], config["conda_config"]["main"])
    resources:
        cpus_per_task=8,
        mem_mb=32000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    priority: 1
    shell:
        """
        mkdir -p {output}
        # subsample the fastq files
        seqtk sample -s100 {input.R1} $(numfmt --from=si  {wildcards.subsample_size})> {input.R1}.subsampled.{wildcards.subsample_size}.fastq
        seqtk sample -s100 {input.R2} $(numfmt --from=si  {wildcards.subsample_size}) > {input.R2}.subsampled.{wildcards.subsample_size}.fastq

        salmon quant -i {input.index} \
                    -l A \
                    -1 {input.R1}.subsampled.{wildcards.subsample_size}.fastq \
                    -2 {input.R2}.subsampled.{wildcards.subsample_size}.fastq \
                    --validateMappings \
                    -o {output} \
                    -p 8 \
                    --numBootstraps 50
        
        rm {input.R1}.subsampled.{wildcards.subsample_size}.fastq
        rm {input.R2}.subsampled.{wildcards.subsample_size}.fastq
        """

rule link_sr_bulk_salmon_full:
    input:
        rules.sr_bulk_salmon_quant.output # directory(join(results_dir, "salmon/salmon_quant/{cell_line}"))
    output:
        directory(join(main_wf_config['output_path'], "rarefraction_analysis/salmon/1.0/{cell_line}"))
    localrule: True
    shell:
        """
        sleep 5
        ln -s {input} {output}
        """


# Entire quantification worflow head
rule rarefraction_analysis:
    input:
        expand(
            [
                rules.salmon_quant_downsample.output[0],
                rules.oarfish_cov_rare_fraction_analysis.output[0],
                rules.lr_bam_downsample.output[0],
                rules.link_lr_bulk_oarfish_cov_full.output[0],
                rules.link_sr_bulk_salmon_full.output[0]
            ],
            cell_line = sr_bulk_config['cell_lines'],
            sample = lr_bulk_config['sample_id'],
            subsample_size =  main_wf_config['rarefaction_size']
        )
    output:
        touch(join(config['flag_dir'], "rarefraction_analysis.done"))
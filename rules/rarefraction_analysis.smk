from os.path import join
lr_bulk_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_bulk"]))
sr_bulk_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_bulk"]))
lr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_sc_sn"]))
sr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_sc_sn"]))
main_wf_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["main_wf"]))
# main_conda = join(config['main_wf_dir'], config['conda_config']['main'])

rule lr_bam_downsample:
    input:
        "{x}.bam"
    output:
        temp("{x}.read_id_subsampled.{subsample_rate}.bam")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    shadow: "full"
    shell:
        """
        samtools view {input} | cut -f1 > {input}.read_ids
        awk 'BEGIN {{ {{srand()}} {{if (rand() < {wildcards.subsample_rate}) print $0}} }}' {input}.read_ids > {input}.read_ids.subsampled
        rm {input}.read_ids
        samtools view -b -N {input}.read_ids.subsampled {input} > {output}
        """


rule oarfish_cov_rare_fraction_analysis:
    input:
        bam = lr_bulk_config['output_path'] + "/TranscriptAlignment/{sample}_{cell_line}.read_id_subsampled.{subsample_rate}.bam",
        ref = lr_bulk_config["reference"]["transcript"]
    output:
        out_dir_cov = directory(os.path.join(main_wf_config['output_path'],"rarefraction_analysis/oarfish/{sample}/{subsample_rate}_{cell_line}"))
    conda:
        join(config['main_wf_dir'], config["conda_config"]["oarfish"])
    resources:
        cpus_per_task=16,
        mem_mb=64000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A"
    priority: 101
    shell:
        """
        mkdir -p {output.out_dir_cov}
        oarfish --alignments {input.bam} --threads {resources.cpus_per_task} --output {output.out_dir_cov}/ --model-coverage  -d . --filter-group no-filters
        """

# SR bulk
## Step 2.1: Index the transcriptome
rule salmon_quant_downsample:
    input:
        R1 = rules.sr_bulk_fastp.output.R1,
        R2 = rules.sr_bulk_fastp.output.R2,
        index = rules.sr_bulk_salmon_index.output
    output:
         directory(join(main_wf_config['output_path'], "rarefraction_analysis/salmon/{subsample_rate}_{cell_line}"))
    conda:
        join(config['main_wf_dir'], config["conda_config"]["main"])
    shadow: "full"
    resources:
        cpus_per_task=8,
        mem_mb=32000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        mkdir -p $(dirname {output})
        # subsample the fastq files
        seqtk sample -s100 {input.R1} {wildcards.subsample_rate} > {input.R1}.subsampled.{wildcards.subsample_rate}.fastq
        seqtk sample -s100 {input.R2} {wildcards.subsample_rate} > {input.R2}.subsampled.{wildcards.subsample_rate}.fastq

        salmon quant -i {input.index} \
                    -l A \
                    -1 {input.R1}.subsampled.{wildcards.subsample_rate}.fastq \
                    -2 {input.R2}.subsampled.{wildcards.subsample_rate}.fastq \
                    --validateMappings \
                    -o {output} \
                    -p 8
        
        rm {input.R1}.subsampled.{wildcards.subsample_rate}.fastq
        rm {input.R2}.subsampled.{wildcards.subsample_rate}.fastq
        """

# Entire quantification worflow head
rule rarefraction_analysis:
    input:
        expand(
            [
                rules.salmon_quant_downsample.output[0],
                rules.oarfish_cov_rare_fraction_analysis.output[0]
            ],
            cell_line = sr_bulk_config['cell_lines'],
            sample = lr_bulk_config['sample_id'],
            subsample_rate = main_wf_config['rarefaction_rate']
        )
    output:
        touch(join(config['flag_dir'], "rarefraction_analysis.done"))
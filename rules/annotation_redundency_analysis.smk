from os.path import join
import glob

lr_bulk_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_bulk"]))
sr_bulk_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_bulk"]))
lr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["lr_sc_sn"]))
sr_sc_sn_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["sr_sc_sn"]))
main_wf_config = config_parser.load_configfile(join(config['main_wf_dir'],config["sub_wf_config"]["main_wf"]))

rule lr_bulk_minimap2_transcript:
    priority: 10
    input:
        fastq = lambda w: glob.glob(os.path.join(lr_bulk_config['samples_fastq_dir'][w.sample],f"{w.cell_line}.fastq*"))[0],
        ref = "/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_isoforms_{completeness}.fa"
    output:
        bam = temp(main_wf_config['tmp_dir'] + "/annotation_redundency_analysis/{sample}_{cell_line}_{completeness, (C|I|O)}.bam")
    resources:
        cpus_per_task=16,
        mem_mb=64000
    params:
        minimap2 = config["software"]["minimap2"],
        minimap2_trans_options = lambda w: lr_bulk_config["minimap2_trans_options"][w.sample]
    shell:
        """
        mkdir -p $(dirname {output})
        {params.minimap2} {params.minimap2_trans_options} -t {resources.cpus_per_task} {input.ref}  {input.fastq}  | samtools view -bS - > {output.bam}
        """


rule annotation_redundency_analysis_oarfish_cov:
    input:
        bam = rules.lr_bulk_minimap2_transcript.output,
    output:
        out_dir_cov = directory(os.path.join(main_wf_config['output_path'],"annotation_redundency_analysis/oarfish/{sample}/{completeness, (C|I|O)}/{cell_line}"))
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

# SR bulk
use rule sr_bulk_salmon_index as annotation_redundency_analysis_salmon_index with:
    input:
        "/vast/projects/LongBench/reference_files/SIRV_Set4_Norm_Sequences_20210507/SIRV_isoforms_{completeness}.fa"
    output:
        directory(join(main_wf_config['output_path'], "annotation_redundency_analysis/salmon/index_{completeness, (C|I|O)}"))

use rule sr_bulk_salmon_quant as annotation_redundency_analysis_salmon_quant with:
    input:
        R1 = rules.sr_bulk_fastp.output.R1,
        R2 = rules.sr_bulk_fastp.output.R2,
        index = rules.annotation_redundency_analysis_salmon_index.output
    output:
        directory(join(main_wf_config['output_path'], "annotation_redundency_analysis/salmon/salmon_quant/{completeness, (C|I|O)}/{cell_line}"))


rule annotation_redundency_analysis:
    input:
        expand(
            [
                rules.annotation_redundency_analysis_oarfish_cov.output[0],
                rules.annotation_redundency_analysis_salmon_quant.output[0]
            ],
            cell_line = sr_bulk_config['cell_lines'],
            sample = lr_bulk_config['sample_id'],
            completeness = ["C", "I", "O"]
        )
    output:
        touch(join(config['flag_dir'], "annotation_redundency_analysis.done"))
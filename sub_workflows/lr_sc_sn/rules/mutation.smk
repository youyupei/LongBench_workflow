# rule _find_mutations:
#     input:
#         bam_path="../results/flames_out/{x}/align2genome.bam",
#         allele_stat_path= lambda wildcards: f"../resources/{wildcards.x.split('_')[0]}_mutation/allele_stat.csv.gz",
#         barcode_file="../results/flames_out/{x}/whitelist.csv"
#     output:
#         output_rds ="../results/mutation_analysis/{x}_allele_stat.rds"
#     resources:
#         cpus_per_task=32,
#         mem_mb=1000000,
#         slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
#     params:
#         script = os.path.join(config['sub_wf_dir'],'scripts/mutation.R')
#     shell:
#         """
#         mkdir -p $(dirname {output.output_rds})
#         Rscript {params.script} {input.bam_path} {input.allele_stat_path} {input.barcode_file} {resources.cpus_per_task} {output.output_rds} 
#         """


rule run_clair3_rna:
    input:
        bam = join(results_dir, "PseudoBulkAlignment/Genome/{sample}_{cell_line}.sorted.bam"),
        bai = join(results_dir, "PseudoBulkAlignment/Genome/{sample}_{cell_line}.sorted.bam.bai"),
        ref_genome = config['reference']['genome']
    output:
        directory(join(results_dir, "clair3_rna/{sample}/{cell_line}"))
    params:
        preset = lambda wildcards: config['clair3_preset'][wildcards.sample],
    container: 
        "docker://hkubal/clair3-rna:latest"
    resources:
        cpus_per_task = 8,
        mem_mb = 64000
    shell:
        """
        mkdir -p {output}        
        /opt/bin/run_clair3_rna \
            -B {input.bam} \
            -R {input.ref_genome} \
            -o {output} \
            -t {resources.cpus_per_task} \
            -p {params.preset} \
            --enable_phasing_model \
            --include_all_ctgs \
            --conda_prefix /opt/conda/envs/clair3_rna
        """


rule variant_calling:
    input:
        expand(
            rules.run_clair3_rna.output,
            sample = config['sample_id'],
            cell_line = config['cell_line_list']
        )
rule _find_mutations:
    input:
        bam_path="../results/flames_out/{x}/align2genome.bam",
        allele_stat_path= lambda wildcards: f"../resources/{wildcards.x.split('_')[0]}_mutation/allele_stat.csv.gz",
        barcode_file="../results/flames_out/{x}/whitelist.csv"
    output:
        output_rds ="../results/mutation_analysis/{x}_allele_stat.rds"
    resources:
        cpus_per_task=32,
        mem_mb=1000000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        mkdir -p $(dirname {output.output_rds})
        Rscript scripts/mutation.R {input.bam_path} {input.allele_stat_path} {input.barcode_file} {resources.cpus_per_task} {output.output_rds} 
        """


rule get_mutations:
    input:
        expand("../results/mutation_analysis/{x}_allele_stat.rds", x=['scm_new', 'cll_new'])
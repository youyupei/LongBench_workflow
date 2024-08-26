rule demuxlet:
    """
    Run demuxlet on the single cell nuclei data.
    """
    input:
        expand(config['output_path'] + '/demuxlet/{sample}.demuxlet.best',
                sample=config['sample_name'])
    output:
        touch(os.path.join(results_dir, ".flag/demuxlet.done"))

rule prepare_input_vcf:
    """
    Take the mutation table and prepare the input VCF file for the cell lines.
    """
    input:
        config['reference']['cellline_mutations']
    output:
        os.path.join(results_dir, "demuxlet/DepmapMutTabConverted.vcf")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    params:
        script = os.path.join(config['sub_wf_dir'], 'scripts/DepmapMutTab2Vcf.py')
    shell:
        """
        OUT=$(dirname {output})/
        mkdir -p $OUT
        python3 {params.script} {input} {output}
        """

rule demuxlet_rule:
    """
    Run the demultiplexing tool on the input VCF file.
    """
    input:
        rules.cellranger.output,
        cr_out= os.path.join(results_dir, 'cellranger/{sample}/'),
        vcf = config['output_path'] + '/demuxlet/DepmapMutTabConverted.vcf'
    output:
        config['output_path'] + '/demuxlet/{sample}.demuxlet.best'
    resources:
        cpus_per_task=4,
        mem_mb=200000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        demuxlet_excutable = config['software']['demuxlet']
    shell:
        """
        # module load samtools
        module load htslib
        filename={output}
        mkdir -p $(dirname $filename)   

        {params.demuxlet_excutable} \
        --sam {input.cr_out}/outs/possorted_genome_bam.bam \
        --tag-group CB \
        --tag-UMI UB \
        --vcf {input.vcf} \
        --out ${{filename%.*}} \
        --field GT \
        --alpha 0 --alpha 0.5 
        touch {output}
        """

# ## Cluster level rules
# rule prepare_input_bam_cluster_level:
#     """
#     Add cluster as CLUSTER tag
#     Add BC_UMI as CB_UB tag
#     """
#     input:
#         inbam = config['output_path'] + '/flames_out/{sample}/align2genome.bam',
#         in_clus_csv = "/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/notebook/misc/{sample}_clusters_res{res}.csv"
#     output:
#         config['output_path'] + '/demuxlet/{sample}.res{res}.cluster_tag_added.bam'
#     resources:
#         cpus_per_task=2,
#         mem_mb=8000,
#         slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
#     threads:
#         16
#     params:
#         script="/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/workflow/scripts/DemuxletPreareBam_clusters.py",
#         python_dir="/stornext/Home/data/allstaff/y/you.yu/.cache/R/basilisk/1.14.0/FLAMES/1.9.1/flames_env/bin/python3"
#     shell:
#         """
#         {params.python_dir}  {params.script} {input.inbam} {input.in_clus_csv} {output} {threads}
#         """

# rule demuxlet_cluster_level:
#     """
#     Run the demultiplexing tool on the input VCF file.
#     """
#     input:
#         bam= config['output_path'] + '/demuxlet/{sample}.res{res}.cluster_tag_added.bam',
#         bam_sorted_idx= config['output_path'] + '/demuxlet/{sample}.res{res}.cluster_tag_added.bam.bai',
#         vcf = '/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/result/demuxlet/DepmapMutTabConverted.vcf'
#     output:
#         config['output_path'] + '/demuxlet/{sample}.res{res}.cluster.demuxlet.best'
#     resources:
#         cpus_per_task=4,
#         mem_mb=200000,
#         slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
#     shell:
#         """
#         module load samtools
#         module load htslib
#         filename={output}
#         mkdir -p $(dirname $filename)    

#         DEMUXLET="/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/software/demuxlet/bin/demuxlet"

#         "$DEMUXLET" \
#         --sam {input.bam} \
#         --tag-group CL \
#         --tag-UMI UB \
#         --vcf {input.vcf} \
#         --out ${{filename%.*}} \
#         --field GT \
#         --alpha 0 --alpha 0.5 
#         """

# ## Cluster level rules v2
# rule prepare_input_bam_cluster_level_contamination_filtered:
#     """
#     Add cluster as CL tag
#     Add BC_UMI as UB tag
#     """
#     input:
#         inbam = config['output_path'] + '/flames_out/{sample}/align2genome.bam',
#         in_clus_csv = "/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/notebook/misc/{sample}_clusters_res5.csv",
#         bc_list = '/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/notebook/misc/{sample}_decontX_contamination0.05.cell.csv'
#     output:
#         config['output_path'] + '/demuxlet/{sample}.cluster_v2_tag_added.bam'
#     resources:
#         cpus_per_task=1,
#         mem_mb=8000,
#         slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
#     params:
#         script="/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/workflow/scripts/DemuxletPreareBam_clusters.py",
#         python_dir="/stornext/Home/data/allstaff/y/you.yu/.cache/R/basilisk/1.14.0/FLAMES/1.9.1/flames_env/bin/python3"
#     threads:
#         16
#     shell:
#         """
#         {params.python_dir}  {params.script} {input.inbam} {input.in_clus_csv} {output} {threads}  {input.bc_list}
#         """
# rule demuxlet_cluster_level_contamination_filtered:
#     """
#     Run the demultiplexing tool on the input VCF file.
#     """
#     input:
#         bam= config['output_path'] + '/demuxlet/{sample}.cluster_v2_tag_added.bam',
#         bam_sorted_idx= config['output_path'] + '/demuxlet/{sample}.cluster_v2_tag_added.bam.bai',
#         vcf = '/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/result/demuxlet/DepmapMutTabConverted.vcf'
#     output:
#         config['output_path'] + '/demuxlet/{sample}.cluster_v2.demuxlet.best'
#     resources:
#         cpus_per_task=4,
#         mem_mb=200000,
#         slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
#     shell:
#         """
#         module load samtools
#         module load htslib
#         filename={output}
#         mkdir -p $(dirname $filename)    

#         DEMUXLET="/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/software/demuxlet/bin/demuxlet"

#         "$DEMUXLET" \
#         --sam {input.bam} \
#         --tag-group CL \
#         --tag-UMI UB \
#         --vcf {input.vcf} \
#         --out ${{filename%.*}} \
#         --field GT \
#         --alpha 0 --alpha 0.5 
#         """

# rule demuxlet_cluster_level_homo:
#     """
#     input VCF file: homo.
#     """
#     input:
#         bam= config['output_path'] + '/demuxlet/{sample}.cluster_v2_tag_added.bam',
#         bam_sorted_idx= config['output_path'] + '/demuxlet/{sample}.cluster_v2_tag_added.bam.bai',
#         vcf = '/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/result/demuxlet/homo_only_filtered.vcf'
#     output:
#         config['output_path'] + '/demuxlet/{sample}.cluster_homo_vcf.demuxlet.best'
#     resources:
#         cpus_per_task=4,
#         mem_mb=200000,
#         slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
#     shell:
#         """
#         module load samtools
#         module load htslib
#         filename={output}
#         mkdir -p $(dirname $filename)    

#         DEMUXLET="/home/users/allstaff/you.yu/LongBench/analysis/single_cell_nuclei/software/demuxlet/bin/demuxlet"

#         "$DEMUXLET" \
#         --sam {input.bam} \
#         --tag-group CL \
#         --tag-UMI UB \
#         --vcf {input.vcf} \
#         --out ${{filename%.*}} \
#         --field GT \
#         --alpha 0 --alpha 0.5 
#         """

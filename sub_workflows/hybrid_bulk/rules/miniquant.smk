from os.path import join


def _fmt(path_template, wildcards):
    return path_template.format(cell_line=wildcards.cell_line)


rule miniquant_hybrid_quant:
    input:
        ref = config["miniquant_config"]["reference_fa"],
        lr = lambda w: _fmt(config["miniquant_config"]["lr_fastq_pattern"][w.sample], w),
        sr1 = lambda w: _fmt(config["miniquant_config"]["sr_r1_pattern"], w),
        sr2 = lambda w: _fmt(config["miniquant_config"]["sr_r2_pattern"], w),
    output:
        out_dir = directory(join(config["output_path"], "miniquant_output/{sample}/{cell_line}"))
    resources:
        cpus_per_task=16,
        mem_mb=32000
    params:
        library_prep = lambda w: config["miniquant_config"]["long_reads_library_prep"][w.sample],
        strandness = config["miniquant_config"]["short_reads_strandness"],
        mem_gb = config["miniquant_config"]["mem_gb"]
    container: 
        "docker://tidesun/miniquant:latest" # March 30, 2026, using v1.4.1
    shell:
        """
        mkdir -p {output.out_dir}

        /app/miniQuant_linux/miniQuant quant \
          -r {input.ref} \
          -l {input.lr} \
          -1 {input.sr1} \
          -2 {input.sr2} \
          --long_reads_library_prep {params.library_prep} \
          --short_reads_strandness {params.strandness} \
          -t {resources.cpus_per_task} \
          --mem {params.mem_gb} \
          -o {output.out_dir}
        """


rule miniquant_dataset_to_dge_bundle:
    input:
        quant_dirs = lambda w: expand(
            rules.miniquant_hybrid_quant.output.out_dir,
            sample=[w.sample],
            cell_line=config["cell_lines"]
        ),
        bulk_meta = config["bulk_meta"],
        tx2gene = config["tx2gene_rds"],
        wrapper = join(config["sub_wf_dir"], "scripts/wrapper_MiniQuant.R")
    output:
        rds = join(config["output_path"], "miniquant_output/miniquant_bulk_dge.{sample}.rds")
    conda:
        config["conda"]["main"]
    resources:
        cpus_per_task=1,
        mem_mb=8000
    params:
        dataset_dir = join(config["output_path"], "miniquant_output/{sample}"),
        dataset_name = "{sample}"
    shell:
        """
        Rscript {input.wrapper} {params.dataset_dir} {params.dataset_name} {input.bulk_meta} {output.rds} {input.tx2gene}
        """


rule miniquant_all_dge:
    input:
        expand(
            rules.miniquant_dataset_to_dge_bundle.output.rds,
            sample=config["sample_id"]
        )
    output:
        touch(join(config["output_path"], ".flag/miniquant_all_dge.done"))

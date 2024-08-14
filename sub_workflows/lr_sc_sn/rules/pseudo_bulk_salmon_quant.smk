rule get_cell_line_bc_list:
    input:
        #expand(os.path.join(results_dir,  "reports/RDS/{sample}_annotated.rds"), sample = config['sample_id']) 
        lambda w: os.path.join(results_dir,  f"reports/RDS/{'sc' if 'sc' in w.sample else 'sn'}/{w.sample}_annotated.rds")
    output:
        [os.path.join(results_dir,  "misc/{sample}/", f"{x}_BC_list.txt") for x in config['cell_line_list']]
    shell:
        """
        module load R
        Rscript -e '
            library(dplyr)
            outdir <- dirname("{output[0]}")
            dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
            so <- readRDS("{input}")
            metadata <- so@meta.data

            # Create a list where each unique cell line has its associated barcodes
            barcode_list <- metadata %>%
                tibble::rownames_to_column(var = "barcode") %>%
                dplyr::rename(cell_lines = cell_lines) %>%
                dplyr::group_by(cell_lines) %>%
                dplyr::summarise(BC_list = list(barcode))

            # Convert the resulting tibble to a named list
            barcode_list <- setNames(barcode_list$BC_list, barcode_list$cell_lines)
            for (cell_line in names(barcode_list)) {{writeLines(barcode_list[[cell_line]], con = paste0(outdir, "/", cell_line, "_BC_list.txt"))}}
            '
        """


rule get_cell_line_pseudo_bulk_fq:
    input:
        bc_list =os.path.join(results_dir,  "misc/{sample}/{cell_line}_BC_list.txt"),
        fastq = os.path.join(results_dir, 'flames_out/{sample}/matched_reads_dedup.fastq')
    output:
        results_dir + "/PseudoBulkAlignment/{sample}_{cell_line}_pseudo_bulk.fastq"
    shell:
        """
        FASTQ_FILE={input.fastq}
        BARCODE_FILE={input.bc_list}
        OUTPUT_FILE={output}

        grep -A 3 -E "^@($(paste -sd '|' $BARCODE_FILE))_" $FASTQ_FILE  | grep -v "^--$" > $OUTPUT_FILE 
        """


rule ont_bulk_minimap2_transcript:
    priority: 10
    input:
        fastq = results_dir + "/PseudoBulkAlignment/{sample}_{cell_line}_pseudo_bulk.fastq",
        ref = config['reference']['transcript']
    output:
        bam = results_dir + "/PseudoBulkAlignment/{sample}_{cell_line}.bam"
    resources:
        cpus_per_task=32,
        mem_mb=64000
    params:
        minimap2 = config["software"]["minimap2"]
    shell:
        """
        module load samtools
        {params.minimap2} -ax lr:hq -t {resources.cpus_per_task} {input.ref}  {input.fastq} | samtools view -bS - > {output.bam}
        """


rule run_salmon:
    input:
        bam = results_dir + "/PseudoBulkAlignment/{sample}_{cell_line}.bam",
        ref = config["reference"]["transcript"]
    output:
        out_dir = directory(os.path.join(results_dir,"PseudoBulkSalmon/{sample}/{cell_line}"))
    conda:
        main_conda
    resources:
        cpus_per_task=32,
        mem_mb=32000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A", # auto
        gibbs_samples = 40
    shell:
        """
        mkdir -p {output.out_dir}
        salmon quant -t {input.ref} -l {params.lib_type} -a {input.bam} -p {resources.cpus_per_task}  -o {output.out_dir} --ont --numBootstraps {params.gibbs_samples}
        """

rule pseudo_bulk_salmon_quant:
    input: 
        expand(os.path.join(results_dir,"PseudoBulkSalmon/{sample}/{cell_line}"), 
                                                    cell_line = config['cell_line_list'],
                                                    sample = config['sample_id']) 
    output:
        touch(os.path.join(results_dir, ".flag/pseudo_bulk_salmon_quant.done"))

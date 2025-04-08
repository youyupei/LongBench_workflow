
# Geme quantification with featureCounts
rule featureCounts_gene_quant:
    input:
        bam = join(config["output_path"], "subjunc/bam/{cell_line}.sorted.bam"),
        gtf = config['reference']['gtf']
    output:
        rds = results_dir + "/featureCounts/{cell_line}_featureCounts.rds"
    resources:
        cpus_per_task=4,
        mem_mb=16000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        mkdir -p $(dirname {output.rds})
        Rscript -e "
            fc_SE <- Rsubread::featureCounts('{input.bam}',annot.ext='{input.gtf}', 
                isGTFAnnotationFile=TRUE, 

                GTF.featureType='gene', 
                GTF.attrType='gene_id', 
                useMetaFeatures=TRUE, ,
                isPairedEnd=TRUE,
                nthreads={resources.cpus_per_task})
            saveRDS(fc_SE, file='{output.rds}')
            "
        """

# Run salmon on the trimmed reads
## Step 2.1: Index the transcriptome
rule salmon_index:
    input: config['reference']['transcript']
    output: directory(join(results_dir, "salmon/index"))
    resources:
        cpus_per_task=16,
        mem_mb=32000
    conda:
        config["conda"]["main"]
    shell:
        """
        salmon index -t {input} -i {output} -k 31
        """

rule salmon_quant:
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
        index = rules.salmon_index.output
    output:
         directory(join(results_dir, "salmon/salmon_quant/{cell_line}"))
    conda:
        config["conda"]["main"]
    shell:
        """
        mkdir -p $(dirname {output})
        salmon quant -i {input.index} \
                    -l A \
                    -1 {input.R1} \
                    -2 {input.R2} \
                    --validateMappings \
                    -o {output} \
                    -p 16 \
                    --numBootstraps 50
        """

# Entire quantification workflow head
rule quantification:
    input:
        expand(
            [
                rules.salmon_quant.output[0],
                rules.featureCounts_gene_quant.output[0],
            ],
            cell_line = config['cell_lines']
        )
    output:
        touch(join(results_dir, "qc/.flag/salmon.done"))
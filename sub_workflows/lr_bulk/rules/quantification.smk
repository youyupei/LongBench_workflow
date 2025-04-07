
from os.path import join
results_dir = config["output_path"]
main_conda = config["conda"]["main"]
kallisto_output_dir = results_dir + "/kallisto_output"

# Geme quantification with featureCounts
rule featureCounts_gene_quant:
    input:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
        gtf = config['reference']['gtf']
    output:
        rds = results_dir + "/featureCounts/{sample}_{cell_line}_featureCounts.rds"
    resources:
        cpus_per_task=16,
        mem_mb=16000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        mkdir -p $(os.path.dirname{output.rds})
        out_dir=$(dirname {input.fastq})
        Rscript -e "
            fc_SE <- featureCounts('{input.bam}',annot.ext='{input.gtf}', 
                isGTFAnnotationFile=TRUE, 
                GTF.featureType='gene', 
                GTF.attrType='gene_id', 
                useMetaFeatures=TRUE, 
                nthreads={resources.cpus_per_task})
            saveRDS(fc_SE, file='{output.rds}')
            "
        """


# Transcritp quantification
## Salmon
rule salmon:
    input:
        bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam",
        ref = config["reference"]["transcript"]
    output:
        out_dir = directory(os.path.join(results_dir,"salmon_output/{sample}/{cell_line}"))
    conda:
        main_conda
    resources:
        cpus_per_task=32,
        mem_mb=100000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A", # auto
        gibbs_samples = 50
    shell:
        "salmon quant -t {input.ref} -l {params.lib_type} -a {input.bam} -p {resources.cpus_per_task}  -o {output.out_dir} --ont --numBootstraps {params.gibbs_samples}"


## Oarfish
### Oarfish no coverage module
rule oarfish_no_cov:
    input:
        bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam",
        ref = config["reference"]["transcript"]
    output:
        out_dir_nocov = directory(os.path.join(results_dir,"oarfish_nocov_output/{sample}/{cell_line}")),
    conda:
        config["conda"]["oarfish"]
    resources:
        cpus_per_task=16,
        mem_mb=64000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A", # auto
        gibbs_samples = 50
    priority: 100
    shell:
        """
        mkdir -p {output.out_dir_nocov}
        oarfish --alignments {input.bam} --threads {resources.cpus_per_task} --output {output.out_dir_nocov}/ -d . --filter-group no-filters --num-bootstraps 50
        """

### Oarfish  coverage module
rule oarfish_cov:
    input:
        bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam",
        ref = config["reference"]["transcript"]
    output:
        out_dir_cov = directory(os.path.join(results_dir,"oarfish_cov_output/{sample}/{cell_line}"))
    conda:
        config["conda"]["oarfish"]
    resources:
        cpus_per_task=16,
        mem_mb=64000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        lib_type = "A", # auto
        gibbs_samples = 50
    priority: 101
    shell:
        """
        mkdir -p {output.out_dir_cov}
        oarfish --alignments {input.bam} --threads {resources.cpus_per_task} --output {output.out_dir_cov}/ --model-coverage  -d . --filter-group no-filters --num-bootstraps 50
        """

## Kallisto

rule generate_t2g_and_cleaned_txfa:
    input:
        gtf = config["reference"]["gtf"],
        fa = config["reference"]["transcript"]
    output:
        out_t2g = join(kallisto_output_dir, "kallisto_bulk.t2g"),
        out_fa = temp(join(scratch_dir, "kallisto_bulk_transcript.fa"))
    params:
        script = join(config['sub_wf_dir'], "scripts/generate_kallisto_input.py")
    resources:
        cpus_per_task=8,
        mem_mb=20000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        "python3 {params.script} {input.gtf} {input.fa} {output.out_t2g} {output.out_fa}"


rule kallisto_index:
    input:
        transcriptome = rules.generate_t2g_and_cleaned_txfa.output.out_fa
    output:
        index_file = join(kallisto_output_dir, "kallisto_bulk.idx")

    resources:
        cpus_per_task=32,
        mem_mb=100000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        kallisto =  config["software"]["kallisto"],
        out_dir = kallisto_output_dir
    shell:
        """
        mkdir -p {params.out_dir}
        {params.kallisto} index -i {output.index_file} {input.transcriptome} -t {resources.cpus_per_task}
        """

rule kallisto_quantification:
    input:
        fastq = lambda w: os.path.join(config['samples_fastq_dir'][w.sample], "{cell_line}.fastq.gz"),
        kallisto_index = rules.kallisto_index.output.index_file,
        t2g = rules.generate_t2g_and_cleaned_txfa.output.out_t2g
    output:
        directory(join(kallisto_output_dir, "{sample}/{cell_line}"))
    params:
        kallisto =  config["software"]["kallisto"],
        bustools =  config["software"]["bustools"],
        out_dir = kallisto_output_dir
    resources:
        cpus_per_task=32,
        mem_mb=100000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        tech='bulk' 
        mkdir -p {output}

        {params.kallisto} bus -t {resources.cpus_per_task} --long --threshold 0.8 -x ${{tech}} -i {input.kallisto_index} -o {output} {input.fastq}

        {params.bustools} sort -t {resources.cpus_per_task} {output}/output.bus \
        -o {output}/sorted.bus; \
        {params.bustools} count {output}/sorted.bus \
        -t {output}/transcripts.txt \
        -e {output}/matrix.ec \
        -o {output}/count --cm -m \
        -g {input.t2g}

        # define the platform based on the sample wildcards
        # if wildcards.sample start with pb, then it is pacbio, else it is ONT
        if [[ {wildcards.sample} == pb* ]]; then
            platform='PacBio'
        else
            platform='ONT'
        fi

        {params.kallisto} quant-tcc -t {resources.cpus_per_task} \
            --long -P $platform -f {output}/flens.txt \
            {output}/count.mtx -i {input.kallisto_index} \
            -e {output}/count.ec.txt \
            -o {output} \
            --bootstrap-samples 50
        """

rule kallisto_convert:
    input:
        join(kallisto_output_dir, "{sample}/{cell_line}")
    output:
        join(kallisto_output_dir, ".{sample}_{cell_line}.converted")
    threads: 1
    resources:
        cpus_per_task=1,
        mem_mb=5000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        script = join(config['sub_wf_dir'], "scripts/lrk_to_abund.py")
    shell:
        """
        python {params.script} {input} && touch {output}
        """

# lr_bulk quantificaiton target rule
rule run_quantification:
    input:
        expand([
            rules.featureCounts_gene_quant.output[0],
            # rules.salmon.output[0],
            rules.oarfish_cov.output[0],
            # rules.kallisto_convert.output[0]]
            ],
            sample = config["sample_id"],
            cell_line = config["cell_lines"]),
    output:
        touch(results_dir + "/.flag/run_quantification.done")

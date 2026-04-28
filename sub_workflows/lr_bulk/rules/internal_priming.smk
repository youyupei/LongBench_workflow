import textwrap, os
from  os.path import  join


rule run_internal_priming_analysis:
    input:
        ".flag/{x}_run_primspotter.done"
    output:
        touch(".flag/{x}_internal_priming.done")

rule _internal_priming_identifier_single_run:
    input: 
        bam =  os.path.join(scratch_dir,"subsample_data/{sample}_{cell_line}/genome_map_3M.F2304.bam"),
        gtf = config['reference']['gtf_human'],
        genome = config['reference']['genome']
    output:
        summary=join(results_dir, "int_prim_analysis/{sample}_{cell_line}_summary.txt"),
        gene_level_counts=join(results_dir, "int_prim_analysis/{sample}_{cell_line}_gene_counts.tsv"),
        bam = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP_tag_added.bam"),
        bai = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP_tag_added.bam.bai")
    resources:
        cpus_per_task=32,
        mem_mb=32000
    params:
        python_script="/home/users/allstaff/you.yu/github/PrimeSpotter/PrimeSpotter/PrimeSpotter.py"
    conda:
        config["conda"]["PrimeSpotter"]
    shell:
        """
        mkdir -p $(dirname {output.summary})
        mkdir -p $(dirname {output.bam})
        # module load samtools
        
        python3 {params.python_script} --bam_file {input.bam} \
                                        --gtf_file {input.gtf} \
                                        --output-summary {output.summary} \
                                        --output-gene-count {output.gene_level_counts} \
                                        --genome-ref {input.genome} \
                                        --processes {resources.cpus_per_task}| samtools view -S -b | samtools sort > {output.bam}
        samtools index {output.bam}
        """


rule internal_priming_identifier:
    input:
        expand([join(results_dir, "int_prim_analysis/{sample}_{cell_line}_summary.txt"),
                join(results_dir, "int_prim_analysis/{sample}_{cell_line}_gene_counts.tsv")],
                sample = config['sample_id'],
                cell_line = config['cell_lines'])



rule internal_priming_split_bam:
    """
    Split the bam file into two: one with internal priming reads and one without.
    """
    input:
        bam = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP_tag_added.bam"),
        bai = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP_tag_added.bam.bai")
    output:
        temp_bam = [
            temp(join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP.bam")),
            temp(join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP.bam.bai")),
            temp(join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_nonIP.bam")),
            temp(join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_nonIP.bam.bai"))
        ]
    resources:
        cpus_per_task=1,
        mem_mb=8000
    shell:
        """
        # module load samtools
        samtools view -b --tag IP:T {input.bam} > {output.temp_bam[0]}
        samtools index {output.temp_bam[0]}
        samtools view -b --tag IP:F {input.bam} > {output.temp_bam[2]}
        samtools index {output.temp_bam[2]}

        """
rule internal_priming_quantification_featureCounts:
    input:
        bam_IP = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP.bam"),
        bai_IP = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_IP.bam.bai"),
        bam_nonIP = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_nonIP.bam"),
        bai_nonIP = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_nonIP.bam.bai"),
        gtf = config['reference']['gtf']
    output:
        rds = join(results_dir , "int_prim_analysis/featureCounts/{sample}_{cell_line}_featureCounts.rds")
    resources:
        cpus_per_task=16,
        mem_mb=16000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
        mkdir -p $(dirname {output.rds})
        Rscript -e "

            fc_gene_IP <- Rsubread::featureCounts('{input.bam_IP}',annot.ext='{input.gtf}',
                isGTFAnnotationFile=TRUE, 
                GTF.featureType='gene', 
                GTF.attrType='gene_id', 
                useMetaFeatures=TRUE, 
                isLongRead=TRUE,
                nthreads={resources.cpus_per_task})
            fc_gene_nonIP <- Rsubread::featureCounts('{input.bam_nonIP}',annot.ext='{input.gtf}',
                isGTFAnnotationFile=TRUE, 
                GTF.featureType='gene', 
                GTF.attrType='gene_id', 
                useMetaFeatures=TRUE, 
                isLongRead=TRUE,
                nthreads={resources.cpus_per_task})
            
            fc_exon_IP <- Rsubread::featureCounts('{input.bam_IP}',annot.ext='{input.gtf}',
                isGTFAnnotationFile=TRUE, 
                GTF.featureType='exon', 
                GTF.attrType='gene_id', 
                useMetaFeatures=TRUE, 
                isLongRead=TRUE,
                nthreads={resources.cpus_per_task})

            fc_exon_nonIP <- Rsubread::featureCounts('{input.bam_nonIP}',annot.ext='{input.gtf}',
                isGTFAnnotationFile=TRUE, 
                GTF.featureType='exon', 
                GTF.attrType='gene_id', 
                useMetaFeatures=TRUE, 
                isLongRead=TRUE,
                nthreads={resources.cpus_per_task})
            
        
            saveRDS(list(
                gene_IP = fc_gene_IP,
                gene_nonIP = fc_gene_nonIP,
                exon_IP = fc_exon_IP,
                exon_nonIP = fc_exon_nonIP
            ), file='{output.rds}')
            "
        """


rule internal_priming_Align_QC:
    input: 
        fa=config['reference']['genome'], 
        anno=config['reference']['gtf_gz'],
        genome_bam = join(scratch_dir, "int_prim_analysis/{sample}_{cell_line}_nonIP.bam"),
    resources:
        cpus_per_task=12,
        mem_mb=500000
    output:
        output = directory(os.path.join(results_dir, "qc/AlignQC/{sample}_{cell_line}_IP_filtered/")),
        tmp_dir = temp(directory(os.path.join(scratch_dir, "alignQC_tmp","{sample}_{cell_line}_IP_filtered")))
    priority: 10
    container: "docker://vacation/alignqc"
    shell:
        """
        mkdir -p $(dirname {output.output})
        mkdir -p {output.tmp_dir}
        alignqc analyze {input.genome_bam} \
            -g {input.fa} \
            --gtf {input.anno} \
            --output_folder {output.output} \
            --threads {resources.cpus_per_task} \
            --specific_tempdir {output.tmp_dir}
        """

rule internal_priming_quantification:
    input:
        # expand(rules.internal_priming_quantification_featureCounts.output.rds,
        #         sample = config['sample_id'],
        #         cell_line = config['cell_lines'])
        expand(rules.internal_priming_Align_QC.output[0],
                sample = config['sample_id'],
                # sample = ['pb_bulk'],
                cell_line = config['cell_lines'])


# ========================================================
# Full-sample internal priming filter pipeline
# ========================================================

rule internal_priming_identifier_full_sample:
    """
    Run PrimeSpotter on the full-sample genome BAM (primary alignments only, -F 2304).
    Outputs the IP-tagged BAM and summary files under int_prim_analysis/full/.
    """
    input:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
        gtf = config['reference']['gtf_human'],
        genome = config['reference']['genome']
    output:
        summary = join(results_dir, "int_prim_analysis/full/{sample}_{cell_line}_summary.txt"),
        gene_level_counts = join(results_dir, "int_prim_analysis/full/{sample}_{cell_line}_gene_counts.tsv"),
        bam = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_IP_tag_added.bam"),
        bai = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_IP_tag_added.bam.bai")
    resources:
        cpus_per_task=32,
        mem_mb=32000
    params:
        python_script="/home/users/allstaff/you.yu/github/PrimeSpotter/PrimeSpotter/PrimeSpotter.py"
    conda:
        config["conda"]["PrimeSpotter"]
    shell:
        """
        mkdir -p $(dirname {output.summary})
        mkdir -p $(dirname {output.bam})
        TMPBAM=$(mktemp --suffix=.bam)
        samtools view -b -F 2304 {input.bam} > $TMPBAM
        samtools index $TMPBAM
        python3 {params.python_script} --bam_file $TMPBAM \
                                        --gtf_file {input.gtf} \
                                        --output-summary {output.summary} \
                                        --output-gene-count {output.gene_level_counts} \
                                        --genome-ref {input.genome} \
                                        --processes {resources.cpus_per_task} | samtools view -S -b | samtools sort > {output.bam}
        samtools index {output.bam}
        rm $TMPBAM $TMPBAM.bai
        """


rule internal_priming_split_bam_full_sample:
    """
    Split the full-sample IP-tagged BAM into nonIP (retained) and IP (discarded).
    nonIP BAM is kept for downstream quantification; IP BAM is marked temp.
    """
    input:
        bam = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_IP_tag_added.bam"),
        bai = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_IP_tag_added.bam.bai")
    output:
        nonIP_bam = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_nonIP.bam"),
        nonIP_bai = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_nonIP.bam.bai"),
        IP_bam = temp(join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_IP.bam")),
        IP_bai = temp(join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_IP.bam.bai"))
    resources:
        cpus_per_task=4,
        mem_mb=8000
    shell:
        """
        samtools view -b --tag IP:T {input.bam} > {output.IP_bam}
        samtools index {output.IP_bam}
        samtools view -b --tag IP:F {input.bam} > {output.nonIP_bam}
        samtools index {output.nonIP_bam}
        """


rule internal_priming_filter_transcript_bam:
    """
    Filter the transcript BAM to retain only reads that are NOT internal-priming.
    Read names are extracted from the nonIP genome BAM and used to subset the
    existing transcript BAM. Output is unsorted (Oarfish does not require a
    sorted/indexed BAM).
    """
    input:
        nonIP_genome_bam = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_nonIP.bam"),
        nonIP_genome_bai = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_nonIP.bam.bai"),
        transcript_bam = results_dir + "/TranscriptAlignment/{sample}_{cell_line}.bam"
    output:
        filtered_bam = join(scratch_dir, "int_prim_analysis/full/{sample}_{cell_line}_nonIP_transcript.bam")
    resources:
        cpus_per_task=8,
        mem_mb=16000
    shell:
        """
        TMPIDS=$(mktemp)
        samtools view -F4 {input.nonIP_genome_bam} | cut -f1 | sort -u > $TMPIDS
        samtools view -b -N $TMPIDS {input.transcript_bam} > {output.filtered_bam}
        rm $TMPIDS
        """


rule internal_priming_full_sample:
    """Target rule: detect internal priming on all full samples."""
    input:
        expand(
            [join(results_dir, "int_prim_analysis/full/{sample}_{cell_line}_summary.txt"),
             join(results_dir, "int_prim_analysis/full/{sample}_{cell_line}_gene_counts.tsv")],
            sample=config['sample_id'],
            cell_line=config['cell_lines']
        )

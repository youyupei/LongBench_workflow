

clair3_tmp_sample_name_map = {
    'ont_bulk': 'ont',
    'pb_bulk': 'pacbio',
    'dRNA_bulk': 'drna',
}
rule  whatshap:
    input:
        ref = config['reference']['genome'],
        bam= results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
        vcf = lambda w: f"/home/users/allstaff/you.yu/LongBench/analysis/variant_allele_specific_analysis/Clair3-RNA/spikeins/all_contigs/{clair3_tmp_sample_name_map[w.sample]}/{w.cell_line}/output_enable_phasing.vcf.gz"
    output:
        vcf = results_dir + "/Mutation/{sample}_{cell_line}.phased.vcf"
    conda:
        config['conda']['whatshap']
    resources:
        cpus_per_task=8,
        mem_mb=64000
    shell:
        """
        whatshap phase -o {output.vcf} --ignore-read-groups --reference={input.ref} {input.vcf} {input.bam} 
        """

rule whatshap_stat:
    input:
        vcf = rules.whatshap.output.vcf
    output:
        stat = results_dir + "/Mutation/{sample}_{cell_line}.phased.vcf.stat"
    conda:
        config['conda']['whatshap']
    resources:
        cpus_per_task=2,
        mem_mb=4000
    shell:
        """
        whatshap stats {input.vcf} > {output.stat}
        """



rule get_genomic_region_bed:
    input:
        gtf = config['reference']['gtf_sorted'],
        genome = config['reference']['genome']
    output:
        exons_bed = results_dir + "/Mutation/exons.bed",
        genes_bed = results_dir + "/Mutation/genes.bed",
        introns_bed = results_dir + "/Mutation/introns.bed",
        intergenic_bed = results_dir + "/Mutation/intergenic.bed",
        two_col_genome =  temp(results_dir + "/Mutation/two_col_genome")
    resources:
        cpus_per_task=4,
        mem_mb=8000
    shell:
        """
        module load bedtools
        # Input
        GTF="{input.gtf}"
        GENOME="{input.genome}"
        # Output
        EXON_BED="{output.exons_bed}"
        GENE_BED="{output.genes_bed}"
        INTRON_BED="{output.introns_bed}"
        INTERGENIC_BED="{output.intergenic_bed}"


        # Make exon BED
        awk '$3 == "exon"' $GTF | \
        awk 'BEGIN{{OFS="\t"}}{{print $1, $4-1, $5}}' | \
        grep -E '^chr([1-9]|1[0-9]|20|X|Y)\s' | \
        bedtools sort -i - > $EXON_BED

        # Make gene BED
        awk '$3 == "gene"' $GTF | \
        awk 'BEGIN{{OFS="\t"}}{{print $1, $4-1, $5}}' | \
        grep -E '^chr([1-9]|1[0-9]|20|X|Y)\s' | \
        bedtools sort -i - > $GENE_BED

        # Make intron BED: gene - exon
        bedtools subtract -a $GENE_BED -b $EXON_BED | \
        bedtools sort -i - | bedtools merge -i - > $INTRON_BED

        # Make intergenic BED: genome - gene
        # First, create a two-column genome file
        cut -f1,2 $GENOME.fai > {output.two_col_genome}
        bedtools sort -i $GENE_BED -g {output.two_col_genome} | \
        bedtools complement -i - -g {output.two_col_genome} | \
        grep -E '^chr([1-9]|1[0-9]|20|X|Y)\s' > $INTERGENIC_BED
        """



rule genomic_coverage_analysis:
    input:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
        genome = config['reference']['genome'],
        exons_bed = results_dir + "/Mutation/exons.bed",
        genes_bed = results_dir + "/Mutation/genes.bed",
        introns_bed = results_dir + "/Mutation/introns.bed",
        intergenic_bed = results_dir + "/Mutation/intergenic.bed"
    output:
        cov = results_dir + "/Mutation/{sample}_{cell_line}.genomic_coverage.txt"
    params:
        min_coverage = 30
    resources:
        cpus_per_task=4,
        mem_mb=8000
    shell:
        """
        BAM="{input.bam}"
        OUT="{output.cov}"

        echo -e "Region\\tCovered\\tTotal" > $OUT
        samtools depth -a -b {input.exons_bed} $BAM  | awk '$3 >= {params.min_coverage} {{cov += 1}} {{total += 1}} END {{print "exon", cov, total}}' OFS="\\t" >> $OUT
        samtools depth -a -b {input.introns_bed} $BAM  | awk '$3 >= {params.min_coverage} {{cov += 1}} {{total += 1}} END {{print "intron", cov, total}}' OFS="\\t" >> $OUT
        samtools depth -a -b {input.intergenic_bed} $BAM  | awk '$3 >= {params.min_coverage} {{cov += 1}} {{total += 1}} END {{print "intergenic", cov, total}}' OFS="\\t" >> $OUT
        """



rule longcallR:
    input:
        bam = results_dir + "/GenomeAlignment/{sample}_{cell_line}.sorted.bam",
        reference = config['reference']['genome'],
        annotation = config['reference']['gtf_human']
    output:
        bam = join(scratch_dir, "LongcallR/{sample}_{cell_line}.longcallR.phased.sorted.bam"),
        bai = join(scratch_dir, "LongcallR/{sample}_{cell_line}.longcallR.phased.sorted.bam.bai"),
        vcf = join(results_dir, "LongcallR/{sample}_{cell_line}.longcallR.vcf")
    resources:
        cpus_per_task=32,
        mem_mb=300000
    params:
        preset = lambda wildcards: {'ont_bulk': "ont-cdna", 'pb_bulk': "hifi-masseq", 'dRNA_bulk': "ont-drna"}[wildcards.sample]
    shell:
        """
        mkdir -p $(dirname {output.bam})
        mkdir -p $(dirname {output.vcf})
        # the longcallR take prefix as output (remove the .bam suffix)
        prefix=$(dirname {output.bam})/$(basename {output.bam} .phased.sorted.bam)

        longcallR \
            --bam-path {input.bam} \
            --ref-path {input.reference} \
            --preset {params.preset} \
            --downsample  --downsample-depth 5000 \
            --min-depth 30 \
            --contigs chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
            --output $prefix    -t {resources.cpus_per_task}

        # after successful run, sort and index the bam file
        samtools sort -@ {resources.cpus_per_task} -o {output.bam} $prefix.phased.bam
        rm $prefix.phased.bam
        samtools index {output.bam}
        # move from scratch to project directory
        mv $prefix.vcf {output.vcf}
        """

rule longcallR_analsyis_ase:
    input:
        bam = rules.longcallR.output.bam,
        annotation = config['reference']['gtf_human']
    output:
        flag = touch(os.path.join(flag_dir, "Bulk_longcall.{sample}_{cell_line}.done"))
    params:
        script = '/stornext/General/data/user_managed/grpu_mritchie_1/Yupei/github/external/longcallR/allele_specific/longcallR-ase.py',    
        min_coverage = 30,
        output_prefix = lambda w: results_dir + f"/LongcallR/{w.sample}_{w.cell_line}.longcallR"
    resources:
        cpus_per_task=8,
        mem_mb=64000
    shell:
        """
        TMPDIR=/tmp python3 {params.script} \
            --bam {input.bam} \
            --annotation {input.annotation} \
            --min_support {params.min_coverage} \
            --output {params.output_prefix} \
            -t {resources.cpus_per_task} 
        touch {output.flag}
        """ 

rule longcallR_analsyis_asj:
    input:
        bam = rules.longcallR.output.bam,
        reference = config['reference']['genome'],
        #vcf = results_dir + "/Mutation/{sample}_{cell_line}.phased.vcf",
        annotation = config['reference']['gtf_human']
    output:
        flag = touch(os.path.join(flag_dir, "Bulk_longcall_asj.{sample}_{cell_line}.done"))
    params:
        script = '/stornext/General/data/user_managed/grpu_mritchie_1/Yupei/github/external/longcallR/allele_specific/longcallR-asj.py',    
        min_coverage = 30,
        output_prefix = lambda w: results_dir + f"/LongcallR/{w.sample}_{w.cell_line}.merged_genes_exons"
    resources:
        cpus_per_task=8,
        mem_mb=500000
    shell:
        """
        python3 {params.script} \
            --reference {input.reference} \
            --bam_file {input.bam} \
            --annotation_file {input.annotation} \
            --output_prefix {params.output_prefix} \
            -t {resources.cpus_per_task} \
            --min_sup  {params.min_coverage} 
        
        touch {output.flag}
        """ 

rule mutation_all:
    input:
        expand(
            [rules.whatshap.output.vcf,
             rules.whatshap_stat.output.stat,
             rules.genomic_coverage_analysis.output.cov,
             #rules.longcallR_analsyis_ase.output.flag,
             rules.longcallR_analsyis_asj.output.flag
             ],
            sample=config['sample_id'],
            cell_line=config['cell_lines']
        )
    output:
        touch(os.path.join(flag_dir, "Bulk_whatshapMutation.done"))
    conda:
        config['conda']['whatshap']
    shell:
        """
        echo "All Whatshap jobs completed successfully."
        """



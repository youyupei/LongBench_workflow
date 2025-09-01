rule Deepvariant_download_rnaseq_models:
    output:
        model = directory(os.path.join(scratch_dir, "Deepvariant_model")),
    shell:
        """
        mkdir -p {output.model}
        curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.data-00000-of-00001 > {output.model}/model.ckpt.data-00000-of-00001
        curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.example_info.json > {output.model}/example_info.json
        curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.index > {output.model}/model.ckpt.index
        curl https://storage.googleapis.com/deepvariant/models/DeepVariant/1.4.0/DeepVariant-inception_v3-1.4.0+data-rnaseq_standard/model.ckpt.meta > {output.model}/model.ckpt.meta
        """
        
rule Deepvariant_celline:
    input: 
        model = rules.Deepvariant_download_rnaseq_models.output.model,
        fa=config['reference']['genome'], 
        anno=config['reference']['gtf_gz'],
        genome_bam =join(results_dir,"subjunc/bam/{cell_line}.sorted.bam"),
        genome_bai =join(results_dir,"subjunc/bam/{cell_line}.sorted.bam.bai")
    # conda:
    #     config['conda']['AlignQC']
    resources:
        cpus_per_task=32,
        mem_mb=400000,
        slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au"
    output:
        directory(os.path.join(results_dir, "Deepvariants/{cell_line}"))
    priority: 10
    container: "docker://google/deepvariant:1.4.0"
    shell:
        """
        mkdir -p  {output}
        # Run DeepVariant.
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=WES \
        --customized_model={input.model}/model.ckpt \
        --vcf_stats_report=true \
        --ref={input.fa} \
        --reads={input.genome_bam} \
        --output_vcf={output}/output.vcf.gz \
        --output_gvcf={output}/output.g.vcf.gz \
        --intermediate_results_dir "{output}/intermediate_results_dir" \
        --make_examples_extra_args="split_skip_reads=true,channels=''" \
        --num_shards=32
        """



rule DeepVariant:
    input:
        expand(
            rules.Deepvariant_celline.output,
            cell_line=config['cell_lines']
        )
    output:
        touch(os.path.join(results_dir, "DeepVariant.done"))




rule  whatshap:
    input:
        ref = config['reference']['genome'],
        bam= join(results_dir,"subjunc/bam/{cell_line}.sorted.bam"),
        vcf = lambda w: f"/vast/projects/LongBench/analysis/variant_allele_specific_analysis/Clair3-illumina/spikeins/all_contigs/{w.cell_line}/merge_output.vcf.gz"
    output:
        vcf = results_dir + "/Mutation/{cell_line}.phased.vcf"
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
        stat = results_dir + "/Mutation/{cell_line}.phased.vcf.stat"
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
        bam= join(results_dir,"subjunc/bam/{cell_line}.sorted.bam"),
        genome = config['reference']['genome'],
        exons_bed = results_dir + "/Mutation/exons.bed",
        genes_bed = results_dir + "/Mutation/genes.bed",
        introns_bed = results_dir + "/Mutation/introns.bed",
        intergenic_bed = results_dir + "/Mutation/intergenic.bed"
    output:
        cov = results_dir + "/Mutation/ill_bulk_{cell_line}.genomic_coverage.txt"
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


rule mutation_all:
    input:
        expand(
            [rules.whatshap.output.vcf,
             rules.whatshap_stat.output.stat,
             rules.genomic_coverage_analysis.output.cov],
            cell_line=config['cell_lines']
        )
    conda:
        config['conda']['whatshap']
    shell:
        """
        echo "All Whatshap jobs completed successfully."
        """

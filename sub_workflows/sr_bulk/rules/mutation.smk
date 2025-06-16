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
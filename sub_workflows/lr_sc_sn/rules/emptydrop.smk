main_conda = config["conda"]["main"]

rule emptydrop:
    input:
        expand(config['output_path'] + '/.flag/{sample}_flames/emptydrop_gene_quantification.done',
                sample = config['sample_id'])
    output:
        touch(config['output_path'] + '/.flag/emptydrop.done')



rule get_empty_drop:
    """
    Get the empty drop cells from BLAZE
    * BC ranked at 20000-30000 but excluding BC that are <4 ED a way from the whitelist
    * Use High quality reads only
    * Assignment ED=0 for fast process. 
    """
    input:
        putative_bc_csv = config['output_path'] + '/flames_out/{sample}/putative_bc.csv'
    output:
        config['output_path'] + '/emptydrop/{sample}.empty_drops_list.csv'
    params:
        bc_rank = "20000-30000",
        minQ = 15
    run:
        import io
        from tqdm import tqdm
        import zipfile
        import pandas as pd
        from collections import Counter

        input_rank = [int(x) for x in params.bc_rank.split('-')]
        full_10x_wl = "/home/users/allstaff/you.yu/github/BLAZE/blaze/10X_bc/3M-february-2018.zip"
        out_emptydrop_fn = output[0]

        # read the putative
        dfs = pd.read_csv(input.putative_bc_csv, chunksize=1_000_000)
        raw_bc_count = Counter()
        for df in tqdm(dfs, desc = 'Counting high-quality putative BC', unit='M reads'):
            raw_bc_count += Counter(df[df.putative_bc_min_q >=params.minQ].putative_bc.value_counts().to_dict())
            
        # filter by 10X list
        whole_whitelist = []
        with zipfile.ZipFile(full_10x_wl) as zf:
            # check if there is only 1 file
            assert len(zf.namelist()) == 1
            with io.TextIOWrapper(zf.open(zf.namelist()[0]), encoding="utf-8") as f:
                for line in f:
                    whole_whitelist.append(line.strip())

        whole_whitelist = set(whole_whitelist)
        raw_bc_count = {k:v for k,v in raw_bc_count.items() if k in whole_whitelist}

        # get list of keys in raw_bc_count sorted by value in descending order
        ranked_bc = [k for k,v in sorted(raw_bc_count.items(), key=lambda x: x[1], reverse=True)]
        ept_bc = ranked_bc[input_rank[0]:input_rank[1]]

        with open(out_emptydrop_fn, 'w') as f:
            for k in ept_bc:
                f.write(k+'\n')

rule filter_empty_list:
    """
    Filter the empty drop list by excluding BC that are <4 ED a way from the whitelist
    """
    input:
        flag = rules.blaze.output,
        emp_list = config['output_path'] + '/emptydrop/{sample}.empty_drops_list.csv',
        whitelist = config['output_path'] + '/flames_out/{sample}/whitelist.csv'
    output:
        config['output_path'] + '/emptydrop/{sample}.empty_drops_list.filtered.csv'
    conda:
        main_conda
    shell:
        """
    python3 -c '
from fast_edit_distance import edit_distance
with open("{input.emp_list}") as el, open("{input.whitelist}") as wl:
    whitelist = [x.strip() for x in wl]
    emp_list = [x.strip() for x in el]
    emp_list = [x for x in emp_list if all([edit_distance(x, y) >= 4 for y in whitelist])]
with open("{output}", "w") as f:
    for x in emp_list:
        f.write(x + "\\n")
    '
        """

rule empty_drop_read_assignment:
    input:
        emp_list = config['output_path'] + '/emptydrop/{sample}.empty_drops_list.filtered.csv',
        putative_bc_csv = config['output_path'] + '/flames_out/{sample}/putative_bc.csv',
        fastq = lambda wildcards: config["samples_fastq_dir"][wildcards.sample]
    output:
        config['output_path'] + '/emptydrop/{sample}.empty_drops.fastq'
    conda:
        main_conda
    resources:
        cpus_per_task=32,
        mem_mb=20000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    shell:
        """
    python3 -c '
from blaze import read_assignment
read_assignment.assign_read(
    fastq_fns=["{input.fastq}"], # has to be a list
    fastq_out="{output[0]}",
    putative_bc_csv="{input.putative_bc_csv}", 
    whitelsit_csv="{input.emp_list}", 
    max_ed=0, 
    n_process={resources.cpus_per_task}, 
    batchsize=1000, 
    minQ=15)
    '
        """

rule empty_drop_run_flames_genome_mapping:
    input:
        config_file = config['flames_config']['genome_mapping'],
        fastq = config['output_path'] + '/emptydrop/{sample}.empty_drops.fastq',
        gtf = config['reference']['gtf'],
        genome = config['reference']['genome']
    output:
        bam = config['output_path'] + '/emptydrop/{sample}_flames/align2genome.bam',
        bai = config['output_path'] + '/emptydrop/{sample}_flames/align2genome.bam.bai'
    resources:
        cpus_per_task=32,
        mem_mb=200000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        minimap2 = config["software"]["minimap2"],
        k8 = os.path.dirname(config["software"]["minimap2"]) + '/k8'
    shell:
        """
        # module load samtools
        # module unload R
        # module load R/4.4.0

        out_dir=$(dirname {output.bam})
        mkdir -p $out_dir

        Rscript -e "
            library(FLAMES)
            config_file <- '{input.config_file}'
            fastq_flames= '{input.fastq}'
            outdir = '$out_dir'
            GTF = '{input.gtf}'
            genome = '{input.genome}'

            sce <- sc_long_pipeline(fastq=fastq_flames, 
                                    outdir=outdir, 
                                    annot=GTF, 
                                    genome_fa=genome, 
                                    config_file=config_file,
                                    barcodes_file='not used',
                                    minimap2='{params.minimap2}',
                                    k8='{params.k8}')
            "
        """


rule empty_drop_run_flames_gene_quantification:
    input:
        config_file = config['flames_config']['gene_quantification'],
        fastq = config['output_path'] + '/emptydrop/{sample}.empty_drops.fastq',
        bam = rules.empty_drop_run_flames_genome_mapping.output.bam,
        gtf = config['reference']['gtf'],
        genome = config['reference']['genome']
    output:
        flag = touch(config['output_path'] + '/.flag/{sample}_flames/emptydrop_gene_quantification.done'),
        gene_count = config['output_path'] + '/emptydrop/{sample}_flames/gene_count.csv'
    resources:
        cpus_per_task=32,
        mem_mb=200000,
        slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au"
    params:
        minimap2 = config["software"]["minimap2"],
        k8 = os.path.dirname(config["software"]["minimap2"]) + '/k8'
    shell:
        """
        # module load samtools
        # module unload R
        # module load R/4.4.0

        out_dir=$(dirname {output.gene_count})
        mkdir -p $out_dir

        Rscript -e "
            library(FLAMES)
            config_file <- '{input.config_file}'
            fastq_flames= '{input.fastq}'
            outdir = '$out_dir'
            GTF = '{input.gtf}'
            genome = '{input.genome}'

            sce <- sc_long_pipeline(fastq=fastq_flames, 
                                    outdir=outdir, 
                                    annot=GTF, 
                                    genome_fa=genome, 
                                    config_file=config_file,
                                    barcodes_file='not used',
                                    minimap2='{params.minimap2}',
                                    k8='{params.k8}')
            "
        """

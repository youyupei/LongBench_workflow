# BAM SAM and fastq
rule utilities_sort_and_index_bam:
    input:
        os.path.join(config['output_path'],"{x}.bam")
    output:
        bam = os.path.join(config['output_path'],"{x}.sorted.bam"),
        bai = os.path.join(config['output_path'],"{x}.sorted.bam.bai")
    resources:
        cpus_per_task=8,
        mem_mb=64000
        #slurm_extra="--mail-type=FAIL --mail-user=you.yu@wehi.edu.au --job-name=coverage_plot"
    shell:
        """
        samtools sort -@ {resources.cpus_per_task} -o {output.bam} {input}
        samtools index -@ {resources.cpus_per_task} {output.bam}
        """

# rule sam_to_bam:
#     input:
#         os.path.join(config['output_path'],"{x}.sam")
#     output:
#         os.path.join(config['output_path'],"{x}.bam")
#     shell:
#         """
#         module load samtools && samtools view -bS {input} > {output}
#         """


rule utilities_subsample_bam:
    input:
        os.path.join(config['output_path'],"{x}.bam")
    output:
        os.path.join(config['output_path'],"{x}.subsampled_{rate}.bam")
    resources:
        cpus_per_task=1,
        mem_mb=32000
    shell:
        "samtools view -s {wildcards.rate} -b {input} > {output}"


# rule bam_to_fastq:
#     input:
#         os.path.join(config['output_path'],"{x}.bam")
#     output:
#         os.path.join(config['output_path'],"{x}.fastq")
#     resources:
#         cpus_per_task=1,
#         mem_mb=32000
#     shell:
#         """
#         module load samtools 
#         samtools fastq {input} > {output}
#         """

# GTF and BED
rule utilities_gtf_to_bed:
    input:
        "{x}.gtf"
    output:
        "{x}.gtf.bed"
    resources:
        cpus_per_task=1,
        mem_mb=32000
    shell:
        """
        gxf2bed --input {input} --output {output}
        """
        
# github
rule utilities_git_check_commit:
    """
    usage: ask for .flag/gitrepo_{git_repo}.commit where git_repo takes the form of "user_repo"
    """
    input:
        os.path.join(config['output_path'],".flag/gitrepo_{git_repo}.exist")
    output:
        touch(os.path.join(config['output_path'],".flag/gitrepo_{git_repo}.commit"))
    params:
        local_dir = config['git_repo_dir'],
        git_repo = lambda w: w.git_repo.replace("_", '/')
    localrule: True
    shell:
        """
        output=$PWD/{output}
        cd {params.local_dir}/$(basename {params.git_repo})

        # Check if the current commit hash file exists, create if not
        if [ ! -f $output ]; then
            echo $(git rev-parse HEAD) > $output
            exit 0
        fi

        last_used_commit=$(cat $output)
        latest_commit=$(git ls-remote origin main | awk '{{print $1}}')

        if [ "$last_used_commit" != "$latest_commit" ]; then
            echo "New commit detected. Pulling latest changes..."
            git pull
            echo "$latest_commit" > $output
        fi
        """

rule utilities__git_clone:
    output:
        touch(os.path.join(config['output_path'],".flag/gitrepo_{git_repo}.exist"))
    params:
        local_dir = config['git_repo_dir'],
        git_repo = lambda w: w.git_repo.replace("_", '/')
    localrule: True
    shell:
        """
        mkdir -p {params.local_dir}
        cd {params.local_dir}

        # check directory exists
        if [ -d $(basename {params.git_repo}) ]; then
            echo "Directory $(basename {params.git_repo}) exists"
            cd $(basename {params.git_repo})
            git pull
        else
            git clone https://github.com/{params.git_repo}.git
        fi
        """

# conda
# rule create_conda_env_from_yml:
#     """
#     usage: create conda environment from yml file
#     """
#     input:
#         "{x}.yml"
#     output:
#         touch(".flag/conda_env_{x}.done")
#     params:
#         env_dir = config['envs_dir']
#     shell:
#         """
#         conda env create -f {input}
#         """
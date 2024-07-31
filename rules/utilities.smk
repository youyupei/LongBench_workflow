from snakemake.utils import min_version, _load_configfile
from smk_utils import *

# Set up config
global_config_path = _load_configfile("config/config.yaml")['global_config_path']
global_config = _load_configfile(global_config_path)
if config:  
    # the subworkflow config exisit and will overwrite the global config if conflicts
    config = update_config(global_config, config)
else:   
    config = global_config

    

# BAM SAM and fastq
rule sort_and_index_bam:
    input:
        "{x}.bam"
    output:
        "{x}.bam.bai"
    resources:
        cpus_per_task=8,
        mem_mb=64000
        #slurm_extra="--mail-type=END,FAIL --mail-user=you.yu@wehi.edu.au --job-name=coverage_plot"
    shell:
        """
        module load samtools/1.20
        mv {input} {input}.tmp.bam
        samtools sort -@ {resources.cpus_per_task} -o {input} {input}.tmp.bam
        samtools index -@ {resources.cpus_per_task} {input}
        rm {input}.tmp.bam
        """

# rule sam_to_bam:
#     input:
#         "{x}.sam"
#     output:
#         "{x}.bam"
#     priority: -100  # Lowest priority, only run if no other rules can generate the output
#     shell:
#         """
#         module load samtools && samtools view -bS {input} > {output}
#         """


rule subsample_bam:
    input:
        "{x}.bam"
    output:
        "{x}.subsampled_{rate}.bam"
    resources:
        cpus_per_task=1,
        mem_mb=32000
    shell:
        "module load samtools && samtools view -s {wildcards.rate} -b {input} > {output}"


rule bam_to_fastq:
    input:
        "{x}.bam"
    output:
        "{x}.fastq"
    resources:
        cpus_per_task=1,
        mem_mb=32000
    shell:
        """
        module load samtools 
        samtools fastq {input} > {output}
        """

# GTF and BED
rule gtf_to_bed:
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
rule git_check_commit:
    """
    usage: ask for .flag/gitrepo_{git_repo}.commit where git_repo takes the form of "user_repo"
    """
    input:
        ".flag/gitrepo_{git_repo}.exist"
    output:
        touch(".flag/gitrepo_{git_repo}.commit")
    params:
        local_dir = config['git_script_dir'],
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

rule _git_clone:
    output:
        touch(".flag/gitrepo_{git_repo}.exist")
    params:
        local_dir = config['git_script_dir'],
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
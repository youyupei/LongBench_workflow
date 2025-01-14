use rule split_fa as masked_genome_split_fa with:
    input:
        genome = config["reference"]["genome_masked"]
    output:
        temp(join(scratch_dir, "subjunc/genome.masked.fa"))

use rule subjunc_index as masked_genome_subjunc_index with:
    input:
        genome = rules.masked_genome_split_fa.output
    output:
        index = temp(directory(join(config["output_path"], "subjunc/genome_masked_index")))
    
use rule subjunc_mapping as masked_genome_subjunc_mapping with:
    input:
        R1 = rules.fastp.output.R1,
        R2 = rules.fastp.output.R2,
        index = rules.masked_genome_subjunc_index.output.index,
        gtf = config["reference"]["gtf"]
    output:
        bam = temp(join(config["output_path"], "subjunc/bam_masked/{cell_line}.bam"))

use rule sort_and_index_bam as masked_genome_sort_and_index_bam with:
    input:
        bam = rules.masked_genome_subjunc_mapping.output.bam
    output:
        sorted_bam = join(config["output_path"], "subjunc/bam_masked/{cell_line}.sorted.bam"),
        bai = join(config["output_path"], "subjunc/bam_masked/{cell_line}.sorted.bam.bai")

rule mapping_masked_genome:
    input:
        expand(
            [
                rules.masked_genome_sort_and_index_bam.output.sorted_bam,
                rules.masked_genome_sort_and_index_bam.output.bai
            ],
            cell_line = config['cell_lines']
        )
    output:
        touch(join(config["output_path"], ".flag/mapping_masked_genome.done"))
process mergeBams {
    publishDir params.project_dir + "/data/bam", mode: 'copy', overwrite:true

    label "small" // cpu: 4, memory: 16GB
    tag "$sample"

    input:
    tuple val(sample), path(bam_files)

    output:
    tuple val(sample), path(sampleBamFile), path(sampleBamIndex)

    script:
    sampleBamFile = sample.toString() + ".bam"
    sampleBamIndex = sampleBamFile + ".bai"
    """
    module load samtools
    samtools cat -o ${sampleBamFile} ${bam_files}
    samtools sort -@ 4 -o ${sampleBamFile} ${sampleBamFile}
    samtools index -@ 4 ${sampleBamFile}
    """
}

process modkitFullExtract {
    publishDir params.project_dir + "/data/modkit_extract", mode: 'copy', overwrite:true, pattern: "*.tsv.bgz"

    label "jumbo"
    tag "$sample $chr"

    input:
    tuple val(sample), path(sampleBamFile), path(sampleBamIndex), val(chr)

    output:
    tuple val(sample), path(sampleModkitExtractFile)

    script:
    sampleModkitExtractFile = sample.toString() + "_" + chr + ".tsv.bgz"
    """
    modkit extract full --bgzf --region ${chr} ${sampleBamFile} ${sampleModkitExtractFile}
    """
}

process modkitPileUp {
    label "jumbo"
    tag "$sample $chr"

    publishDir params.project_dir + "/data/modkit_pileup", mode: 'copy', overwrite:true, pattern: "*.bed"

    input:
    tuple val(sample), path(sampleBamFile), path(sampleBamIndex), val(chr)

    output:
    tuple val(sample), path(sampleModkitPileupFile)

    script:
    sampleModkitPileupFile = sample.toString() + "_" + chr + ".bed"
    """
    modkit pileup --region ${chr} ${sampleBamFile} ${sampleModkitPileupFile}
    """
}

workflow {
    // input files
    samples = Channel.fromPath("bam_samples.tsv")
        .splitCsv(sep: "\t", header: true)
        .map { row -> tuple(row.sample_name, file(row.bam_file)) }
        .groupTuple()

    // list of chromosomes that have > 100 counts in run1_H69
    chromosomes = Channel.of(
        "chr1",
        "chr2",
        "chr3",
        "chr4",
        "chr5",
        "chr6",
        "chr7",
        "chr8",
        "chr9",
        "chr10",
        "chr11",
        "chr12",
        "chr13",
        "chr14",
        "chr15",
        "chr16",
        "chr17",
        "chr18",
        "chr19",
        "chr20",
        "chr21",
        "chr22",
        "chrX",
        "chrY",
        "chrM",
        "ERCC-00002",
        "ERCC-00074",
        "ERCC-00096",
        "ERCC-00130",
        "SIRV1",
        "SIRV2",
        "SIRV3",
        "SIRV4",
        "SIRV5",
        "SIRV6",
        "SIRV7",
        "SIRV10001",
        "SIRV10003",
        "SIRV12001",
        "SIRV12003",
        "SIRV6001",
        "SIRV6002",
        "SIRV6003",
        "SIRV8001",
        "SIRV8003"
    )

    // main pipeline
    merged_bams = samples
        | mergeBams

    merged_bams.combine(chromosomes)
        | modkitFullExtract

    merged_bams.combine(chromosomes)
        | modkitPileUp
}
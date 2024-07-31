rule gtf_to_bed:
    input:
        {x}.gtf
    output:
        {x}.gtf.bed
    shell:
        """
        module load bedops
        gtf2bed --input={input} > {output}
        """

rule gff_to_gtf:
    input:
        {x}.gff
    output:
        {x}.gtf
    run:
        """
        from bioinfokit.analys import gff
        gff.gff_to_gtf(file="Athaliana_167_TAIR10.gene_chr1.gff3")
        """
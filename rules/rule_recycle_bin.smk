# Do not import this file directly. This is a recycle bin for unused rules.
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
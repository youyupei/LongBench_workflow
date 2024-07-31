import pandas as pd
import sys

# Define the file paths
input_file, output_file = sys.argv[1], sys.argv[2]

# Function to format the INFO field
def format_info(row):
    info_fields = [
        f"AF={row['af']}",
        f"DP={row['ref_count'] + row['alt_count']}",
        f"REF_COUNT={row['ref_count']}",
        f"ALT_COUNT={row['alt_count']}",
        f"GENE={row['gene']}",
        f"GENE_ID={row['ensembl_gene_id']}",
        f"HGNC_NAME={row['hgnc_name']}",
        f"PROTEIN_CHANGE={row['protein_change']}",
        f"MUTATION_EFFECT={row['mutation_effect']}",
    ]
    return ";".join(info_fields)

# Function to create VCF header
def create_vcf_header(df, sample_names):
    contigs = df['chrom'].unique()
    contig_headers = [f'##contig=<ID={contig}>' for contig in contigs]
    vcf_header = [
        '##fileformat=VCFv4.2',
        '##source=DepMap Portal',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        '##INFO=<ID=REF_COUNT,Number=1,Type=Integer,Description="Reference Count">',
        '##INFO=<ID=ALT_COUNT,Number=1,Type=Integer,Description="Alternate Count">',
        '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene Name">',
        '##INFO=<ID=GENE_ID,Number=1,Type=String,Description="Ensembl Gene ID">',
        '##INFO=<ID=HGNC_NAME,Number=1,Type=String,Description="HGNC Gene Name">',
        '##INFO=<ID=PROTEIN_CHANGE,Number=1,Type=String,Description="Protein Change">',
        '##INFO=<ID=MUTATION_EFFECT,Number=1,Type=String,Description="Mutation Effect">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'] + contig_headers + [
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(sample_names)
    ]
    return vcf_header

# Function to write VCF file
def write_vcf(df, output_file, sample_names):
    # sort the dataframe by chrom and pos, make chr10 comes after chr9
    df['chrom_n'] = df['chrom'].str.replace('chr', '')
    df['chrom_n'] = df['chrom_n'].replace('X', 23)
    df['chrom_n'] = df['chrom_n'].replace('Y', 24)
    df['chrom_n'] = df['chrom_n'].replace('M', 25)
    df['chrom_n'] = df['chrom_n'].astype(int)
    df = df.sort_values(['chrom_n', 'pos'])

    grouped = df.groupby(['chrom_n', 'pos'])
    vcf_header = create_vcf_header(df, sample_names)
    with open(output_file, 'w') as vcf:
        for line in vcf_header:
            vcf.write(line + '\n')
        for (chrom, pos), group in grouped:
            row = group.iloc[0]
            info = format_info(row)
            gt_data = {name: '0|0' for name in sample_names}  # Initialize all GT fields as missing
            for _, sub_row in group.iterrows():
                gt_data[sub_row['cell_line_display_name']] = sub_row['gt']
            
            genotype = "\t".join([gt_data[name] for name in sample_names])
            vcf.write(
                f"{row['chrom']}\t{pos}\t{row['mutation_id']}\t{row['ref']}\t{row['alt']}\t.\t.\t{info}\tGT\t{genotype}\n"
            )
    print(f"VCF file has been written to {output_file}")


if __name__ == '__main__':
    df = pd.read_csv(input_file)
    sample_names = df['cell_line_display_name'].unique()
    write_vcf(df, output_file, sample_names)
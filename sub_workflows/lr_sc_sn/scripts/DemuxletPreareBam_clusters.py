import pysam
import sys
from tqdm import tqdm
import pandas as pd
from collections import defaultdict


# Usage: ./<script name> <input bam file>  <cluster_csv> <output bam file name> <threads> <BC file, a BC per line (optional)>
bam_fn = sys.argv[1]
cluster_csv = sys.argv[2]
out_bam_fn = sys.argv[3]
threads = int(sys.argv[4])

# read in the cluster csv file
cluster_df = pd.read_csv(cluster_csv, index_col=0, header=0)
cluster_dict = cluster_df.iloc[:,0].to_dict()
cluster_dict = defaultdict(lambda: "NA", cluster_dict)

filter_bc = True if len(sys.argv) == 6 else False
if filter_bc:
    bc_fn = sys.argv[5]
    bcset = set()
    with open(bc_fn, "r") as f:
        for line in f:
            bcset.add(line.strip())

# read in the bam file
bamfile = pysam.AlignmentFile(bam_fn, "rb", threads=threads)
out_b = pysam.AlignmentFile(out_bam_fn, "wb", header=bamfile.header, threads=threads) 

# iterate over the reads, parse the reads name to get the tag after # and add them and a tag to the read  and write to a new bam file
for read in tqdm(bamfile):
    read_name = read.query_name
    bc, umi = read_name.split("#")[0].split("_")

    if filter_bc and bc not in bcset:
        continue

    read.set_tag("CL", str(cluster_dict[bc]), value_type='Z')
    read.set_tag("UB", bc+'_'+umi, value_type='Z')
    out_b.write(read)
out_b.close()


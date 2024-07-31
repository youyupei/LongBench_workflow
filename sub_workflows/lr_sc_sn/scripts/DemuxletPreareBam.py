import pysam
import sys
from tqdm import tqdm

# Usage: ./add_tag.py <input bam file>  <output bam file name> <thread> <read id file, a read id per line (optional)>
bam_fn = sys.argv[1]
out_bam_fn = sys.argv[2]
threads = int(sys.argv[3])

# read in the set of reads from txt file
filter_bc = True if len(sys.argv) == 5 else False
if filter_bc:
    bc_fn = sys.argv[4]
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
    
    read.set_tag("CB", bc, value_type='Z')
    read.set_tag("UB", umi, value_type='Z')
    out_b.write(read)
out_b.close()
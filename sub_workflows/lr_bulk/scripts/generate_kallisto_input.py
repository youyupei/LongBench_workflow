import pandas as pd
import sys, os
# usage python3 gtf2tx2gene.py <input.gtf> <output.tx2gene>
gtf, fa, out_t2g, out_fa = sys.argv[1:5]
buffer_size = 10000  # Define a reasonable buffer size


        

# Load GTF file
df = pd.read_csv(
    gtf,
    sep="\t",
    comment="#",
    header=None,
    usecols=[0, 2, 8],
    names=["chr", "feature", "attributes"],
)

# Filter for 'transcript' feature
df = df[df["feature"] == "transcript"]
# Extract transcript ID and gene ID from attributes
df["transcript_id"] = df["attributes"].str.extract('transcript_id "([^"]+)"')
df["gene_id"] = df["attributes"].str.extract('gene_id "([^"]+)"')
# remove .version from transcript_id
df["transcript_id"] = df["transcript_id"].str.split('.').str[0]

# Save to file
df = df[["transcript_id", "gene_id"]].dropna()


with open(fa, 'r') as f, open(out_fa, 'w') as fout:
    buffer = []  # Create a buffer to accumulate output
    lines = f.readlines()
    
    tx_id_lst = []
    for line in lines:
        if line.startswith('>'):
            tx_id = line.strip().split('.')[0][1:]
            tx_id_lst.append(tx_id)
            buffer.append('>' + tx_id + '\n')  # Add to buffer
        else:
            buffer.append(line)  # Add to buffer
        
        # Write buffer to file when it reaches the specified size
        if len(buffer) >= buffer_size:
            if set(tx_id_lst) - set(df.transcript_id):
                print(set(tx_id_lst) - set(df.transcript_id.values))
                os.r
            fout.write(''.join(buffer))
            buffer = []  # Clear buffer

    # Write any remaining lines in the buffer
    if buffer:
        fout.write(''.join(buffer))

# Save to file
df.to_csv(
    out_t2g, sep="\t", index=False, header=False
)

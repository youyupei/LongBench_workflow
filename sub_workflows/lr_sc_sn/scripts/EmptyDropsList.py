from tqdm import tqdm
import zipfile
import pandas as pd

input_rank = [int(x) for x in snakemake.params.bc_rank.split('-')]
full_10x_wl = "home/users/allstaff/you.yu/github/BLAZE/blaze/10X_bc/3M-february-2018.zip"
out_emptydrop_fn = snakemaker.output[0]

# read the putative
dfs = pd.read_csv(args.out_raw_bc_fn, chunksize=1_000_000)
raw_bc_count = Counter()
for df in tqdm(dfs, desc = 'Counting high-quality putative BC', unit='M reads'):
    raw_bc_count += Counter(df[
        df.putative_bc_min_q >=args.minQ].putative_bc.value_counts().to_dict())


# filter by 10X list
whole_whitelist = []
with zipfile.ZipFile("full_10x_wl") as zf:
    # check if there is only 1 file
    assert len(zf.namelist()) == 1
    with io.TextIOWrapper(zf.open(zf.namelist()[0]), 
                                            encoding="utf-8") as f:
        for line in f:
            whole_whitelist.append(line.strip())

whole_whitelist = set(whole_whitelist)
raw_bc_count = {k:v for k,v in raw_bc_count.items() if k in whole_whitelist}

# get list of keys in raw_bc_count sorted by value in descending order
ranked_bc = [k for k,v in sorted(raw_bc_count.items(), key=lambda x: x[1], reverse=True)]
ept_bc = ranked_bc[input_rank[0]:input_rank[1]]

with open(out_emptydrop_fn, 'w') as f:
    for k in ept_bc:
        f.write(k+'\n')

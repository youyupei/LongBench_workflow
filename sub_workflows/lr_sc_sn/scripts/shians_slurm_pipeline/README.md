# Dorado SLURM pipeline scripts

Scripts to create Modbam files from Nanopore POD5 output. The main benefit to this pipeline is that it breaks the dorado basecalling jobs into smaller blocks such that they can be efficiently parallelised onto the SLURM queue and not be killed by the time limit. Splitting POD5 files into blocks is done using symbolic links instead of copies, which is much faster and doesn't take up significant space.

The pipeline assumes that `dorado_basecall_block.sh` is in the `scripts` folder and can be called as `scripts/dorado_basecall_block.sh`.

## Usage
`sh run_dorado_basecall_6mA.sh <dir> <block_size>` where `dir` is the directory of the pod5 files and `block_size` is the number of pod5 files to proces in each batch. In testing it took  around 25 minutes per 100 pod5 file to call all context 6mA with sup accuracy, so setting `block_size` to around 100-200 is advised.

**Variables at the top of `run_dorado_basecall_6mA.sh` should be changed acoording to your needs**

Suppose you had the following directory structure
```
.
├── data
│   └── pod5
│       ├── folders full of pod5s...
│       └── ...
└── scripts
    ├── dorado_basecall_block.sh
    └── run_dorado_basecall_6mA.sh
```

Then you can run
```
sh run_dorado_basecall_6mA.sh data/pod5 200
```

If `dorado_basecall_block.sh` has been moved then please edit the `run_dorado_basecall_6mA.sh` with its location.

## Pipeline

The pipeline performs the following tasks:

1. Creates symlinks to blocks of pod5 files under `pod5_tmp`.
2. Submits a SLURM array job to process all blocks using `dorado` to create bam files in the `bam_tmp` folder.
3. Merges the output BAM files into a single file using `samtools merge`.
4. Demultiplexes using `dorado demux` into the `demux` folder.
5. (Optional) if `anno/sample_anno.csv` is detected as a dorado-compatible sample sheet then demultiplexed bam files will be renamed to their alias.
6. Index the bam files using `samtools index`.
7. Cleans up the temporary folders `*_tmp`.

This places the final resultant bam files in the `demux` folder.

## TODO
* Current pipeline is hard-coded for 6mA, need to generalise.
    * Generalisation is held back by dorado naming schemes, ideally we can define model and mod, then automatically download the correct model using `${MODEL}_${MOD}`, however there is not an easy way to figure out the latest available version. Could request user input, but this is cumbersome.

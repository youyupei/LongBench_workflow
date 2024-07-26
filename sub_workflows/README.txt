This is for storing workflows for processing individual dataset

Note: The Snakefile or rules defined in each dataset folder should be self-contained and should not depend on any other dataset. This is to ensure that the workflow can be run independently.
Universal rules that can be shared across datasets should be stored in the `../rules` folder.
├── dataset1/
│   ├── Snakefile
|   └── .../
├── dataset2/
│   ├── Snakefile
|   └── .../
└── .../
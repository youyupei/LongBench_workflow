# Code for  LongBench Project

The workflow was structure as follow:

* The processing and analysis of each individual dataset are in `sub_workflows/`. All of the will be automatically run when excuting the main workflow. 
* The config file for each sub workflow is stored in `config/`.
* Common rules that can be shared across the sub workflows were stored in `modules/`
* Analysis scripts can be find in `rmarkdown`


## Workflows:

* Main workflow:
    Mainly for combined analysis for different sub workflows.

* Sub workflows:
    * lr_bulk: preprocessing and analysis for long-read bulk samples, including the ONT and pacbio data
    * lr_sc_sn: preprocessing and analysis for long-read Single-cell and Single-nuclei samples, including the ONT and pacbio data
    * sr_sc_sn: preprocessing and analysis for short-read Single-cell and Single-nuclei samples
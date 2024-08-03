from snakemake.utils import min_version, _load_configfile
min_version("6.0")
configfile: "config/config.yaml"

#Sub workflows
module lr_sc_sn:
    snakefile: 'sub_workflows/lr_sc_sn/Snakefile'
    config:config
use rule * from lr_sc_sn as lr_sc_sn_*

module lr_bulk:
    snakefile: f'sub_workflows/lr_bulk/Snakefile'
    config:config
use rule * from lr_bulk as lr_bulk_*

# rule all:
#     input:
#         rules.lr_bulk_all.input,
#         rules.lr_sc_sn_all.input

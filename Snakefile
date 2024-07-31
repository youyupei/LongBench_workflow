from snakemake.utils import min_version, _load_configfile
min_version("6.0")
configfile: "config/config.yaml"

# # Sub workflows
ALL_SUBWFS = ['lr_sc_sn']


rule all:
    input:
        

for MODULE in ALL_SUBWFS:
    module MODULE:
        snakefile: f'sub_workflows/{MODULE}/Snakefile'
        config: _load_configfile(f'sub_workflows/{MODULE}/config/config.yaml')
    use rule * from MODULE

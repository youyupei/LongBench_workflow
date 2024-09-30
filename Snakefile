from snakemake.utils import min_version, _load_configfile
min_version("6.0")
configfile: "config/config.yaml"
sys.path.insert(0, config['main_wf_dir'] + "/config/")
import config_parser


EMAIL = config["email"]
REPORT = config["main_wf_dir"] + "/reports/snakemake_report.html"
RULEGRAPH = config["main_wf_dir"] + "/reports/rulegraph.pdf"
onsuccess:
    shell(
        """
        mail -s 'DONE: Snakemake main workflow' {EMAIL} < {log}
        # snakemake --unlock
        # snakemake --report {REPORT} all
        # snakemake --rulegraph all | dot -Tpdf > {RULEGRAPH}
        """)
onerror:
    shell("mail -s 'ERROR: Snakemake main workflow' {EMAIL} < {log}")



#Sub workflows
rule NoRuleSpecified:
    localrule: True
    shell:
        "echo 'There is no default rule for this workflow, please specify a rule to run, e.g. snakemake ... all'"

module lr_sc_sn:
    snakefile: 'sub_workflows/lr_sc_sn/Snakefile'
    config:config
use rule * from lr_sc_sn exclude utilities_gtf_to_bed as lr_sc_sn_*

module lr_bulk:
    snakefile: f'sub_workflows/lr_bulk/Snakefile'
    config:config
use rule * from lr_bulk as lr_bulk_*

module sr_sc_sn:
    snakefile: f'sub_workflows/sr_sc_sn/Snakefile'
    config:config
use rule * from sr_sc_sn as sr_sc_sn_*


# get sub-workflow config:
sub_wf_config = {}
for sub_wf in ['lr_sc_sn', 'lr_bulk']:
    sub_wf_config[sub_wf] = config_parser.sub_wf_config_parser(
        main_cfg_fn = config['global_config_path'], 
        sub_wf_name = sub_wf)


include: "rules/qc_plot.smk"
    
rule all:
    input:
        rules.lr_bulk_all.input,
        rules.lr_sc_sn_all.input,
        rules.sr_sc_sn_all.input,
        rules.combined_qc_plot.output

# Main workflow is mainly for combining the results from the sub-workflows and plotting

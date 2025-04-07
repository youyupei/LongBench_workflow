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


# Config the values for main workflow
scratch_dir = config['scratch_dir']

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

module sr_bulk:
    snakefile: f'sub_workflows/sr_bulk/Snakefile'
    config:config
use rule * from sr_bulk as sr_bulk_*

# get sub-workflow config:
sub_wf_config = {}
for sub_wf in ['lr_sc_sn', 'lr_bulk', 'sr_bulk', "sr_sc_sn"]:
    sub_wf_config[sub_wf] = config_parser.sub_wf_config_parser(
        main_cfg_fn = config['global_config_path'], 
        sub_wf_name = sub_wf)

include: "rules/qc_plot.smk"
# include: "rules/DE_analysis.smk"
include: "rules/sc_cell_line_anno.smk"
include: "rules/rarefraction_analysis.smk"
include: "rules/rarefraction_analysis_sc.smk"
include: "rules/annotation_redundency_analysis.smk"
include: "rules/rmarkdown.smk"


rule all:
    input:
        rules.lr_bulk_all.input,
        rules.lr_sc_sn_all.input,
        rules.sr_sc_sn_all.input,
        rules.sr_bulk_all.input,
        rules.main_qc_plot.output,
        rules.run_cell_line_annotation_pipeline.output,
        # rules.rmd_bulk_de_human.output,
        rules.knit_rmarkdown.input,
        rules.lr_sc_rarefraction_analysis.output
    default_target: True
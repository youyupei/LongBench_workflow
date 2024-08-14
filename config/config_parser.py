from snakemake.utils import min_version, _load_configfile
import os


# enable command line input. Not sure snakemake provides this functionality but I can do it by checking the modification
# between the config and configfile, any difference should be pase to the final config
# TOdo

def global_config_parser(main_cfg_fn):
    # parse the main config
    main_cfg = load_configfile(main_cfg_fn)
    root_dir = main_cfg['main_wf_dir']
    # convert relative paths to absolute paths
    main_cfg['git_repo_dir'] = os.path.join(root_dir, main_cfg['git_repo_dir'])
    
    for k in main_cfg['conda_config'].keys():
        main_cfg['conda'][k] = os.path.join(root_dir, main_cfg['conda_config'][k])
        main_cfg['conda_config'][k] = os.path.join(root_dir, main_cfg['conda_config'][k])
    return main_cfg


def sub_wf_config_parser(main_cfg_fn, sub_wf_name):
    # parse the main config
    main_cfg = global_config_parser(main_cfg_fn)
    root_dir = main_cfg['main_wf_dir']
    main_cfg['sub_wf_dir'] = os.path.join(root_dir, main_cfg['sub_wf_dir'][sub_wf_name])

    # parse the sub workflow config
    wf_cfg_fn = os.path.join(root_dir, main_cfg['sub_wf_config'][sub_wf_name])
    wf_cfg = load_configfile(wf_cfg_fn)

    config = update_config(main_cfg, wf_cfg)

    return config


def load_configfile(configfile):
    """Load a config file."""
    config = _load_configfile(configfile)
    return config

def update_config(config, new_config):
    """Update a config dictionary with a new one."""
    for key, value in new_config.items():
        if isinstance(value, dict):
            config[key] = update_config(config.get(key, {}), value)
        else:
            config[key] = value
    return config


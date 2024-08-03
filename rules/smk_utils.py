from snakemake.utils import min_version, _load_configfile


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


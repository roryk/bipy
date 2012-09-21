from bipy.log import logger, setup_logging
from bipy.cluster import start_cluster, stop_cluster
import yaml
import sys
from bipy.toolbox import fastqc


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    start_cluster(config)
    from bipy.cluster import view

    input_files = config["input"]
    for stage in config["run"]:
        if config["stage"][stage]["program"] == "fastqc":
            fastqc_config = config["stage"][stage]
            view.map(fastqc.run, input_files,
                     [fastqc_config] * len(input_files),
                     [config] * len(input_files))
    stop_cluster()

if __name__ == "__main__":
    main(*sys.argv[1:])

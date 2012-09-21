from bipy.toolbox import blastn
from bipy import cluster
from bipy.log import logger, setup_logging
import yaml
import sys


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    c = cluster.spawn_cluster(config)
    refs = config.get("ref", [])
    for ref in config.get("ref", []):
        blastn.run(config.get("query"), ref, config)
    c.stop()

if __name__ == "__main__":
    main(*sys.argv[1:])

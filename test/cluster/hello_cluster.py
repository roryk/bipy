#!/usr/bin/env python
import sys
import yaml
from bipy.cluster import start_cluster, stop_cluster
from bipy.log import setup_logging


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    start_cluster(config)
    from bipy.cluster import view
    setup_logging(config)
    from bipy.log import logger
    serial_result = map(lambda x: x**10, range(32))
    parallel_result = view.map(lambda x: x**10, range(32))
    logger.info(parallel_result)
    logger.info(serial_result == parallel_result)
    stop_cluster()

if __name__ == "__main__":
    main(sys.argv[1])

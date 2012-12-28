#!/usr/bin/env python
import sys
import yaml
from bipy.cluster import start_cluster, stop_cluster
from bipy.cluster import cluster_test
from bipy.log import setup_logging,  launch_iowatcher, setup_ipython_logging
import os
import sh
import subprocess
from IPython.parallel.util import interactive
import time

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    from bipy.log import logger
    start_cluster(config)
    from bipy.cluster import view, client
    connect_to = setup_ipython_logging()
    sh.python("/Users/rory/cache/bipy/scripts/bipy_logger.py", connect_to, "/Users/rory/cache/bipy/log/bipy_test.log", _bg=True)
    time.sleep(6)

    #launch_iowatcher(config)
    logger.info(client[:])
    client[:]['connect_to'] = connect_to
    client[:].execute('from bipy.log import setup_ipython_logging, setup_logging')
    client[:]['config'] = config
    client[:].execute('setup_ipython_logging(connect_to)')
    client[:].execute('from bipy.log import logger')
    #client[:].execute('setup_logging(config)')

    serial_result = map(lambda x: x**10, range(32))
    parallel_result = view.map(cluster_test, range(32), block=True)
    logger.info(parallel_result)
    logger.info(serial_result == parallel_result)
    stop_cluster()

if __name__ == "__main__":
    sys.path.append(os.path.split(sys.argv[0])[0])
    main(sys.argv[1])

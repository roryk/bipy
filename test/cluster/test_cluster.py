#!/usr/bin/env python
import yaml
from bipy.cluster import start_cluster, stop_cluster
from bipy.log import setup_logging
import unittest

CONFIG_FILE = "test/cluster/test_cluster.yaml"


def mappable_function(x):
    return x ** 10


def test_cluster():
    with open(CONFIG_FILE) as in_handle:
        config = yaml.load(in_handle)
    setup_logging(config)
    start_cluster(config)

    from bipy.cluster import view
    serial_result = map(mappable_function, range(32))
    parallel_result = view.map(mappable_function, range(32))
    assert(serial_result == parallel_result)

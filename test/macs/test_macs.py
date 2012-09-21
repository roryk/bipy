from bipy.log import setup_logging
from bipy.log import logger
from bipy.toolbox import macs
import sys
import yaml
import unittest
import os
from bcbio.utils import safe_makedir

STAGENAME = "macs"


class TestMacs(unittest.TestCase):

    def setUp(self):
        self.config_file = "test/macs/test_macs.yaml"
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)

        self.in_file = self.config["input"]
        self.stage_config = self.config["stage"][STAGENAME]
        self.options = self.stage_config["options"]

    def test_run(self):
        out_dir = os.path.join(self.config["dir"]["results"], STAGENAME)
        safe_makedir(out_dir)
        out_file = macs.run(self.in_file, self.options, None, out_dir)
        self.assertTrue(all(map(os.path.exists, out_file)))
        self.assertTrue(all([os.path.getsize(x) > 0 for x in out_file]))

    def test_run_with_config(self):
        out_file = macs.run_with_config(self.in_file, self.config,
                                        None,
                                        STAGENAME)
        self.assertTrue(all(map(os.path.exists, out_file)))
        self.assertTrue(all([os.path.getsize(x) > 0 for x in out_file]))

from itertools import repeat
from bipy.log import setup_logging
from bipy.utils import validate_config
from bipy.log import logger
from bipy.toolbox import cutadapt_tool
from bcbio.utils import file_exists
import sys
import yaml
import unittest

CONFIG_FILE = "test/cutadapt/test_cutadapt.yaml"


class TestCutadapt(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)
        self.input_files = self.config["input"]
        self.stage_config = self.config["stage"]["cutadapt"]

    def test_run(self):
        out_files = [cutadapt_tool.run(in_file, self.stage_config,
                                       self.config) for
                     in_file in self.input_files]
        self.assertTrue = all(map(file_exists, out_files))

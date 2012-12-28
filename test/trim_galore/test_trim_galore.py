import yaml
import unittest
from bipy.toolbox import trim
import os
import filecmp

STAGENAME = "trim_galore"


class TestTrimgalore(unittest.TestCase):

    def setUp(self):
        self.config_file = "test/trim_galore/test_trim_galore.yaml"
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)

        self.input_file = self.config["input"]
        self.gtf = self.config["annotation"]["file"]
        self.stage_config = self.config["stage"][STAGENAME]
        self.options = self.stage_config["options"]
        self.correct = self.config["correct"]

    def test_trim_galore(self):
        trim_galore = trim.TrimGalore(self.config)
        out_file = trim_galore(self.input_file)
        self.assertTrue(filecmp.cmp(out_file, self.correct))

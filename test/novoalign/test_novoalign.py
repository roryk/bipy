from bipy.toolbox import novoindex, novoalign
from bipy.utils import replace_suffix
import yaml
import sys
import unittest
import os

CONFIG_FILE = "test/novoalign/test_novoalign.yaml"


class TestNovoalign(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)
        self.input_files = self.config["input"]
        self.db = os.path.basename(replace_suffix(self.config["ref"], "nix"))
        self.db = os.path.join(self.config["dir"]["ref"], self.db)

    def test_novoindex(self):
        stage_config = self.config["stage"]["novoindex"]
        self.db = novoindex.run(self.config["ref"],
                                stage_config,
                                self.config)
        self.assertTrue(os.path.exists(self.db))
        self.assertTrue(os.path.getsize(self.db) > 0)

    def test_novoalign(self):
        stage_config = self.config["stage"]["novoalign"]
        for input_file in self.input_files:
            output_file = novoalign.run(input_file,
                                        self.db,
                                        stage_config,
                                        self.config)
            self.assertTrue(os.path.exists(output_file))
            self.assertTrue(os.path.getsize(output_file) > 0)

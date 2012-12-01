from bipy.toolbox.trim import Cutadapt
import yaml
import unittest
import filecmp

CONFIG_FILE = "test/cutadapt/test_cutadapt.yaml"


class TestCutadapt(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)

        self.stage_config = self.config["stage"]["cutadapt"]

    def test_single(self):
        correct_file = self.config["input_single_correct"]
        single = self.config["input_single"]
        cutadapt = Cutadapt(self.config)
        out_file = cutadapt(single)
        self.assertTrue(filecmp.cmp(correct_file, out_file))

    def test_pairedend(self):
        correct_files = self.config["input_paired_correct"]
        paired = self.config["input_paired"]
        cutadapt = Cutadapt(self.config)
        out_files = map(cutadapt, paired)
        self.assertTrue(all(map(filecmp.cmp, correct_files, out_files)))

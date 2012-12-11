from bipy.toolbox.tophat import Bowtie
import yaml
import unittest
import filecmp
import os
from bcbio.utils import file_exists

CONFIG_FILE = "test/bowtie/test_bowtie.yaml"


class TestBowtie(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)
        self.stage_config = self.config["stage"]["bowtie"]

    def test_single(self):
        correct_file = self.config["input_single_correct"]
        single = self.config["input_single"]
        bowtie = Bowtie(self.config)
        out_file = bowtie(single)
        self.assertTrue(file_exists(out_file))
        #self.assertTrue(filecmp.cmp(correct_file, out_file))
        os.remove(out_file)

    def test_pairedend(self):
        correct_file = self.config["input_paired_correct"]
        paired = self.config["input_paired"]
        bowtie = Bowtie(self.config)
        out_file = bowtie(paired)
        #self.assertTrue(filecmp.cmp(correct_file, out_file))
        self.assertTrue(file_exists(out_file))
        os.remove(out_file)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBowtie)
    unittest.TextTestRunner(verbosity=2).run(suite)

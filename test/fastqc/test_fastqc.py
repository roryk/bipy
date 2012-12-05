import yaml
from bipy.toolbox import fastqc
import unittest
import os
import filecmp
import shutil

cur_dir = os.path.dirname(os.path.realpath(__file__))
CONFIG_FILE = os.path.join(cur_dir, "test_fastqc.yaml")


class TestFastqc(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)
        self.stage = fastqc.FastQCStage(self.config)
        self.input = self.config["input"]

    def test_fastqc(self):
        correct_file = os.path.join(cur_dir, "data", "correct_fastqc.txt")
        run_result = self.stage(self.input)
        out_table = os.path.join(os.path.dirname(run_result),
                                 "test_fastqc_fastqc",
                                 "fastqc_data.txt")
        self.assertTrue(filecmp.cmp(correct_file, out_table))
        shutil.rmtree(os.path.join(cur_dir, "results"))

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFastqc)
    unittest.TextTestRunner(verbosity=2).run(suite)

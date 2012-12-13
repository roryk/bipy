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
        self.stage = fastqc.FastQC(self.config)

    def _get_result(self, run_result):
        return os.path.join(os.path.dirname(run_result),
                            "test_fastqc_fastqc",
                            "fastqc_data.txt")

    def test_fastqc(self):
        input_single = self.config["input_single"]
        correct_file = os.path.join(cur_dir, "data", "correct_fastqc.txt")
        run_result = self.stage(input_single)
        out_table = os.path.join(os.path.dirname(run_result),
                                 "test_fastqc_fastqc",
                                 "fastqc_data.txt")
        self.assertTrue(filecmp.cmp(correct_file, out_table))
        shutil.rmtree(os.path.join(cur_dir, "results"))

    def test_fastqc_paired(self):
        input_paired = self.config["input_paired"]
        correct_file = [os.path.join(cur_dir, "data", "correct_fastqc.txt"),
                        os.path.join(cur_dir, "data", "correct_fastqc.txt")]
        run_result = self.stage(input_paired)
        out_table = map(self._get_result, run_result)
        self.assertTrue(all(map(filecmp.cmp, correct_file, out_table)))
        shutil.rmtree(os.path.join(cur_dir, "results"))

    def test_fastqc_threads(self):
        config = self.config
        config["stage"]["fastqc"]["options"] = ["--threads", 2]
        stage = fastqc.FastQC(config)
        input_file = self.config["input_single"]
        correct_file = os.path.join(cur_dir, "data", "correct_fastqc.txt")
        run_result = stage(input_file)
        out_table = self._get_result(run_result)
        self.assertTrue(filecmp.cmp(correct_file, out_table))
        shutil.rmtree(os.path.join(cur_dir, "results"))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFastqc)
    unittest.TextTestRunner(verbosity=2).run(suite)

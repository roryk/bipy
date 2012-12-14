from bipy.toolbox.trim import Cutadapt
from bipy.toolbox.fastq import filter_reads_by_length
import yaml
import unittest
import filecmp
import os
from bipy.utils import append_stem
import shutil

CONFIG_FILE = "test/cutadapt/test_cutadapt.yaml"
CORRECT_DIR = "test/cutadapt/data/correct"


class TestCutadapt(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)

        self.stage_config = self.config["stage"]["cutadapt"]

    def _find_length_filter_correct(self, out_file):
        correct_dir = os.path.join(os.path.dirname(out_file), "correct",
                                   "length_filter")
        correct_file = os.path.join(correct_dir, os.path.basename(out_file))
        return correct_file

    def _cutadapt_single_correct(self, out_file):
        correct_dir = os.path.join(CORRECT_DIR, "cutadapt_single")
        return os.path.join(correct_dir, os.path.basename(out_file))

    def _cutadapt_paired_correct(self, out_file):
        correct_dir = os.path.join(CORRECT_DIR, "cutadapt_paired")
        return os.path.join(correct_dir, os.path.basename(out_file))

    def test_pairedend(self):
        paired = self.config["input_paired"]
        cutadapt = Cutadapt(self.config)
        out_files = cutadapt(paired)
        correct_files = map(self._cutadapt_paired_correct, out_files)
        self.assertTrue(all(map(filecmp.cmp, correct_files, out_files)))
        shutil.rmtree(os.path.dirname(out_files[0]))

    def test_single(self):
        single = self.config["input_single"]
        cutadapt = Cutadapt(self.config)
        out_file = cutadapt(single)
        correct_file = self._cutadapt_single_correct(out_file)
        self.assertTrue(filecmp.cmp(correct_file, out_file))
        os.remove(out_file)

    def test_length_filter(self):
        paired = self.config["input_paired"]
        out_files = filter_reads_by_length(paired[0], paired[1], min_length=20)
        correct_files = map(self._find_length_filter_correct, out_files)
        self.assertTrue(all(map(filecmp.cmp, correct_files, out_files)))
        map(os.remove, out_files)
        map(os.remove, [append_stem(x, "singles") for x in paired])

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCutadapt)
    unittest.TextTestRunner(verbosity=2).run(suite)

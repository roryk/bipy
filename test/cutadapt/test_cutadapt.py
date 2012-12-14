from bipy.toolbox.trim import Cutadapt
from bipy.toolbox.fastq import filter_reads_by_length
import yaml
import unittest
import filecmp
import os
from bipy.utils import append_stem

CONFIG_FILE = "test/cutadapt/test_cutadapt.yaml"


class TestCutadapt(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)

        self.stage_config = self.config["stage"]["cutadapt"]

    def _find_correct_file(self, out_file):
        correct_dir = os.path.join(os.path.dirname(out_file), "correct")
        correct_file = os.path.join(correct_dir, os.path.basename(out_file))
        return correct_file

    def test_single(self):
        correct_file = self.config["input_single_correct"]
        single = self.config["input_single"]
        cutadapt = Cutadapt(self.config)
        out_file = cutadapt(single)
        self.assertTrue(filecmp.cmp(correct_file, out_file))
        os.remove(out_file)

    def test_pairedend(self):
        correct_files = self.config["input_paired_correct"]
        paired = self.config["input_paired"]
        cutadapt = Cutadapt(self.config)
        out_files = cutadapt(paired)
        self.assertTrue(all(map(filecmp.cmp, correct_files, out_files)))
        map(os.remove, out_files)

    def test_length_filter(self):
        paired = self.config["input_paired"]
        out_files = filter_reads_by_length(paired[0], paired[1], min_length=20)
        correct_files = map(self._find_correct_file, out_files)
        self.assertTrue(all(map(filecmp.cmp, correct_files, out_files)))
        map(os.remove, out_files)
        map(os.remove, [append_stem(x, "singles") for x in paired])


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCutadapt)
    unittest.TextTestRunner(verbosity=2).run(suite)

from bipy import utils
import yaml
import unittest
import os
import tempfile
import sh

CONFIG_FILE = "test/utils/test_utils.yaml"
PAIRED_FILES = ["FA14dayrep1_1.fastq", "FA14dayrep1_2.fastq",
                "FA1dayrep1_1.fastq", "FA1dayrep1_2.fastq",
                "FA28dayrep1_1.fastq", "FA28dayrep1_2.fastq",
                "FA3dayrep3_1.fastq", "FA3dayrep3_2.fastq",
                "FA2dayrep1_1.fastq", "FA2dayrep1_2.fastq",
                "FA3dayrep1_1.fastq", "FA3dayrep1_2.fastq",
                "FA7dayrep1_1.fastq", "FA7dayrep1_2.fastq",
                "FA14dayrep2_1.fastq", "FA14dayrep2_2.fastq",
                "FA1dayrep2_1.fastq", "FA1dayrep2_2.fastq",
                "FA28dayrep2_1.fastq", "FA28dayrep2_2.fastq",
                "FA2dayrep2_1.fastq", "FA2dayrep2_2.fastq",
                "FA3dayrep2_1.fastq", "FA3dayrep2_2.fastq",
                "FA7dayrep2_1.fastq", "FA7dayrep2_2.fastq",
                "FA1dayrep3_1.fastq", "FA1dayrep3_2.fastq",
                "FA28dayrep3_1.fastq", "FA28dayrep3_2.fastq",
                "FA2dayrep3_1.fastq", "FA2dayrep3_2.fastq",
                "FA3dayrep3_1.fastq", "FA3dayrep3_2.fastq",
                "FA7dayrep3_1.fastq", "FA7dayrep3_2.fastq"]

CORRECT_DIR = "test/cutadapt/data/correct"


class TestUtils(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)

    def _touch_files(self, touched_files):
        test_dir = tempfile.mkdtemp()
        out_files = [os.path.join(test_dir, x) for x in touched_files]
        [open(x, "w").close() for x in out_files]
        return test_dir, out_files

    def test_locate(self):
        test_dir, out_files = self._touch_files(PAIRED_FILES)
        located_files = utils.locate("*.fastq", test_dir)
        self.assertTrue(all([x in out_files for x in located_files]))

    def test_combine_pairs(self):
        combined_files = utils.combine_pairs(PAIRED_FILES)
        correct_files = [list(x) for x in zip(PAIRED_FILES[0::2], PAIRED_FILES[1::2])]
        print combined_files
        print correct_files
        self.assertTrue(all([x in correct_files for x in combined_files]))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestUtils)
    unittest.TextTestRunner(verbosity=2).run(suite)

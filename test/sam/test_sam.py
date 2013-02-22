import yaml
from bipy.toolbox import sam
from bcbio.utils import file_exists, flatten
import unittest
import os
import hashlib
import tempfile


cur_dir = os.path.dirname(os.path.realpath(__file__))
CONFIG_FILE = os.path.join(cur_dir, "test_sam.yaml")


class TestSam(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)

    def test_disambiguate(self):
        in_files = self.config["input_bamdiff"]
        disambiguate = sam.Disambiguate(self.config)
        output = list(flatten(disambiguate(in_files)))
        out_md5 = map(self._get_md5, output)
        correct_files = self._correct_files(output)
        correct_md5 = map(self._get_md5, correct_files)
        self.assertTrue(out_md5 == correct_md5)
        #map(os.remove, output)

    def _get_md5(self, out_file):
        with open(out_file, "rb") as out_handle:
            data = out_handle.read()
        md5 = hashlib.md5(data)
        return md5.digest()

    def _correct_files(self, out_files):
        dirs = map(os.path.dirname, out_files)
        dirs = [os.path.join(x, "correct") for x in dirs]
        basenames = map(os.path.basename, out_files)
        return [os.path.join(d, f) for d, f in zip(dirs, basenames)]

    def test_get_reads_in_bamfile(self):
        bam_file = self.config["input_bamdiff"][0]
        reads = sam._get_reads_in_bamfile(bam_file)
        self.assertTrue(type(reads) == int)

    def test_downsample_bam_with_outfile(self):
        out_handle = tempfile.NamedTemporaryFile(suffix=".bam")
        bam_file = self.config["input_bamdiff"][0]
        target_reads = 2
        out_file = sam.downsample_bam(bam_file, target_reads, out_handle.name)
        self.assertEquals(sam._get_reads_in_bamfile(out_file), target_reads)

    def test_downsample_bam_with_memoize(self):
        bam_file = self.config["input_bamdiff"][0]
        target_reads = 2
        out_file = sam.downsample_bam(bam_file, target_reads)
        self.assertEquals(sam._get_reads_in_bamfile(out_file), target_reads)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSam)
    unittest.TextTestRunner(verbosity=2).run(suite)

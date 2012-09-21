from bipy.toolbox import bedtools
from bipy.utils import replace_suffix
from bcbio.utils import safe_makedir
import os
import unittest


class TestBedtools(unittest.TestCase):

    def setUp(self):
        self.bam_file = "test/data/s_1_1_10k.bam"
        self.gtf = "test/data/E_coli_k12.ASM584v1.15.gtf"

    def test_count_overlap(self):
        out_dir = "results/count_overlaps"
        safe_makedir(out_dir)
        out_file = os.path.join(out_dir, "count_overlap_test.counts")

        out_file = bedtools.count_overlaps(self.bam_file, self.gtf,
                                           out_file)
        self.assertTrue(os.path.exists(out_file))

    def test_multi_intersect(self):
        out_dir = "results/multi_intersect"
        safe_makedir(out_dir)
        out_file = os.path.join(out_dir, "multi_intersect_test.bed")
        bed_files = ["test/data/a.bed", "test/data/b.bed", "test/data/c.bed"]
        options = ["-header", "-cluster"]
        out_file = bedtools.multi_intersect(bed_files, options, out_file)

        self.assertTrue(os.path.exists(out_file))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBedtools)
    unittest.TextTestRunner(verbosity=2).run(suite)

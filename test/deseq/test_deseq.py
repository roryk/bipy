from bipy.toolbox import deseq
from bipy.utils import replace_suffix
from bcbio.utils import safe_makedir
import os
import unittest


class TestDeseq(unittest.TestCase):

    def setUp(self):
        self.count_file = "test/data/pasilla_gene_counts.tsv"
        self.conds = ["untreat", "untreat", "untreat", "untreat",
                      "treat", "treat", "treat"]

    def test_run(self):
        out_file = "results/deseq/test_deseq.txt"
        result = deseq.run(self.count_file, self.conds, out_file=out_file)
        self.assertTrue(out_file == result)
        self.assertTrue(os.path.getsize(result) > 0)

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDeseq)
    unittest.TextTestRunner(verbosity=2).run(suite)

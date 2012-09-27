from bipy.toolbox import dss
from bcbio.utils import safe_makedir, file_exists
import os
import unittest


class TestDss(unittest.TestCase):

    def setUp(self):
        self.count_file = "test/data/pasilla_gene_counts.tsv"
        self.conds = ["untreat", "untreat", "untreat", "untreat",
                      "treat", "treat", "treat"]

    def test_run(self):
        out_prefix = "results/tests/dss/test_dss"
        safe_makedir(os.path.dirname(out_prefix))
        result = dss.run(self.count_file, self.conds, out_prefix=out_prefix)
        self.assertTrue(file_exists(result))

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDss)
    unittest.TextTestRunner(verbosity=2).run(suite)

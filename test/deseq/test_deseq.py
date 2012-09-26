from bipy.toolbox import deseq
from bipy.utils import replace_suffix
from bcbio.utils import safe_makedir, file_exists
import os
import unittest


class TestDeseq(unittest.TestCase):

    def setUp(self):
        self.count_file = "test/data/pasilla_gene_counts.tsv"
        self.conds = ["untreat", "untreat", "untreat", "untreat",
                      "treat", "treat", "treat"]

    def test_run(self):
        out_prefix = "results/tests/deseq/test_deseq"
        safe_makedir(os.path.dirname(out_prefix))
        result = deseq.run(self.count_file, self.conds, out_prefix=out_prefix)
        self.assertTrue(file_exists(result))

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDeseq)
    unittest.TextTestRunner(verbosity=2).run(suite)

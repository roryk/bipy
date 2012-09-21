import yaml
import unittest
from bipy.toolbox import htseq_count
import os

STAGENAME = "htseq-count"


class TestHtseqcount(unittest.TestCase):

    def setUp(self):
        self.config_file = "test/htseq-count/test_htseq_count.yaml"
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)

        self.input_file = self.config["input"]
        self.gtf = self.config["annotation"]["file"]
        self.stage_config = self.config["stage"][STAGENAME]
        self.options = self.stage_config["options"]

    def test_run(self):
        out_file = "results/htseq-count/test_run.count"
        run_result = htseq_count.run(self.input_file,
                                     self.gtf,
                                     self.options,
                                     out_file)
        self.assertTrue(run_result == out_file)
        self.assertTrue(os.path.exists(run_result))
        self.assertTrue(os.path.getsize(run_result) > 0)

    def test_run_with_config(self):
        run_result = htseq_count.run_with_config(self.input_file,
                                                 self.config,
                                                 STAGENAME)
        self.assertTrue(os.path.exists(run_result))
        self.assertTrue(os.path.getsize(run_result) > 0)

    def test_combine(self):
        to_combine = self.config["to_combine"]
        out_file = "results/%s/combined_counts.counts" % (STAGENAME)
        df = htseq_count.combine_counts(to_combine, None, out_file=out_file)
        self.assertTrue(os.path.exists(out_file))
        self.assertTrue(os.path.getsize(out_file) > 0)

    def test_rpkm(self):
        rpkms = htseq_count.calculate_rpkm(self.input_file, self.gtf)
        rpkms.to_csv("results/%s/rpkms.txt" % (STAGENAME))



if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestHtseqcount)
    unittest.TextTestRunner(verbosity=2).run(suite)

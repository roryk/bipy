import yaml
import unittest
from bipy.toolbox import rseqc
from bcbio.utils import safe_makedir, file_exists
import os

STAGENAME = "rseqc"


class TestRseqc(unittest.TestCase):

    def setUp(self):
        self.config_file = os.path.join("test", STAGENAME,
                                        "test_" + STAGENAME + ".yaml")
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)

        self.input_file = self.config["input"]
        self.gtf = self.config["annotation"]["file"]
        self.stage_config = self.config["stage"][STAGENAME]

    def test_bam2bigwig(self):
        out_file = rseqc.bam2bigwig(self.input_file, self.config)
        self.assertTrue(file_exists(out_file))

    def test_bamstat(self):
        out_file = rseqc.bam_stat(self.input_file, self.config)
        self.assertTrue(file_exists(out_file))

    def test_clipping_profile(self):
        out_file = rseqc.clipping_profile(self.input_file, self.config)
        self.assertTrue(file_exists(out_file))

    def test_genebody_coverage(self):
        # test on the bam file
        out_file_bam = rseqc.genebody_coverage(self.input_file, self.config)
        print out_file_bam
        self.assertTrue(file_exists(out_file_bam))
        # test on the bigwig file
        #bigwig = rseqc.bam2bigwig(self.input_file, self.config)
        #out_file_bw = reseqc.genebody_coverage(bigwig, self.config)
        #self.assertTrue(file_exists(out_file_bw))

    def test_junction_annotation(self):
        """ XXX this needs an organism with actual splice juntions """
        out_file_junction = rseqc.junction_annotation(self.input_file,
                                                      self.config)
        print out_file_junction
        self.assertTrue(file_exists(out_file_junction))

    def test_junction_saturation(self):
        """ XXX this needs an organism with actual splice junctions """
        out_file_saturation = rseqc.junction_saturation(self.input_file,
                                                        self.config)
        print out_file_saturation
        self.assertTrue(file_exists(out_file_saturation))

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestRseqc)
    unittest.TextTestRunner(verbosity=2).run(suite)

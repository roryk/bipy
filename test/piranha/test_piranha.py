import yaml
import unittest
from bipy.toolbox import piranha
from bcbio.utils import safe_makedir, file_exists
import os
from bipy.toolbox import reporting
from bipy.utils import in2out

STAGENAME = "piranha"


class TestPiranha(unittest.TestCase):

    def setUp(self):
        self.config_file = os.path.join("test", STAGENAME, "test_" +
                                        STAGENAME + ".yaml")
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)

        self.in_file = self.config["input"]
        self.stage_config = self.config["stage"][STAGENAME]
        self.bin_size = self.stage_config.get("bin_size", 30)
        self.covariate = self.stage_config.get("covariate", None)
        self.out_dir = os.path.join("results", "test", STAGENAME)

    def test_run(self):

        exp_out = in2out(self.in_file, "piranha.bed", transform=True,
                         out_dir = self.out_dir)
        print exp_out
        out_file = piranha.run(self.in_file, bin_size=self.bin_size,
                               covariate=self.covariate, out_file=exp_out)

        print out_file
        self.assertTrue(file_exists(out_file))
        self.assertTrue(exp_out == out_file)
        os.remove(out_file)

    def test_stage_wrapper(self):
        exp_out = in2out(self.in_file, "piranha.bed", transform=True,
                         out_dir = self.out_dir)
        pstage = piranha.PiranhaStage(self.config)
        out_file = pstage(self.in_file)
        self.assertTrue(file_exists(out_file))
        self.assertTrue(exp_out == out_file)
        os.remove(out_file)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPiranha)
    unittest.TextTestRunner(verbosity=2).run(suite)

import yaml
import unittest
from bipy.toolbox import tophat
from bcbio.utils import file_exists
import os

STAGENAME = "tophat"


class TestTophat(unittest.TestCase):

    def setUp(self):
        config_file = "test/tophat/test_tophat.yaml"
        with open(config_file) as in_handle:
            self.config = yaml.load(in_handle)

        self.input_pairs = self.config["input"]

    def test_run_with_config(self):

        for input_files in self.input_pairs:
            if len(input_files) == 2:
                out_file = tophat.run_with_config(input_files[0],
                                                  input_files[1],
                                                  self.config["ref"],
                                                  "tophat",
                                                  self.config)
                print out_file
                self.assertTrue(file_exists(out_file))
            else:
                out_file = tophat.run_with_config(input_files[0],
                                                  None,
                                                  self.config["ref"],
                                                  "tophat",
                                                  self.config)
                print out_file
                self.assertTrue(file_exists(out_file))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTophat)
    unittest.TextTestRunner(verbosity=2).run(suite)

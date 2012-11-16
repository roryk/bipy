import yaml
from bipy.pipeline import GenericPipeline
import unittest
from bcbio.utils import file_exists
import os
from itertools import takewhile

BASE_DIR = os.path.dirname(__file__)
CONFIG_FILE = os.path.join(BASE_DIR, "test_pipeline.yaml")


class TestPipeline(unittest.TestCase):

    def setUp(self):
        with open(CONFIG_FILE) as in_handle:
            self.config = yaml.load(in_handle)
            self.input_files = self.config["input"]
            self.pipeline = GenericPipeline(self.input_files,
                                            self.config)

    def test_pipeline(self):
        for out_files in self.pipeline.process():
            print out_files

"""
    def test_pipeline_manually(self):
        stage = self.pipeline.to_run.get()
        print self.pipeline.view
        self.pipeline.view.map(stage, self.input_files)
        """
if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPipeline)
    unittest.TextTestRunner(verbosity=2).run(suite)

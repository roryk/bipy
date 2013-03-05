import yaml
import unittest
from bipy.toolbox import tophat
from bcbio.utils import file_exists
import os
import tempfile
import shutil

STAGENAME = "tophat"

TOPHAT_1_CONFIG = "test/tophat/test_tophat.yaml"
TOPHAT_2_CONFIG = "test/tophat/test_tophat2_bowtie2.yaml"

def _load_config(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    return config

def _run_fixture(input_files, config):
    if len(input_files) == 2:
        out_file = tophat.run_with_config(input_files[0], input_files[1],
                                          config["ref"], "tophat", config)
    else:
        out_file = tophat.run_with_config(input_files[0], None,
                                          config["ref"], "tophat", config)
    return out_file


class TestTophat(unittest.TestCase):

    def setUp(self):
        pass

    def test_run_with_config_tophat1(self):
        config = _load_config(TOPHAT_1_CONFIG)
        for input_files in config["input"]:
            out_file = _run_fixture(input_files, config)
            self.assertTrue(file_exists(out_file))

    def test_run_with_config_tophat2(self):
        config = _load_config(TOPHAT_2_CONFIG)
        for input_files in config["input"]:
            out_file = _run_fixture(input_files, config)
            self.assertTrue(file_exists(out_file))

    def test_run_with_prebuilt_transcripts_tophat2(self):
        temp_index = tempfile.NamedTemporaryFile()
        config = _load_config(TOPHAT_2_CONFIG)
        config["stage"]["tophat"]["options"] = {"transcriptome-index": temp_index.name}
        for input_files in config["input"]:
            out_file = _run_fixture(input_files, config)
            self.assertTrue(file_exists(out_file))

        # check to see that it runs with the index file too
        shutil.rmtree(os.path.dirname(config["dir"]["results"]))
        for input_files in config["input"]:
            out_file = _run_fixture(input_files, config)
            self.assertTrue(file_exists(out_file))


    def tearDown(self):
        config = _load_config(TOPHAT_2_CONFIG)
        out_dir = os.path.dirname(config["dir"]["results"])
        if os.path.exists(out_dir):
            shutil.rmtree(out_dir)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestTophat)
    unittest.TextTestRunner(verbosity=2).run(suite)

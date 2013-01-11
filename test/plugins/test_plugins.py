import yaml
import unittest
from bipy.pipeline.stages import AbstractStage
from bipy.plugins import StageRepository
from bcbio.utils import safe_makedir
import os
import tempfile

STAGENAME = "plugins"


class TestPlugins(unittest.TestCase):

    def setUp(self):
        self.config_file = os.path.join("test", STAGENAME,
                                        "test_" + STAGENAME + ".yaml")
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)
        self.repository = StageRepository(self.config)
        self.out_dir = os.path.join(self.config["dir"]["results"], "test",
                                    STAGENAME)
        safe_makedir(self.out_dir)

    def test_repository(self):
        """
        test that the repository object gets made

        """
        repository = StageRepository(self.config)
        self.assertTrue(isinstance(repository, StageRepository))

    def test_builtin_plugins(self):
        """
        test that we can find some of the built in plugins

        """
        fastqc = self.repository["fastqc"]
        self.assertTrue(issubclass(fastqc, AbstractStage))

    def test_custom_plugins(self):
        """
        test that the custom plugin directory is working

        """
        plugin_file = self._make_test_class()
        plugin_dir = os.path.dirname(plugin_file.name)
        repo = StageRepository({"dir": {"plugins": plugin_dir}})
        testplugin = repo["test_plugin"](self.config)
        self.assertTrue(testplugin("Test") == "Test")

    def _make_test_class(self):
        temp = tempfile.NamedTemporaryFile(suffix=".py")

        test_class = """
from bipy.pipeline.stages import AbstractStage
class TestPlugin(AbstractStage):

    stage = "test_plugin"

    def __init__(self, config):
        self.config = config

    def __call__(self, string):
        return string
"""
        temp.write(test_class)
        temp.flush()
        return temp

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPlugins)
    unittest.TextTestRunner(verbosity=2).run(suite)

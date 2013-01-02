import yaml
import unittest
from bipy.pipeline.stages import AbstractStage
from bipy.plugins import StageRepository
import os

STAGENAME = "plugins"


class TestPlugins(unittest.TestCase):

    def setUp(self):
        self.config_file = os.path.join("test", STAGENAME,
                                        "test_" + STAGENAME + ".yaml")
        with open(self.config_file) as in_handle:
            self.config = yaml.load(in_handle)
        self.repository = StageRepository(self.config)

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
        testplugin = self.repository["test_plugin"](self.config)
        self.assertTrue(testplugin("Test") == "Test")


class TestPlugin(AbstractStage):

    stage = "test_plugin"

    def __init__(self, config):
        self.config = config

    def __call__(self, string):
        return string

if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPlugins)
    unittest.TextTestRunner(verbosity=2).run(suite)

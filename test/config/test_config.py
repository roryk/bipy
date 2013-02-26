import unittest
import yaml
from bcbio.utils import file_exists
from bipy.utils import which
import logging
import os

logger = logging.getLogger('config_test')
logger.setLevel(logging.INFO)


class ConfigValidator(object):

    def __init__(self):
        pass

    def check_config(self, config):
        program_ok = self.check_program(config)
        return program_ok

    def check_exists(self, field, config):
        filename = config.get(field, None)
        if not filename:
            self._field_missing_warning(field)
            return True
        elif file_exists(filename):
            logger.info("%s location: %s" % (field, filename))
            return True
        else:
            logger.error("%s: %s was not found." % (field, filename))
            return False

    def check_program(self, config):
        need_dirs = ["picard"]
        programs = config.get("program", None)
        if not programs:
            self._field_missing_warning("program")
            return True

        for program, loc in programs.items():
            if program in need_dirs:
                if os.path.isdir(loc):
                    logger.info("%s lives in %s." % (program, loc))
                else:
                    self._directory_missing_error(program, loc)
                    return False
            else:
                if which(loc):
                    logger.info("%s will be run as %s." % (program, loc))
                else:
                    self._executable_missing_error(program, loc)
                    return False

        logger.info("The 'program' section of configuration seems valid.")
        return True

    def _directory_missing_error(self, program, loc):
        logger.error("%s is not a valid directory for %s." % (loc, program))

    def _executable_missing_error(self, program, loc):
        logger.error("%s is not a valid path for running %s." % (loc, program))

    def _field_missing_warning(self, field):
        logger.warning("%s field does not exist in the config file "
                       "and may cause downstream errors.")


class ConfigValidatorTests(unittest.TestCase):
    cur_dir = os.path.dirname(os.path.abspath(__file__))
    config_file = os.path.join(cur_dir, "test_config.yaml")

    def setUp(self):
        with open(self.config_file) as in_handle:
            self.correct_config = yaml.load(in_handle)
        self.validator = ConfigValidator()

    def test_config(self):
        self.assertTrue(self.validator.check_config(self.correct_config))

    def test_gtf(self):
        self.assertTrue(self.validator.check_exists("gtf", self.correct_config))

    def test_program(self):
        self.assertTrue(self.validator.check_program(self.correct_config))

    def test_wrong_program_dir(self):
        in_file = os.path.join(self.cur_dir, "test_wrong_program_dir.yaml")
        with open(in_file) as in_handle:
            config = yaml.load(in_handle)
        self.assertFalse(self.validator.check_config(config))

    def test_wrong_program_file(self):
        in_file = os.path.join(self.cur_dir, "test_wrong_program_file.yaml")
        with open(in_file) as in_handle:
            config = yaml.load(in_handle)
        self.assertFalse(self.validator.check_config(config))



def main():
    suite = unittest.TestLoader().loadTestsFromTestCase(ConfigValidatorTests)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == "__main__":
    main()

"""
general functions for handling common pipeline steps
"""

from bcbio.utils import safe_makedir, file_exists
from bipy.log import logger
import datetime
import os
import yaml
import abc


def setup_pipeline(config):
    """
    creates output directories
    and performs some minor validation on the configuration file
    """
    # make initial directories
    if "dir" not in config:
        logger.error("'dir' must be in config file, see example "
                     " configurations.")
        exit(-1)
    config = _setup_config(config)
    map(safe_makedir, config["dir"].values())

    _write_config(config)

    return config

def _write_config(config):
    """
    output a yaml file of the configuration for this run in the results
    directory
    """
    yaml_file = os.path.join(config["dir"].get("results", "results"),
                             "config.yaml")
    with open(yaml_file, "w") as out_handle:
        out_handle.write(yaml.dump(config))

def _setup_config(config):
    """
    configuration file setup
    """
    if "pipeline" in config:
        # add timestamp for the run
        if config["pipeline"].get("timestamp", None):
            time_fmt = '%Y-%m-%d-%H-%M-%S'
            now = datetime.datetime.now().strftime(time_fmt)
            stamped_dir = os.path.join(config["dir"]["results"], now)
            config["dir"]["results"] = stamped_dir

    return config


class AbstractRunner(object):

    @abc.abstractmethod
    def run(self):
        """Runs the command"""
        return

    @abc.abstractmethod
    def output(self):
        """Returns summary output of running the command in a dictionary"""
        return

    @abc.abstractmethod
    def command(self):
        """Returns the actual command line for running the command as a string"""
        return

    @abc.abstractmethod
    def accepts(self):
        """Returns a list of extensions the command accepts for input files"""
        return

    @abc.abstractmethod
    def next(self):
        """Returns the files that would be used post this stage of processing"""
        return


class AbstractStage(object):

    def __init__(self, config):
        self.config = config
        self._validate_config()

    def _validate_config(self):
        if "stage" not in self.config:
            raise ValueError('Could not find "stage" in the config file.')

    @property
    def stages(self):
        return self.config["stage"]

    def run_start_message(self, in_file):
        print "Starting lol. %s" % (self.stage)
        #logger.info("Starting %s on %s" % (self.stage, in_file))

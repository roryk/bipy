from bipy.toolbox import fastqc
from bipy.log import logger
from bcbio.utils import file_exists


class AbstractStage(object):
    """
    Wrapper around a stage for use with the pipeline manager. Subclasses should
    implement their own __init__ in addition to calling the super
    __init__. __init__ should parse the appropriate information from the config
    dict and save it in the objects state. __call__ should be implemented as a
    pure function with no closures.

    """

    def __init__(self, config):
        self.config = config
        self._validate_config()

    def __call__(self, in_file):
        pass

    def _validate_config(self):
        if "stage" not in self.config:
            raise ValueError('Could not find "stage" in the config file.')


class FastQCStage(AbstractStage):

    stage = "fastqc"

    def __init__(self, config):
        self.config = config
        super(FastQCStage, self).__init__(self.config)
        self.fastqc_config = config["stage"]["fastqc"]

    def _start_message(self, in_file):
        logger.info("Starting %s on %s" % (self.stage, in_file))

    def _end_message(self, in_file):
        logger.info("%s complete on %s." % (self.stage, in_file))

    def _check_run(self, in_file):
        if not file_exists(in_file):
            raise IOError('%s not found.' % (in_file))

    def __call__(self, in_file):
        self._start_message(in_file)
        self._check_run(in_file)
        out_file = fastqc.run(in_file, self.fastqc_config, self.config)
        self._end_message(in_file)
        return out_file


STAGE_LOOKUP = {"fastqc": FastQCStage}

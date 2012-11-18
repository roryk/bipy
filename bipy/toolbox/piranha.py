from bipy.pipeline.stages import AbstractStage
from bipy.utils import replace_suffix, in2out
from bcbio.utils import file_exists, safe_makedir
from sh import Piranha
from bipy.log import logger
import os


def run(in_file, bin_size=30, covariate=None, out_file=None):
    """
    takes a sorted BAM input file and runs Piranha on it
    with the specified bin size

    """
    if not out_file:
        out_file = replace_suffix(in_file, "piranha.bed")

    if file_exists(out_file):
        return out_file

    if covariate and file_exists(covariate):
        print "%s, %s, %s, %s" % (in_file, covariate, str(bin_size), out_file)
        Piranha("-s", in_file, covariate, b=bin_size, o=out_file)
    else:
        print "%s, %s, %s" % (in_file, str(bin_size), out_file)
        Piranha("-s", in_file, b=bin_size, o=out_file)

    return out_file


class PiranhaStage(AbstractStage):

    stage = "piranha"

    def __init__(self, config):
        self.config = config
        super(PiranhaStage, self).__init__(self.config)
        self.stage_config = config["stage"][self.stage]
        self.covariate = self.stage_config.get("covariate", None)
        self.bin_size = self.stage_config.get("bin_size", None)
        self.out_dir = os.path.join(config["dir"]["results"],
                                    self.stage)
        safe_makedir(self.out_dir)


    def _start_message(self, in_file):
        logger.info("Starting %s on %s." % (self.stage, in_file))

    def _end_message(self, in_file):
        logger.info("%s complete on %s." % (self.stage, in_file))

    def _check_run(self, in_file):
        if not file_exists(in_file):
            raise IOError("%s not found." % (in_file))

    def __call__(self, in_file):
        self._start_message(in_file)
        self._check_run(in_file)
        out_file = in2out(in_file, "piranha.bed", transform=True,
                          out_dir=self.out_dir)
        run(in_file, self.bin_size, self.covariate, out_file)
        self._end_message(in_file)
        return out_file

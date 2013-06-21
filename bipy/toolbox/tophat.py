from bcbio.ngsalign import tophat
import os
from bcbio.utils import file_exists, safe_makedir
from bipy.utils import remove_suffix, is_pair, get_in
from bipy.toolbox import fastqc
from bipy.pipeline.stages import AbstractStage
from bcbio.distributed.transaction import file_transaction
import sh
#from bipy.log import logger
from bcbio.log import logger
import sys


FASTQ_FORMAT_TO_BCBIO = {"sanger": None,
                         "illumina_1.3+": "illumina",
                         "illumina_1.5+": "illumina",
                         "illumina_1.8+": None,
                         "solexa": "solexa"}


def run_with_config(fastq_file, pair_file, ref_file,
                    stage_name, config):
    if not any_bowtie_reference_exists(ref_file):
        bowtie_reference_not_found_error()

    out_file = _bcbio_tophat_wrapper(fastq_file, pair_file,
                                     ref_file, stage_name,
                                     config)
    return out_file


def _bcbio_tophat_wrapper(fastq_file, pair_file, ref_file,
                          stage_name, config):
    bcbio_config = {}
    stage_config = config["stage"][stage_name]
    cores = stage_config.get("cores", 1)
    # use the listed quality format, if there isn't one, try to figure
    # out what format it is
    quality_format = stage_config.get("quality_format", None)
    if quality_format is None:
        fastq_format = fastqc.detect_fastq_format(fastq_file)
        quality_format = FASTQ_FORMAT_TO_BCBIO[fastq_format]

    max_errors = stage_config.get("max_errors", None)
    options = stage_config.get("options", {})
    tophat_loc = config["program"].get("tophat", "tophat")
    bowtie_loc = config["program"].get("bowtie", "bowtie")
    out_base = remove_suffix(os.path.basename(fastq_file))
    align_dir = os.path.join(config["dir"]["results"], stage_name)

    bcbio_config["resources"] = {"tophat": {"cores": cores,
                                            "options": options}}
    bcbio_config["algorithm"] = {}
    bcbio_config["program"] = {}
    bcbio_config["algorithm"]["quality_format"] = quality_format
    bcbio_config["algorithm"]["max_errors"] = max_errors
    bcbio_config["gtf"] = config.get("gtf", None)
    if bcbio_config["gtf"]:
        if not file_exists(bcbio_config["gtf"]):
            raise ValueError("GTF file does not exist. Please check to make sure that "
                             "the value of gtf is set corretly in the configuration file.")
            sys.exit(1)
    bcbio_config["program"]["tophat"] = tophat_loc
    bcbio_config["program"]["bowtie"] = bowtie_loc
    bcbio_config["program"]["picard"] = config["program"]["picard"]
    bcbio_config["program"]["gatk"] = {"dir": ""}

    out_file = tophat.align(fastq_file, pair_file, ref_file, out_base,
                            align_dir, bcbio_config)
    return out_file

class Tophat(AbstractStage):

    stage = "tophat"

    def __init__(self, config):
        self.config = config
        self.ref = config.get("ref", None)

    def __call__(self, in_file):
        self._start_message(in_file)
        if isinstance(in_file, basestring):
            logger.info("Detected %s as non-paired." % in_file)
            out_file = run_with_config(in_file, None, self.ref,
                                       self.stage, self.config)
        elif is_pair(in_file):
            logger.info("Detected %s as a pair." % in_file)
            out_file = run_with_config(in_file[0], in_file[1],
                                       self.ref, self.stage, self.config)
        else:
            logger.info("Detected %s as non-paired." % in_file)
            out_file = run_with_config(in_file[0], None, self.ref,
                                       self.stage, self.config)
        self._end_message(in_file)
        return out_file


class Bowtie(AbstractStage):

    stage = "bowtie"

    def __init__(self, config):
        self.config = config
        self.stage_config = config["stage"][self.stage]
        # defaults = {"q": True, "n": 2, "k": 1,
        #            "X": 2000, "best": True,
        #            "sam": True,
        #            "phred64-quals": True}
        #self.options = dict(defaults.items() +
        #                    self.stage_config.get("options", {}).items())
        self.options = self.stage_config.get("options", {})
        self.bowtie = sh.Command(self.stage_config.get("program", "bowtie"))
        self.out_prefix = os.path.join(get_in(self.config,
                                              ("dir", "results"), "results"),
                                              self.stage)
        self.ref_file = self.config["ref"]
        if not any_bowtie_reference_exists(self.ref_file):
            bowtie_reference_not_found_error()

    def _bowtie_se(self, in_file, out_file):
        self.bowtie(self.options, self.ref_file, in_file, out_file)

    def _bowtie_pe(self, in_file, out_file):
        self.bowtie(self.options, self.ref_file,
                    "-1", in_file[0], "-2", in_file[1], out_file)

    def _get_out_file(self, in_file):
        base, _ = os.path.splitext(os.path.basename(in_file))
        out_prefix = os.path.join(get_in(self.config,
                                         ("dir", "results"), "results"),
                                         self.stage)
        out_dir = os.path.join(out_prefix, base)
        out_file = os.path.join(out_dir, base + ".sam")
        return out_file

    def out_file(self, in_file):
        if is_pair(in_file):
            return self._get_out_file(in_file[0])
        else:
            return self._get_out_file(in_file)

    def __call__(self, in_file):
        self._start_message(in_file)
        out_file = self.out_file(in_file)

        if file_exists(out_file):
            return out_file

        with file_transaction(out_file) as tmp_out_file:
            if is_pair(in_file):
                self._bowtie_pe(in_file, tmp_out_file)
            else:
                self._bowtie_se(in_file, tmp_out_file)
        self._end_message(in_file)

        return out_file

def any_bowtie_reference_exists(prefix):
    return bowtie_1_reference_exists(prefix) or bowtie_2_reference_exists(prefix)

def bowtie_1_reference_exists(prefix):
    return file_exists(prefix + ".1.ebwt")
def bowtie_2_reference_exists(prefix):
    return file_exists(prefix + ".1.bt2")

def bowtie_reference_not_found_error():
    raise ValueError("Bowtie 1 or Bowtie 2 reference not found. Please make "
                     "sure that ref: is set correctly in the config file.")
    sys.exit(1)

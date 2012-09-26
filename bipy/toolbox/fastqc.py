""" tool to run fastqc on FASTQ/SAM/BAM files from HTS experiments """
import subprocess
from bipy.utils import flatten, remove_suffix
from bcbio.utils import safe_makedir
import os
import logging

logger = logging.getLogger("bipy")

_FASTQ_RANGES = {"sanger": [33, 73],
                 "solexa": [59, 104],
                 "illumina_1.3+": [64, 104],
                 "illumina_1.5+": [66, 104],
                 "illumina_1.8+": [33, 74]}


def detect_fastq_format(in_file, MAX_RECORDS=1000000):
    """
    detects the format of a fastq file
    will return multiple formats if it could be more than one
    """
    logger.info("Detecting FASTQ format on %s." % (in_file))
    kept = list(_FASTQ_RANGES.keys())
    with open(in_file) as in_handle:
        records_read = 0
        for i, line in enumerate(in_handle):
            # get the quality line
            if records_read >= MAX_RECORDS:
                break
            if i % 4 is 3:
                records_read += 1
                for c in line:
                    formats = kept
                    if len(formats) == 1:
                        return formats
                    for form in formats:
                        if (_FASTQ_RANGES[form][0] > ord(c) or
                            _FASTQ_RANGES[form][1] < ord(c)):
                            kept.remove(form)

    return formats


# list of module names for parsing the output files from fastqc
MODULE_NAMES = ["Basic Statistics", "Per base sequence quality",
                "Per sequence quality scores",
                "Per base sequence content",
                "Per base GC content",
                "Per sequence GC content",
                "Per base N content",
                "Sequence Length Distribution",
                "Overrepresented sequences"]


def _make_outdir(config):
    """ make the output directory "fastqc" where the data files live """
    outdir = os.path.join(config["dir"]["results"], "fastqc")
    safe_makedir(outdir)
    return outdir


def _make_outfile(input_file, config):
    outdir = _make_outdir(config)
    outfile = "".join([remove_suffix(os.path.basename(input_file)),
                       "_fastqc.zip"])
    return os.path.join(outdir, outfile)


def _build_command(input_file, fastqc_config, config):
    program = fastqc_config["program"]
    options = list(flatten(fastqc_config["options"]))
    outdir = _make_outdir(config)
    options += ["--outdir", outdir]
    cmd = list(flatten([program, options, input_file]))
    return cmd


def run(input_file, fastqc_config, config):
    outfile = _make_outfile(input_file, config)
    # if it is already done skip it
    if os.path.exists(outfile):
        return outfile

    cmd = _build_command(input_file, fastqc_config, config)
    subprocess.check_call(cmd)

    return outfile

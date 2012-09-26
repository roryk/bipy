import subprocess
from bipy.utils import append_stem
import os
import sh
from bipy.toolbox.fastqc import detect_fastq_format
import logging

logger = logging.getLogger("bipy")

_FASTQ_TYPE_TO_FLAG = {"sanger": "sanger",
                       "illumina_1.3+": "illumina",
                       "illumina_1.5+": "illumina",
                       "illumina_1.8+": "sanger",
                       "solexa": "solexa"}


def _get_quality_type(in_file):
    """ get fastq quality format. if multiple types are detected,
    pick the first one. no quality type is found assume sanger """
    fastq_format = detect_fastq_format(in_file)
    return _FASTQ_TYPE_TO_FLAG.get(fastq_format[0], "sanger")


def _get_length_cutoff(config):
    return config["stage"]["sickle"].get("length_cutoff", 20)


def _get_quality_cutoff(config):
    return config["stage"]["sickle"].get("quality_cutoff", 20)


def run_with_config(first, second=None, config=None):
    first_out = append_stem(first, "sickle")
    second_out = None
    if second:
        out_files = run_as_pe(first, second, config)
        return out_files

    else:
        out_file = run_as_se(first, config)
        return out_file


def run_as_pe(first, second, config):
    first_out = append_stem(first, "sickle")
    second_out = append_stem(second, "sickle")
    single_out = append_stem(first, "single")
    quality_type = _get_quality_type(first)
    length_cutoff = _get_length_cutoff(config)
    quality_cutoff = _get_quality_cutoff(config)
    if all(map(os.path.exists, [first_out, second_out, single_out])):
        return (first_out, second_out)
    sh.sickle("pe", f=first, r=second, l=length_cutoff, q=quality_cutoff,
              t=quality_type, o=first_out, p=second_out, s=single_out)
    return (first_out, second_out)


def run_as_se(first, config):
    first_out = append_stem(first, "sickle")
    pass


def run(in_file, end="se", qual="sanger", l="20", out_file=None):
    if not out_file:
        out_file = append_stem(in_file, "trimmed")

    if os.path.exists(out_file):
        return out_file

    cmd = ["sickle", end, "-f", in_file, "-o", out_file,
           "-t", qual, "-l", l, "-q", qual]

    subprocess.check_call(cmd)
    return out_file

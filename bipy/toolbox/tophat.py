from bcbio.ngsalign import tophat
import os
from bipy.utils import remove_suffix
from bipy.toolbox import fastqc

FASTQ_FORMAT_TO_BCBIO = {"sanger": None,
                         "illumina_1.3+": "illumina",
                         "illumina_1.5+": "illumina",
                         "illumina_1.8+": None,
                         "solexa": "solexa"}


def run_with_config(fastq_file, pair_file, ref_file,
                    stage_name, config):
    out_file = _bcbio_tophat_wrapper(fastq_file, pair_file,
                                     ref_file, stage_name,
                                     config)
    return out_file


def _bcbio_tophat_wrapper(fastq_file, pair_file, ref_file,
                          stage_name, config):
    bcbio_config = {}
    stage_config = config["stage"][stage_name]
    cores = config["cluster"].get("cores", None)
    # use the listed quality format, if there isn't one, try to figure
    # out what format it is
    quality_format = stage_config.get("quality_format", None)
    if quality_format is None:
        fastq_format = fastqc.detect_fastq_format(fastq_file)
        quality_format = FASTQ_FORMAT_TO_BCBIO[fastq_format]

    max_errors = stage_config.get("max_errors", None)
    tophat_loc = config["program"].get("tophat", "tophat")
    bowtie_loc = config["program"].get("bowtie", "bowtie")
    out_base = remove_suffix(os.path.basename(fastq_file))
    align_dir = os.path.join(config["dir"]["results"], stage_name)

    bcbio_config["resources"] = {"tophat": {"cores": cores}}
    bcbio_config["algorithm"] = {}
    bcbio_config["program"] = {}
    bcbio_config["algorithm"]["quality_format"] = quality_format
    bcbio_config["algorithm"]["max_errors"] = max_errors
    bcbio_config["gtf"] = config.get("gtf", None)
    bcbio_config["program"]["tophat"] = tophat_loc
    bcbio_config["program"]["bowtie"] = bowtie_loc

    out_file = tophat.align(fastq_file, pair_file, ref_file, out_base,
                            align_dir, bcbio_config)
    os.remove(out_file)

    out_dir = os.path.dirname(out_file)
    out_file_fixed = os.path.join(out_dir, out_base + ".sam")
    os.symlink("accepted_hits.sam", out_file_fixed)

    return out_file_fixed

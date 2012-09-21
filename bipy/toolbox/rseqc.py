import os
from bipy.utils import replace_suffix, which, remove_suffix, append_stem
import sh
from bipy.log import logger
from bcbio.utils import file_exists, safe_makedir
from os.path import basename


def _results_dir(config, prefix=""):

    def _make_dir(base):
        out_dir = os.path.join(base, "rseqc",
                               prefix)
        safe_makedir(out_dir)
        return out_dir

    if "dir" not in config:
        return _make_dir("results")
    if "results" not in config["dir"]:
        return _make_dir("results")
    else:
        return _make_dir(config["dir"]["results"])


def _fetch_chrom_sizes(config):
    if "annotation" not in config:
        logger.error("'annotation' must be in the yaml file. See example "
                     " configuration files")
        exit(1)
    if "genome" not in config["annotation"]:
        logger.error("'genome_name' must be in the yaml file under  "
                     " 'annotation'. See example configuration files.")
        exit(1)
    genome = config["annotation"]["genome"]
    chrom_size_file = os.path.join(_results_dir(config),
                                   genome + ".sizes")
    if file_exists(chrom_size_file):
        return chrom_size_file

    sh.fetchChromSizes(genome, _out=chrom_size_file)
    if not file_exists(chrom_size_file):
        logger.error("chromosome size file does not exist. Check "
                     "'annotation': 'genome' to make sure it is valid.")
        exit(1)
    return chrom_size_file


def bam2bigwig(in_file, config, out_file=None):
    """
    assumes the library preparation was not strand specific for now
    """
    prefix = "bigwig"
    chrom_size_file = config["annotation"].get("chrom_size_file", None)
    if not chrom_size_file:
        chrom_size_file = _fetch_chrom_sizes(config)
    out_prefix = os.path.join(_results_dir(config, prefix),
                              basename(in_file))
    wiggle_file = out_prefix + ".wig"
    bam2wig = sh.Command(which("bam2wig.py"))
    bam2wig(i=in_file, s=chrom_size_file, o=out_prefix)

    if out_file is None:
        out_file = out_prefix + ".bw"

    sh.wigToBigWig(wiggle_file, chrom_size_file, out_file)
    return out_file


def bam_stat(in_file, config, out_file=None):
    """
    dump read maping statistics from a SAM or BAM file to out_file
    """
    if not out_file:
        out_file = os.path.join(_results_dir(config, "bam_stat"),
                                replace_suffix(os.path.basename(in_file),
                                               "bam_stat.txt"))
    bam_stat_run = sh.Command(which("bam_stat.py"))
    bam_stat_run(i=in_file, _err=out_file)

    return out_file


def clipping_profile(in_file, config, out_prefix=None):
    """
    estimate the clipping profile of the reads
    """
    if not out_prefix:
        out_prefix = os.path.join(_results_dir(config, "clipping"),
                                  append_stem(os.path.basename(in_file),
                                              "clipping"))
    clip_run = sh.Command(which("clipping_profile.py"))
    clip_run(i=in_file, o=out_prefix)

    return out_prefix + ".pdf"

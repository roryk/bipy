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
    clip_plot_file = out_prefix + ".pdf"
    if file_exists(clip_plot_file):
        return clip_plot_file

    clip_run = sh.Command(which("clipping_profile.py"))
    clip_run(i=in_file, o=out_prefix)
    # hack to get around the fact that clipping_profile saves the file in
    # the script execution directory
    sh.mv("clipping_profile.pdf", clip_plot_file)

    return clip_plot_file


def genebody_coverage(in_file, config, out_prefix=None):
    """
    used to check the 5'/3' bias across transcripts
    """
    prefix = "coverage"
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
    coverage_plot_file = out_prefix + ".geneBodyCoverage.pdf"
    if file_exists(coverage_plot_file):
        return coverage_plot_file

    gtf = _get_gtf(config)
    bed = _gtf2bed(gtf)
    coverage_run = sh.Command(which("geneBody_coverage.py"))
    coverage_run(i=in_file, r=bed, o=out_prefix)
    return coverage_plot_file


def junction_annotation(in_file, config, out_prefix=None):
    """
    compile novel/known information about splice junctions
    """
    prefix = "junction"
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
    junction_file = out_prefix + ".splice_events.pdf"
    if file_exists(junction_file):
        return junction_file
    junction_run = sh.Command(which("junction_annotation.py"))
    gtf = _get_gtf(config)
    bed = _gtf2bed(gtf)
    junction_run(i=in_file, o=out_prefix, r=bed)
    return junction_file


def junction_saturation(in_file, config, out_prefix=None):
    """
    check if splicing is deep enough to perform alternative splicing
    analysis
    """
    prefix = "saturation"
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
    saturation_file = out_prefix + ".junctionSsaturation_plot.pdf"
    saturation_run = sh.Command(which("junction_saturation.py"))
    gtf = _get_gtf(config)
    bed = _gtf2bed(gtf)
    saturation_run(i=in_file, o=out_prefix, r=bed)
    return saturation_file


def RPKM_count(in_file, config, out_prefix=None):
    """
    produce RPKM
    """
    prefix = "RPKM_count"
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
    rpkm_count_file = out_prefix + "_read_count.xls"
    gtf = _get_gtf(config)
    bed = _gtf2bed(gtf)
    RPKM_count_run = sh.Command(which("RPKM_count.py"))
    RPKM_count_run(i=in_file, r=bed, o=out_prefix)
    return rpkm_count_file


def RPKM_saturation(in_file, config, out_prefix=None):
    """
    estimate the precision of RPKM calculation by resampling
    """
    prefix = "RPKM_saturation"
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
    rpkm_saturation_file = out_prefix + "saturation.pdf"
    gtf = _get_gtf(config)
    bed = _gtf2bed(gtf)
    RPKM_saturation_run = sh.Command(which("RPKM_saturation.py"))
    RPKM_saturation_run(i=in_file, r=bed, o=out_prefix)
    return rpkm_saturation_file



def _get_out_prefix(in_file, config, out_prefix, prefix):
    if not out_prefix:
        out_prefix = os.path.join(_results_dir(config, prefix),
                                  os.path.basename(in_file))
    return out_prefix

def _get_gtf(config):
    gtf = config["annotation"].get("file", None)
    if not gtf or not file_exists(gtf):
        logger.error("genebody_coverage needs a GTF file passed to it.")
        exit(1)
    return gtf


def _gtf2bed(gtf):
    bed = replace_suffix(gtf, "bed")
    if not file_exists(bed):
        sh.gtf2bed(gtf, _out=bed)
    return bed

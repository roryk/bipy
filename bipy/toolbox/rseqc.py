import os
from bipy.utils import replace_suffix, which, remove_suffix, append_stem
import sh
from bipy.log import logger
from bcbio.utils import file_exists, safe_makedir
from os.path import basename
import glob
import pandas as pd
from math import sqrt
import abc
from mako.template import Template
from bipy.toolbox.reporting import LatexReport


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


def bam2bigwig(in_file, config, out_prefix=None):
    """
    assumes the library preparation was not strand specific for now
    """
    prefix = "bigwig"
    chrom_size_file = config["annotation"].get("chrom_size_file", None)
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
    if not chrom_size_file:
        chrom_size_file = _fetch_chrom_sizes(config)
    wiggle_file = out_prefix + ".wig"

    if not file_exists(wiggle_file):
        bam2wig = sh.Command(which("bam2wig.py"))
        bam2wig(i=in_file, s=chrom_size_file, o=out_prefix)

    bigwig_file = out_prefix + ".bw"

    if file_exists(bigwig_file):
        return bigwig_file

    sh.wigToBigWig(wiggle_file, chrom_size_file, bigwig_file)
    return bigwig_file


def bam_stat(in_file, config, out_prefix=None):
    """
    dump read maping statistics from a SAM or BAM file to out_file
    """
    prefix = "bam_stat"
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
    out_file = out_prefix + ".txt"
    if file_exists(out_file):
        return out_file

    bam_stat_run = sh.Command(which("bam_stat.py"))
    bam_stat_run(i=in_file, _err=out_file)

    return out_file


def clipping_profile(in_file, config, out_prefix=None):
    """
    estimate the clipping profile of the reads
    """
    prefix = "clipping"
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
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
    junction_file = out_prefix + ".splice_junction.pdf"
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
    saturation_file = out_prefix + ".junctionSaturation_plot.pdf"
    if file_exists(saturation_file):
        return saturation_file

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
    if file_exists(rpkm_count_file):
        return rpkm_count_file
    RPKM_count_run = sh.Command(which("RPKM_count.py"))
    RPKM_count_run(i=in_file, r=bed, o=out_prefix)
    return rpkm_count_file


def merge_RPKM(in_dir):
    """
    reads in all RPKM_count files in a directory and combines them into
    one file
    """
    HEADER = ["chrom", "st", "end", "accession" "score",
              "gene_strand", "tag_count", "RPKM"]
    NAME_COL = 3
    HEADER_START = 0
    HEADER_END = 5
    RPKM_COLUMN = 6
    COUNT_COLUMN = 5

    merged_file = os.path.join(in_dir, "RPKM.combined.txt")
    if file_exists(merged_file):
        return merged_file

    RPKM_files = glob.glob(os.path.join(in_dir, "*_read_count.xls"))
    col_names = map(os.path.basename, map(remove_suffix, RPKM_files))
    dataframes = [pd.DataFrame.from_csv(f, sep="\t", index_col=NAME_COL)
                  for f in RPKM_files]
    # make sure the data frames are sorted in the same order
    dataframes = [d.sort_index() for d in dataframes]
    # get the header rows (columns 0-6)
    common = dataframes[0].ix[:, HEADER_START:HEADER_END]
    items = zip(col_names, dataframes)
    for item in items:
        common[item[0] + "_count"] = item[1].ix[:, COUNT_COLUMN]
        common[item[0] + "_RPKM"] = item[1].ix[:, RPKM_COLUMN]

    count_columns = [item[0] + "_count" for item in items]
    rpkm_columns = [item[0] + "_RPKM" for item in items]

    common["count_total"] = common.ix[:, count_columns].sum(axis=1)
    common["RPKM_mean"] = common.ix[:, rpkm_columns].mean(axis=1)
    if len(rpkm_columns) > 1:
        common["RPKM_sd"] = common.ix[:, rpkm_columns].std(axis=1)
        common["RPKM_sem"] = common["RPKM_sd"] / sqrt(len(rpkm_columns))
    else:
        common["RPKM_sd"] = float('NaN')
        common["RPKM_sem"] = float('NaN')

    common.to_csv(merged_file, sep="\t")
    return merged_file


def RPKM_saturation(in_file, config, out_prefix=None):
    """
    estimate the precision of RPKM calculation by resampling
    """
    prefix = "RPKM_saturation"
    out_prefix = _get_out_prefix(in_file, config, out_prefix, prefix)
    rpkm_saturation_file = out_prefix + "RPKM_saturation.pdf"
    gtf = _get_gtf(config)
    bed = _gtf2bed(gtf)

    if file_exists(rpkm_saturation_file):
        return rpkm_saturation_file

    RPKM_saturation_run = sh.Command(which("RPKM_saturation.py"))
    RPKM_saturation_run(i=in_file, r=bed, o=out_prefix)
    return rpkm_saturation_file

def _get_out_dir(in_file, config, out_prefix, prefix):
    if not out_prefix:
        out_dir = os.path.join(_results_dir(config),
                               os.path.basename(in_file),
                               prefix)
        safe_makedir(out_dir)
    return out_dir


def _get_out_prefix(in_file, config, out_prefix, prefix):
    if not out_prefix:
        out_dir = os.path.join(_results_dir(config),
                               os.path.basename(in_file),
                               prefix)
        safe_makedir(out_dir)
        out_prefix = os.path.join(out_dir, prefix)

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


class RseqcParser(object):
    """
    parse a directory full of rseqc results

    """
    DIRS = ["bam_sta", "clipping", "coverage", "junction",
            "saturation", "RPKM_count", "RPKM_saturation"]

    GRAPHS = ((os.path.join("RPKM_saturation", "RPKM_saturation.pdf"),
               "", 1.0),
              (os.path.join("clipping", "clipping.pdf"), "", 1.0),
              (os.path.join("coverage", "coverage.geneBodyCoverage.pdf"),
               "", 1.0),
              (os.path.join("saturation", "junction_saturation.pdf"),
               "", 1.0),
              (os.path.join("junction", "junction.splice_junction.pdf"),
               "", 1.0))

    INFO = (os.path.join("bam_stat", "bam_stat.txt"),
            os.path.join("RPKM_count", "RPKM_count.txt"))

    def __init__(self, base_dir):
        self._dir = base_dir

    def get_rseqc_graphs(self):
        final_graphs = []
        for f, caption, size in self.GRAPHS:
            final_f = os.path.join(self._dir, f)
            if file_exists(final_f):
                final_graphs.append((final_f, caption, size))
        return final_graphs


class RseqcReport(LatexReport):

    CAPTIONS = {"RPKM_saturation.pdf": "",
                "coverage.geneBodyCoverage.pdf": "",
                "clipping.pdf": "",
                "coverage.pdf": "",
                "junction_saturation.pdf": "",
                "splice_events.pdf": ""}

    def __init__(self):
        pass

    def template(self):
        return self._template

    def _add_captions(self, figures):
        new_figures = []
        for figure in figures:
            filename = os.path.basename(figure)
            caption = self.CAPTIONS.get(filename, "")
            new_figures.append((figure[0], caption, figure[2]))
        return new_figures

    def generate_report(self, name, figures=None):
        template = Template(self._template)
        section = template.render(name=name, figures=figures)

    _template = r"""
\subsection*{Rseqc report for ${name}}

% if figures:
    % for i, (figure, caption, size) in enumerate (figures):
    \begin{figure}[htbp]
        \centering
        \includegraphics[width=${size}\linewidth] {${figure}}
        \caption{${caption}}
    \end{figure}
    % endfor
% endif
"""

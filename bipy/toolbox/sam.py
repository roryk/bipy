"""
wrappers around samtools
"""
import sh
from bipy.utils import replace_suffix, append_stem
from bcbio.utils import file_exists


def sam2bam(in_file, out_file=None):
    """ convert a SAM file to a BAM file """
    if out_file is None:
        out_file = replace_suffix(in_file, "bam")

    if file_exists(out_file):
        return out_file

    sh.samtools.view("-Sb", in_file, "-o", out_file)
    return out_file


def bamsort(in_file, out_prefix=None):
    """ sort a BAM file """
    if out_prefix is None:
        out_prefix = replace_suffix(in_file, "sorted")

    out_file = out_prefix + ".bam"

    if file_exists(out_file):
        return out_file

    sh.samtools.sort(in_file, out_prefix)
    return out_file


def bamindex(in_file):
    out_file = in_file + ".bai"

    if file_exists(out_file):
        return out_file

    sh.samtools.index(in_file)
    return out_file

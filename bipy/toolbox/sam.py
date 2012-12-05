"""
wrappers around samtools
"""
import sh
from bipy.utils import replace_suffix, append_stem, is_pair
from bcbio.utils import file_exists
import os
import pysam
from itertools import izip


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


def bamdiff(pair, out_prefix=None):
    """
    takes two coordinate sorted BAM files and outputs only the records unique
    to each one. requires bamUtil to be installed:

    https://github.com/statgen/bamUtil
    """
    if not out_prefix:
        out_prefix = "diff.bam"

    def get_out_file(out_prefix, in_file):
        base, _ = os.path.splitext(out_prefix)
        in_base, in_ext = os.path.splitext(in_file)
        out_file = "_".join([base, "only1", in_base]) + in_ext
        return out_file

    if not is_pair(pair):
        raise ValueError("bamdiff needs to be run on a pair of input BAM "
                         "files.")
    out_files = [get_out_file(out_prefix, x) for x in pair]
    if all(file_exists(out_files)):
        return out_files

    sh.bam.diff(in1=pair[0], in2=pair[1], out=out_prefix)

    return out_files


def disambiguate(pair, cutoff=20):
    """
    takes two coordinate sorted BAM files and returns the reads unique
    or poorly mapping to each file. cutoff specifies the difference in
    mapping quality to call a read as being specific in one file

    """
    # first find the reads unique to file1 and file2
    diff_pair = bamdiff(pair, "diff.bam")

    # find the reads not unique between the two files
    new_pairs = [(pair[0], diff_pair[0]), (pair[1], diff_pair[1])]
    in_both = [bamdiff(x, out_prefix="contested.bam")[0] for x in new_pairs]

    # step through these contested reads and assign them to one file or
    # the other based on passing the quality cutoff
    assigned_files = _assign_by_quality(in_both, cutoff)

    #XXX  merge the files together!
    return assigned_files


def _assign_by_quality(pair, cutoff=20):
    in_handle1 = pysam.Samfile(pair[0], "rb")
    in_handle2 = pysam.Samfile(pair[1], "rb")
    out_files = [x.replace("contested", "assigned") for x in pair]
    out_handle1 = pysam.Samfile(out_files[0], "wb"), template=in_handle1)
    out_handle2 = pysam.Samfile(out_files[1], "wb"), template=in_handle2)
    ambig_files = [x.replace("contested", "ambiguous") for x in pair]
    ambig_handle1 = pysam.Samfile(ambig_files[0], "wb"), template=in_handle1)
    ambig_handle2 = pysam.Samfile(ambig_files[1], "wb"), template=in_handle2)

    for read1, read2 in izip(in_handle1.fetch(), in_handle2.fetch()):
        if (read1.mapq - read2.mapq) > cutoff:
            out_handle1.write(read1)
        elif (read2.mapq - read1.mapq) > cutoff:
            out_handle2.write(read2)
        else:
            # output the reads to the ambiguous files
            ambig_handle1.write(read1)
            ambig_handle2.write(read2)

    return out_files

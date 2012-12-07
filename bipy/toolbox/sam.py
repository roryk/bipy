"""
wrappers around samtools
"""
import sh
from bipy.utils import replace_suffix, append_stem, is_pair
from bcbio.utils import file_exists
from bcbio.distributed.transaction import file_transaction, _flatten_plus_safe
import os
import pysam
from itertools import izip
from bipy.pipeline.stages import AbstractStage


def bam2sam(in_file, out_file=None):
    """ convert a BAM file to a SAM file """
    if is_sam(in_file):
        return in_file

    if out_file is None:
        out_file = replace_suffix(in_file, "sam")

    if file_exists(out_file):
        return out_file

    with file_transaction(out_file) as tmp_out_file:
        sh.samtools.view(in_file, _out=tmp_out_file)

    return out_file


def is_sam_or_bam(in_file):
    return is_bam(in_file) or is_sam(in_file)


def is_bam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext is ".bam":
        return True
    else:
        return False


def is_sam(in_file):
    _, ext = os.path.splitext(in_file)
    if ext is ".sam":
        return True
    else:
        return False


def sam2bam(in_file, out_file=None):
    """ convert a SAM file to a BAM file. if the file is already a
    BAM file, return the BAM file name """

    if is_bam(in_file):
        return in_file

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
    out_dir = os.path.dirname(pair[0])
    if not out_prefix:
        out_prefix = os.path.join(out_dir, "diff.bam")
    else:
        out_prefix = os.path.join(out_dir, out_prefix)

    def get_out_file(out_prefix, in_file, number):
        base, _ = os.path.splitext(out_prefix)
        in_base, in_ext = os.path.splitext(os.path.basename(in_file))
        out_file = "_".join([base, "only" + str(number), in_base]) + in_ext
        return out_file

    if not is_pair(pair):
        raise ValueError("bamdiff needs to be run on a pair of input BAM "
                         "files.")
    out_files = [get_out_file(out_prefix, y, x + 1)
                 for x, y in enumerate(pair)]
    if all(map(file_exists, out_files)):
        return out_files

    sh.bam.diff("--in1", pair[0], "--in2", pair[1], "--out", out_prefix)

    return out_files


class Disambiguate(AbstractStage):
    """
    takes two read name sorted BAM files and returns the reads unique
    or poorly mapping to each file. cutoff specifies the difference in
    mapping quality to call a read as being specific in one file

    """

    def __init__(self, config):
        self.config = config
        self.stage_config = config["stage"]["disambiguate"]
        self.cutoff = self.stage_config.get("cutoff", 20)

    def __call__(self, pair):
        unique_files = [append_stem(x, "unique") for x in pair]
        ambig_files = [append_stem(x, "ambiguous") for x in pair]
        if all(map(os.path.exists, unique_files + ambig_files)):
            return [unique_files, ambig_files]

        handles_0 = self._get_handles(pair[0])
        handles_1 = self._get_handles(pair[1])
        self._process_reads(handles_0, handles_1, None, None)
        [x.close() for x in handles_0]
        [x.close() for x in handles_1]
        return [unique_files, ambig_files]

    def _get_handles(self, in_file):
        assigned_name = append_stem(in_file, "unique")
        ambiguous_name = append_stem(in_file, "ambiguous")

        in_handle = pysam.Samfile(in_file, "rb")
        assigned = pysam.Samfile(assigned_name, "wb", template=in_handle)
        ambiguous = pysam.Samfile(ambiguous_name, "wb", template=in_handle)

        return (in_handle, assigned, ambiguous)

    def _dump_rest(self, handles, read):
        (infile, assigned, ambiguous) = handles
        if read:
            assigned.write(read)
        for read in infile:
            assigned.write(read)

    def _process_reads(self, handles_0, handles_1, read0, read1):
        (in0, assigned0, ambiguous0) = handles_0
        (in1, assigned1, ambiguous1) = handles_1
        if read0 is None:
            try:
                read0 = in0.next()
            except StopIteration:
                self._dump_rest(handles_1, read1)
                return True
        if read1 is None:
            try:
                read1 = in1.next()
            except StopIteration:
                self._dump_rest(handles_0, read0)
                return True

        if read0.qname < read1.qname:
            assigned0.write(read0)
            read0 = None
        elif read1.qname < read0.qname:
            assigned1.write(read1)
            read1 = None
        else:
            score = self._score_read_pair(read0, read1)
            if score == 1:
                assigned0.write(read0)
            elif score == -1:
                assigned1.write(read1)
            else:
                ambiguous0.write(read0)
                ambiguous1.write(read1)
            read0 = None
            read1 = None

        self._process_reads(handles_0, handles_1, read0, read1)

    def _score_read_pair(self, read0, read1):
        if (read0.mapq - read1.mapq) > self.cutoff:
            return 1
        elif (read1.mapq - read0.mapq) > self.cutoff:
            return -1
        else:
            return 0

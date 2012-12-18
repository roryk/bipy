"""
functions to work with fastq files
mate pair fixing script lifted from Peter Cock
# based on this SeqAnswers thread
# http://seqanswers.com/forums/showthread.php?t=6140
# original written by Peter Cock
"""
from Bio import SeqIO
from bipy.utils import append_stem
from bcbio.utils import file_exists
from bipy.log import logger
from bcbio.distributed.transaction import file_transaction
from bipy.pipeline.stages import AbstractStage
import os
from itertools import izip


QUALITY_OFFSETS = {"sanger": 33,
                   "illumina_1.3+": 64,
                   "illumina_1.5+": 66,
                   "illumina_1.8+": 33}

QUALITY_TYPE = {"sanger": "fastq-sanger",
                "solexa": "fastq-solexa",
                "illumina_1.3+": "fastq-illumina",
                "illumina_1.5+": "fastq-illumina",
                "illumina_1.8+": "fastq-sanger"}



def fix_mate_pairs_with_config(fq1, fq2, config):
    if "pair_info" not in config:
        f_suffix = None
        r_suffix = None
    else:
        f_suffix = config["pair_info"].get("forward_read_suffix", None)
        r_suffix = config["pair_info"].get("reverse_read_suffix", None)
    out_files = fix_mate_pairs(fq1, fq2, f_suffix, r_suffix)
    return out_files


def get_read_name_function(suffix):
    if suffix:
        suffix_crop = -len(suffix)

        def read_name_function(read_name):
            """Remove the suffix from a forward read name."""
            assert read_name.endswith(suffix), read_name
            return read_name[:suffix_crop]
    else:
        read_name_function = None
    return read_name_function


def _trim_read(record, bases=8, right_side=True):
    if bases >= len(record):
        record.seq = ""
        record.letter_annotations = {}
        return record

    if right_side:
        return record[:-bases]
    else:
        return record[:bases]


def hard_clip(in_file, bases=8, right_side=True, out_file=None):
    """
    hard clip a fastq file by removing N bases from each read
    bases is the number of bases to clip
    right_side is True to trim from the right side, False to trim from
    the left

    example: hard_clip(fastq_file, bases=4, end="5prime")

    """
    if right_side:
        logger.info("Hard clipping %d bases from the right side of "
                    "reads in %s." % (bases, in_file))
    else:
        logger.info("Hard clipping %d bases from the left side of "
                    "reads in %s." % (bases, in_file))

    quality_type = QUALITY_TYPE[DetectFastqFormat.run(in_file)]
    out_file = append_stem(in_file, "clip")
    if file_exists(out_file):
        return out_file
    in_iterator = SeqIO.parse(in_file, quality_type)

    out_iterator = (_trim_read(record, bases, right_side) for
                    record in in_iterator)
    with file_transaction(out_file) as tmp_out_file:
        with open(tmp_out_file, "w") as out_handle:
            SeqIO.write(out_iterator, out_handle, quality_type)
    return out_file


def filter_single_reads_by_length(in_file, min_length=30):
    """
    removes reads from a fastq file which are below a min_length in bases

    """
    logger.info("Removing reads in %s thare are less than %d bases."
                % (in_file, min_length))
    quality_type = QUALITY_TYPE[DetectFastqFormat.run(in_file)[0]]
    out_file = append_stem(in_file, "fixed")
    if file_exists(out_file):
        return out_file
    in_iterator = SeqIO.parse(in_file, quality_type)
    out_iterator = (record for record in in_iterator if
                    len(record.seq) > min_length)
    with file_transaction(out_file) as tmp_out_file:
        with open(tmp_out_file, "w") as out_handle:
            SeqIO.write(out_iterator, out_handle, quality_type)
    return out_file


def filter_reads_by_length(fq1, fq2, min_length=30):
    """
    removes reads which are empty a pair of fastq files

    """

    logger.info("Removing reads in %s and %s that "
                "are less than %d bases." % (fq1, fq2, min_length))
    # just pick the first one if it can be multiple types
    quality_type = QUALITY_TYPE[DetectFastqFormat.run(fq1)[0]]
    fq1_out = append_stem(fq1, "fixed")
    fq2_out = append_stem(fq2, "fixed")
    fq1_single = append_stem(fq1, "singles")
    fq2_single = append_stem(fq2, "singles")
    if all(map(file_exists, [fq1_out, fq2_out, fq2_single, fq2_single])):
        return [fq1_out, fq2_out]

    fq1_in = SeqIO.parse(fq1, quality_type)
    fq2_in = SeqIO.parse(fq2, quality_type)

    with open(fq1_out, 'w') as fq1_out_handle, open(fq2_out, 'w') as fq2_out_handle, open(fq1_single, 'w') as fq1_single_handle, open(fq2_single, 'w') as fq2_single_handle:
        for fq1_record, fq2_record in izip(fq1_in, fq2_in):
            if len(fq1_record.seq) >= min_length and len(fq2_record.seq) >= min_length:
                fq1_out_handle.write(fq1_record.format(quality_type))
                fq2_out_handle.write(fq2_record.format(quality_type))
            else:
                if len(fq1_record.seq) > min_length:
                    fq1_single_handle.write(fq1_record.format(quality_type))
                if len(fq2_record.seq) > min_length:
                    fq2_single_handle.write(fq2_record.format(quality_type))

    return [fq1_out, fq2_out]


def fix_mate_pairs(fq1, fq2, f_suffix="/1", r_suffix="/2"):
    """
    takes two FASTQ files (fq1 and fq2) of paired end sequencing data
    and filters out reads without a mate pair.
    """
    fq1_out = append_stem(fq1, "fixed")
    fq2_out = append_stem(fq2, "fixed")
    fq1_single = append_stem(fq1, "singles")
    fq2_single = append_stem(fq2, "singles")

    if all(map(file_exists, [fq1_out, fq2_out, fq2_single, fq2_single])):
        return [fq1_out, fq2_out]

    f_dict = SeqIO.index(fq1, "fastq",
                         key_function=get_read_name_function(f_suffix))
    r_dict = SeqIO.index(fq2, "fastq",
                         key_function=get_read_name_function(r_suffix))

    with open(fq1_out, 'w') as fq1_out_handle, open(fq2_out, 'w') as fq2_out_handle, open(fq1_single, 'w') as fq1_single_handle, open(fq2_single, 'w') as fq2_single_handle:
        for key in f_dict:
            if key in r_dict:
                fq1_out_handle.write(f_dict.get_raw(key))
                fq2_out_handle.write(r_dict.get_raw(key))
            else:
                fq1_single_handle.write(f_dict.get_raw(key))
        for key in r_dict:
            if key not in f_dict:
                fq2_single_handle.write(r_dict.get_raw(key))

    return [fq1_out, fq2_out]


class DetectFastqFormat(object):

    _FASTQ_RANGES = {"sanger": [33, 73],
                     "solexa": [59, 104],
                     "illumina_1.3+": [64, 104],
                     "illumina_1.5+": [66, 104],
                     "illumina_1.8+": [33, 74]}

    def __init__(self):
        pass

    @classmethod
    def run(self, in_file, MAX_RECORDS=1000000):
        """
        detects the format of a fastq file
        will return multiple formats if it could be more than one
        """
        kept = list(self._FASTQ_RANGES.keys())
        with open(in_file) as in_handle:
            records_read = 0
            for i, line in enumerate(in_handle):
                # get the quality line
                if records_read >= MAX_RECORDS:
                    break
                if i % 4 is 3:
                    records_read += 1
                    for c in line.strip():
                        formats = kept
                        if len(formats) == 1:
                            return formats
                        for form in formats:
                            if (self._FASTQ_RANGES[form][0] > ord(c) or
                                self._FASTQ_RANGES[form][1] < ord(c)):
                                kept.remove(form)

        return formats

    def __call__(self, in_file, MAX_RECORDS=1000000):
        logger.info("Detecting format of %s" % (in_file))
        quality = self.run(in_file, MAX_RECORDS)
        logger.info("Detected quality format of %s in %s." % (quality, in_file))
        return self.run(in_file, MAX_RECORDS)


class FastqGroomer(AbstractStage):
    """
    Grooms a FASTQ file from its input format to sanger format

    """
    stage = "groom"

    def out_file(self, in_file):
        """
        returns the expected output file name from the in_file

        example: "control_1.fastq" -> "control_1.groom.fastq"

        """
        results_dir = self.config["dirs"].get("results", "results")
        stage_dir = os.path.join(results_dir, self.stage)
        out_file = append_stem(os.path.basename(in_file), "groom")
        return os.path.join(stage_dir, out_file)

    def __init__(self, config):
        self.config = config

    def _run(self, in_file):
        pass

    def _detect_format(self, in_file):
        quality_format = DetectFastqFormat(in_file)
        return QUALITY_TYPE[quality_format]

    def __call__(self, in_file):
        out_file = self.out_file(in_file)
        if file_exists(out_file):
            return out_file
        raise NotImplementedError


class HardClipper(AbstractStage):
    """
    Hard clips records in a FASTQ file, clipping off bases from either
    end

    """
    stage = "hard_clip"

    def __init__(self, config):
        self.config = config
        self.bases = config[self.stage].get("bases", 8)
        self.right_side = config[self.stage].get("right_side", True)

    def out_file(self, in_file):
        return append_stem(in_file, "clip")

    def __call__(self, in_file):
        out_file = self.out_file(in_file)
        if file_exists(out_file):
            return out_file
        hard_clip(in_file, self.bases, self.right_side, out_file)
        return out_file

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

QUALITY_OFFSETS = {"sanger": 33,
                   "illumina_1.3+": 64,
                   "illumina_1.5+": 66,
                   "illumina_1.8+": 33}


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

def _format_fastq_record(record, quality_offset):
    x = "@" + record.description + "\n"
    x += str(record.seq) + "\n"
    x += "+" + record.description + "\n"
    x += "".join(map(chr, [q + quality_offset for q in
                           record.letter_annotations['phred_quality']]))
    x += "\n"
    return x


def filter_reads_by_length(fq1, fq2, min_length=20):
    """
    removes reads which are empty a pair of fastq files

    """

    logger.info("Removing empty reads on %s and %s." % (fq1, fq2))
    # just pick the first one if it can be multiple types
    quality_offset = QUALITY_OFFSETS[DetectFastqFormat.run(fq1)[0]]
    fq1_out = append_stem(fq1, "fixed")
    fq2_out = append_stem(fq2, "fixed")
    fq1_single = append_stem(fq1, "singles")
    fq2_single = append_stem(fq2, "singles")
    if all(map(file_exists, [fq1_out, fq2_out, fq2_single, fq2_single])):
        return [fq1_out, fq2_out]

    fq1_in = SeqIO.parse(fq1, "fastq")
    fq2_handle = open(fq2)
    fq2_in = SeqIO.parse(fq2, "fastq")

    with open(fq1_out, 'w') as fq1_out_handle, open(fq2_out, 'w') as fq2_out_handle, open(fq1_single, 'w') as fq1_single_handle, open(fq2_single, 'w') as fq2_single_handle:
        for fq1_record, fq2_record in zip(fq1_in, fq2_in):
            if len(fq1_record.seq) >= min_length and len(fq2_record.seq) >= min_length:
                fq1_out_handle.write(_format_fastq_record(fq1_record,
                                                          quality_offset))
                fq2_out_handle.write(_format_fastq_record(fq2_record,
                                                          quality_offset))
            else:
                if len(fq1_record.seq) > 0:
                    fq1_single_handle.write(_format_fastq_record(fq1_record,
                                                                 quality_offset))
                if len(fq2_record.seq) > 0:
                    fq2_single_handle.write(_format_fastq_record(fq2_record,
                                                                 quality_offset))

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
        return self.run(in_file, MAX_RECORDS)

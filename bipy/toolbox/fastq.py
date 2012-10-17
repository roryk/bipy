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


def fix_mate_pairs_with_config(fq1, fq2, config):
    if "pair_info" not in config:
        f_suffix = None
        r_suffix = None
    f_suffix = config["pair_info"].get("forward_read_suffix", None)
    r_suffix = config["pair_info"].get("reverse_read_suffix", None)
    fix_mate_pairs(fq1, fq2, f_suffix, r_suffix)


def get_read_name_function(suffix):
    if suffix:
        suffix_crop = -len(suffix)

        def f_name(read_name):
            """Remove the suffix from a forward read name."""
            assert read_name.endswith(suffix), read_name
            return read_name[:suffix_crop]
    else:
        read_name = None
    return read_name


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

"""This module provides an interface to cutadapt with a set of commonly
used adapters for trimming
"""
from bipy.utils import flatten_options, append_stem, flatten, which
from operator import itemgetter
import subprocess
import os
from bcbio.utils import safe_makedir, file_exists
import difflib
import sh


# adapter sequences for various commonly used systems
ADAPTERS = {}
ADAPTERS["illumina"] = [
    ["ACACTCTTTCCCTACACGACGCTCTTCCGATCT", "-a", "ill_pe_adapter1"],
    ["TGTGAGAAAGGGATGTGCTGCGAGAAGGCTAG", "-a", "ill_pe_adapter1_rc"],
    ["GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG", "-a", "ill_pe_adapter2"],
    ["TCTAGCCTTCTCGCCAAGTCGTCCTTACGGCTC", "-a", "ill_pe_adapter2_rc"]]
ADAPTERS["nextera"] = [
    ["AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAG", "-a",
     "nex_pe_adapter1"],
    ["CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATT", "-a",
     "nex_pe_adapter1_rc"],
    ["CAAGCAGAAGACGGCATACGAGATCGGTCTGCCTTGCCAGCCCGCTCAG",
     "-a", "nex_pe_adapter2_nobc"],
    ["CTGAGCGGGCTGGCAAGGCAGACCGATCTCGTATGCCGTCTTCTGCTTG",
     "-a", "nex_pe_adapter2_nobc_rc"],
    ["CTGATGGCGCGAGGGAGGCGTGTAGATCTCGGTGGTCGCCGTATCATTCTGTCTCTTATACACATCT",
     "-a", "nex_transposon_pe_adapter1_rc"],
    ["AGATGTGTATAAGAGACAGAATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAG",
     "-a", "nex_transposon_pe_adapter1"],
    ["AGATGTGTATAAGAGACAGCAAGCAGAAGACGGCATACGAGATCGGTCTGCCTTGCCAGCCCGCTCAG",
     "-a", "nex_tranposon_pe_adapter2"]]
ADAPTERS["polya"] = [
    ["AAAAAAAAAAAAAAAAAAAAAAAAAAA", "-a", "polyA tail"],
    ["TTTTTTTTTTTTTTTTTTTTTTTTTTT", "-a", "polyT tail"]]
ADAPTERS["iontorrent"] = [
    ["CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT", "-a",
     "ion_5_prime_adapter"],
    ["CTGAGTCGGAGACACGCAGGGATGAGATGG", "-a", "3_prime_adapter"],
    ["ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGG", "-a",
     "5_prime_adapter_rc"],
    ["CCATCTCATCCCTGCGTGTCTCCGACTCAG", "-a", "3_prime_adapter_rc"]]

TRUSEQ_BARCODES = {"ATCACG": 1, "AGTCAA": 13, "ACTGAT": 25, "CGGAAT": 37,
                   "CGATGT": 2, "AGTTCC": 14, "ATGAGC": 26, "CTAGCT": 38,
                   "TTAGGC": 3, "ATGTCA": 15, "ATTCCT": 27, "CTATAC": 39,
                   "TGACCA": 4, "CCGTCC": 16, "CAAAAG": 28, "CTCAGA": 40,
                   "ACAGTG": 5, "GTAGAG": 17, "CAACTA": 29, "GACGAC": 41,
                   "GCCAAT": 6, "GTCCGC": 18, "CACCGG": 30, "TAATCG": 42,
                   "CAGATC": 7, "GTGAAA": 19, "CACGAT": 31, "TACAGC": 43,
                   "ACTTGA": 8, "GTGGCC": 20, "CACTCA": 32, "TATAAT": 44,
                   "GATCAG": 9, "GTTTCG": 21, "CAGGCG": 33, "TCATTC": 45,
                   "TAGCTT": 10, "CGTACG": 22, "CATGGC": 34, "TCCCGA": 46,
                   "GGCTAC": 11, "GAGTGG": 23, "CATTTT": 35, "TCGAAG": 47,
                   "CTTGTA": 12, "GGTAGC": 24, "CCAACA": 36, "TCGGCA": 48}

VALID_TRUSEQ_RNASEQ = {k: v for (k, v) in TRUSEQ_BARCODES.items() if v < 13}


def truseq_barcode_lookup(barcode, small=False):
    """
    looks up a truseq adapter sequence by inserting the barcode in the
    correct sequence. throws an exception if the barcode does not match
    known barcodes

    """
    prefix = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
    suffix = "ATCTCGTATGCCGTCTTCTGCTTG"
    if small:
        raise NotImplementedError("Small RNA barcodes not implemented. Need "
                                  "to check to make sure the prefix and "
                                  "suffix sequences are the same as the "
                                  "RNA-seq barcodes.")
    if barcode not in VALID_TRUSEQ_RNASEQ:
        raise ValueError("Barcode not found in TruSeq barcodes. Might need "
                         "to implement v1 and v2 versions.")


    return prefix + barcode + suffix


def _get_adapter(adapter):
    return [adapter[1], adapter[0]]

def _get_platform_adapters(platform):
    platform_adapters = ADAPTERS.get(platform, [])
    adapters = map(_get_adapter, platform_adapters)
    return adapters

def _parse(config):
    # handle the adapters, defaulting to illumina and a poly-a trimmer
    # if none are provided
    adapters = []
    adapters += flatten(map(_get_adapter,
                            config.get("adapters", [])))
    # add built in platform if available
    platform = config.get("platform", None)
    if platform:
        adapters += flatten(map(_get_platform_adapters,
                                [p for p in platform if p in ADAPTERS]))
    # default to illumina and poly A
    if not adapters:
        adapters += flatten(map(_get_platform_adapters,
                        [p for p in ["illumina", "polya"]]))

    arguments = []
    arguments += adapters
    # grab everything else
    arguments += config.get("options", [])
    return map(str, list(flatten(arguments)))


def run(in_file, stage_config, config):
    arguments = [stage_config["program"]]
    arguments += _parse(stage_config)
    results_dir = config["dir"].get("results", None)
    if results_dir:
        out_dir = os.path.join(results_dir, "cutadapt")
        safe_makedir(out_dir)
        out_file = os.path.join(out_dir,
                                os.path.basename(append_stem(in_file,
                                                             "trimmed")))
    else:
        out_file = append_stem(in_file, "trimmed")

    if file_exists(out_file):
        return out_file

    arguments.extend(["--output", out_file, in_file])
    subprocess.check_call(arguments)
    return out_file


def fix_paired_files(first, second, out_file=None):
    if out_file is None:
        common = _common_prefix(first, second) + "fixed"
    mergeFastq = sh.Command(which("mergeShuffledFastqSeqs.pl"))
    #cmd = ["mergeShuffledFastqSeqs.pl", "-f1", first, "-f2", second,
    #       "-o", common, "-t", "-r", "^@(\S+)/[1|2]$"]
           #subprocess.check_call(cmd)
    mergeFastq(f1=first, f2=second, o=common, t=True, r=r'^@(\S+)/[1|2]$')
    # XXX these might not be right.
    return (common + "_1.fastq", common + "_2.fastq")


def _common_prefix(first, second):
    for i, (x, y) in enumerate(zip(first, second)):
        if x != y:
            break
    return first[:i]

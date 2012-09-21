"""This module provides an interface to cutadapt with a set of commonly
used adapters for trimming
"""
from bipy.utils import flatten_options, append_stem, flatten
from operator import itemgetter
import subprocess
import os
from bcbio.utils import safe_makedir, file_exists

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

"""
runs jellyfish with various options
"""
from bipy.utils import flatten, replace_suffix, build_results_dir, append_stem
import os
import subprocess


def _build_output_prefix(input_file, jellyfish_config, config):
    out_dir = build_results_dir(jellyfish_config, config)
    out_prefix = os.path.join(out_dir, replace_suffix(input_file,
                                                      "count"))
    #out_prefix = "_".join([jellyfish_config["name"],
    #                       remove_suffix(input_file)])
    return out_prefix


def _build_command(input_file, out_prefix, jellyfish_config):
    cmd = ["jellyfish", jellyfish_config["task"],
           jellyfish_config["options"]]
    cmd += ["-o", out_prefix, input_file]
    return list(flatten(cmd))


def _build_merge_command(out_prefix, out_file):
    cmd = ["jellyfish", "merge", "-o", out_file, out_prefix + "_*"]
    return(cmd)


def run(input_file, jellyfish_config, config):
    # run the jellyfish counting, this produces a set of files identified
    # by out_prefix
    out_prefix = _build_output_prefix(input_file, jellyfish_config, config)
    cmd = _build_command(input_file, out_prefix, config)
    subprocess.check_call(cmd)

    # combine the output files into one merged file and return that
    out_file = append_stem(out_prefix, "combined")
    merge_cmd = _build_merge_command(out_prefix, out_file)
    subprocess.check_call(merge_cmd)
    # find all of the output files and merge them into one file
    return out_file

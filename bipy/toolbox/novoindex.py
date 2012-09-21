""" handle running novoindex, making a new genome """
import os
from bipy.utils import replace_suffix, flatten
from bcbio.utils import safe_makedir
import subprocess


def _build_output_file(input_file, config):
    safe_makedir(config["dir"]["ref"])
    return os.path.join(config["dir"]["ref"],
                        os.path.basename(replace_suffix(input_file, "nix")))


def _build_command(input_file, novoindex_config, output_file):
    options = novoindex_config["options"].items()
    cmd = map(str, (flatten(["novoindex", options, output_file,
                             input_file])))
    return cmd


def run(input_file, novoindex_config, config):
    output_file = _build_output_file(input_file, config)

    # skip if already made
    if os.path.exists(output_file):
        return output_file

    cmd = _build_command(input_file, novoindex_config, output_file)

    subprocess.check_call(cmd)
    return output_file

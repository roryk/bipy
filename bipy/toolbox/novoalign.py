""" handling running novoalign """
from bipy.utils import (build_results_dir, flatten, which, replace_suffix,
                         flatten_options)
from bcbio.utils import safe_makedir
import os
import subprocess
import logging

logger = logging.getLogger(__name__)


def _build_output_file(input_file, novoalign_config, config):
    outdir = build_results_dir(novoalign_config, config)
    safe_makedir(outdir)
    return os.path.join(outdir,
                        os.path.basename(replace_suffix(input_file, "sam")))


def _build_command(input_file, ref, novoalign_config):
    cmd = [which("novoalign"), flatten_options(novoalign_config),
           "-o", "SAM", "-d", ref, "-f", input_file]
    return list(map(str, flatten(cmd)))


def run(input_file, ref, novoalign_config, config):
    output_file = _build_output_file(input_file, novoalign_config, config)
    logger.info("Running novoalign on %s "
                "and outputting to %s." % (input_file, output_file))
    # skip if we already did this
    if os.path.exists(output_file):
        logger.info("Skipping %s, already complete." % (input_file))
        return output_file

    cmd = _build_command(input_file, ref, novoalign_config)

    # capture stdout from novoalign and redirect to the output file
    with open(output_file, "w") as out_handle:
        subprocess.check_call(cmd, stdout=out_handle)

    logger.info("Novoalign complete. Output in %s " % (output_file))
    return output_file

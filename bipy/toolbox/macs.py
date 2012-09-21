"""
launcher to run macs
"""

"""
 macs14 --treatment (treatment file, bam file from alignment)
--control (control file from scrambles (if present)
--format (BAM or SAM)
--name experiment_name
--gsize effective genome size (hs for human, mm for mouse, ce, dm or an int)
--pvalue cutoff (1e-5 is the default)
--mfold rangeof high confidence ratio to background (default is 10,30)
"""

from bipy.utils import flatten, remove_suffix, replace_suffix
import subprocess
from bipy.log import logger
import os
from bcbio.utils import safe_makedir


def _build_command(input_file, options, control_file=None, out_dir=None):
    name = remove_suffix(os.path.basename(input_file))
    if out_dir:
        name = os.path.join(out_dir, name)

    options = ["=".join(map(str, x)) for x in options]

    cmd = ["macs14", "--treatment=" + input_file, flatten(options),
           "--name=" + name]
    if control_file:
        cmd += ["--control=" + control_file]

    return map(str, flatten(cmd))


def run(input_file, options, control_file=None, out_dir=None):
    out_files = (remove_suffix(input_file) + "_peaks.bed",
                 remove_suffix(input_file) + "_summits.bed")
    if out_dir:
        out_files = [os.path.join(out_dir, os.path.basename(x)) for
                     x in out_files]
    cmd = _build_command(input_file, options, control_file, out_dir)
    subprocess.check_call(cmd)
    return out_files


def run_with_config(input_file, config, control_file=None, stage=None):

    if stage is None:
        stage = "macs"

    if stage not in config["stage"]:
        logger.info("Cannot find the the stage %s in the config." % (stage))

    stage_config = config["stage"][stage]
    options = stage_config.get("options", [])
    out_dir = os.path.join(config["dir"].get("results", None), stage)
    safe_makedir(out_dir)
    out_files = run(input_file, options, control_file, out_dir)
    print out_files
    return out_files

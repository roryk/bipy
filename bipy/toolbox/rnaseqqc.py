"""
rna_seqc wrapper
java -jar RNA-SeQC.jar -o output_dir -r reference.fa -t reference.gtf
-s "sampleid|bam file|notes"

notes = "rna_seqc"
bam_file = in_file
sample_id = os.path.basename(remove_suffix(in_file))

"""
from bipy.log import logger
from bcbio.utils import file_exists, safe_makedir
from bipy.utils import get_stem, flatten
import os
import subprocess


def _validate_config(in_file, stage_config, config):
    """ validates that a set of assumptions about the config file
    needed to run the program are true """
    if "ref" not in config:
        logger.error("ref: must appear in the config file")
        exit(1)
    if not file_exists(config["ref"] + ".fa"):
        logger.error("%s not found, aborting." % (config["ref_fasta"]))
    if not file_exists(in_file):
        logger.error("%s not found, aborting." % (in_file))
    if not file_exists(config["gtf"]):
        logger.error("%s not found, aborting." % (config["gtf"]))
    if not file_exists(stage_config["program"]):
        logger.error("%s not found, aborting." % (stage_config["program"]))


def _build_command(in_file, stage_config, config):
    cmd = ["java", "-jar", stage_config["program"]]
    out_dir = os.path.join(config["dir"]["results"],
                           stage_config.get("name", "rna_seqc"),
                           get_stem(in_file))
    safe_makedir(out_dir)
    cmd += ["-o", out_dir]
    sample = "|".join([get_stem(in_file), in_file, "rna_seqc"])
    cmd += ["-s", sample]
    cmd += ["-r", config["ref_fasta"]]
    cmd += ["-t", config["gtf"]]
    cmd += [stage_config.get("options", [])]
    return list(flatten(cmd))


def run_with_config(in_file, stage_config, config):
    _validate_config(in_file, stage_config, config)
    cmd = _build_command(in_file, stage_config, config)
    subprocess.check_call(cmd)

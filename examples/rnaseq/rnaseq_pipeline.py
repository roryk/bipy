"""
example script for running a RNA-seq analysis

python rnaseq_pipeline.py rnaseq_pipeline.yaml

you will have to write a couple of functions to group the input
data in useful ways

"""
from bipy.cluster import start_cluster, stop_cluster
import sys
import yaml
from bipy.log import setup_logging, logger
from bcbio.utils import safe_makedir, file_exists
import os
from bipy.utils import (combine_pairs, flatten, append_stem,
                        prepare_ref_file, replace_suffix)
from bipy.toolbox import (htseq_count, deseq, annotate, rseqc, sam)
from bcbio.broad import BroadRunner, picardrun
from bipy.toolbox.trim import Cutadapt
from bipy.toolbox.fastqc import FastQC
from bipy.toolbox.tophat import Tophat

import glob
from itertools import product, repeat
import sh


def find_files(in_dir):
    """
    returns a list of the sequence files in a directory recursively

    """

    FASTQ_EXTENSIONS = [".fq", ".fastq"]
    files = [sh.find(in_dir, "-name", "*" + x) for x in FASTQ_EXTENSIONS]
    return files


def make_test(in_file, lines=1000000):
    """
    take a small subset of the input files for testing. only makes sense for
    text files where lines gives an appopriate number of records, for example,
    FASTQ files should be a multiple of 4.

    """
    out_dir = os.path.join(os.path.dirname(in_file), "test")
    safe_makedir(out_dir)
    out_file = os.path.join(out_dir,
                            append_stem(os.path.basename(in_file), "test"))
    sh.head("-" + str(lines), in_file, _out=out_file)
    return out_file


def _get_stage_config(config, stage):
    return config["stage"][stage]


def _get_program(config, stage):
    return config["stage"][stage]["program"]


def _emit_stage_message(stage, curr_files):
    logger.info("Running %s on %s" % (stage, curr_files))


def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)

    # make the needed directories
    map(safe_makedir, config["dir"].values())

    # specific for project
    input_dir = config["input_dir"]
    input_files = find_files(input_dir)

    results_dir = config["dir"].get("results", "results")
    if config.get("test_pipeline", False):
        curr_files = map(make_test, input_files)
        results_dir = os.path.join(results_dir, "test_pipeline")
        config["dir"]["results"] = results_dir
    else:
        curr_files = input_files

    for stage in config["run"]:
        if stage == "fastqc":
            logger.info("Running fastqc on %s." % (curr_files))
            stage_runner = FastQC(config)
            view.map(stage_runner, curr_files, block=False)

        if stage == "cutadapt":
            curr_files = combine_pairs(curr_files)
            logger.info("Running cutadapt on %s." % (curr_files))
            stage_runner = Cutadapt(config)
            curr_files = view.map(stage_runner, curr_files)

        if stage == "tophat":
            logger.info("Running tophat on %s." % (curr_files))
            stage_runner = Tophat(config)
            tophat_outputs = view.map(stage_runner, curr_files)
            bamfiles = view.map(sam.sam2bam, tophat_outputs)
            bamsort = view.map(sam.bamsort, bamfiles)
            view.map(sam.bamindex, bamsort)
            final_bamfiles = bamsort
            curr_files = tophat_outputs

        if stage == "htseq-count":
            logger.info("Running htseq-count on %s." % (curr_files))
            htseq_args = zip(*product(curr_files, [config], [stage]))
            htseq_outputs = view.map(htseq_count.run_with_config,
                                     *htseq_args)

        if stage == "coverage":
            logger.info("Calculating RNASeq metrics on %s." % (curr_files))
            nrun = len(curr_files)
            ref = prepare_ref_file(config["stage"][stage]["ref"], config)
            ribo = config["stage"][stage]["ribo"]
            picard = BroadRunner(config["program"]["picard"])
            out_dir = os.path.join(results_dir, stage)
            safe_makedir(out_dir)
            out_files = [replace_suffix(os.path.basename(x),
                                        "metrics") for x in curr_files]
            out_files = [os.path.join(out_dir, x) for x in out_files]
            out_files = view.map(picardrun.picard_rnaseq_metrics,
                                 [picard] * nrun,
                                 curr_files,
                                 [ref] * nrun,
                                 [ribo] * nrun,
                                 out_files)

        if stage == "rseqc":
            _emit_stage_message(stage, curr_files)
            rseq_args = zip(*product(curr_files, [config]))
            view.map(rseqc.bam_stat, *rseq_args)
            view.map(rseqc.genebody_coverage, *rseq_args)
            view.map(rseqc.junction_annotation, *rseq_args)
            view.map(rseqc.junction_saturation, *rseq_args)
            RPKM_args = zip(*product(final_bamfiles, [config]))
            RPKM_count_out = view.map(rseqc.RPKM_count, *RPKM_args)
            RPKM_count_fixed = view.map(rseqc.fix_RPKM_count_file,
                                        RPKM_count_out)
            """
                            annotate_args = zip(*product(RPKM_count_fixed,
                                         ["gene_id"],
                                         ["ensembl_gene_id"],
                                         ["human"]))
            view.map(annotate.annotate_table_with_biomart,
                     *annotate_args)
                     """
            view.map(rseqc.RPKM_saturation, *rseq_args)
            curr_files = tophat_outputs

    # end gracefully
    stop_cluster()


if __name__ == "__main__":
    # read in the config file and perform initial setup
    main_config_file = sys.argv[1]
    with open(main_config_file) as config_in_handle:
        startup_config = yaml.load(config_in_handle)
    setup_logging(startup_config)
    start_cluster(startup_config)
    from bipy.cluster import view

    main(main_config_file)

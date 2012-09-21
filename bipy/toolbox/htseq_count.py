""" wrapper for using the htseq-count script """
import os
from bipy.utils import replace_suffix, flatten, prepare_ref_file
import subprocess
from bcbio.utils import safe_makedir
import pandas as pd
from bipy.log import logger
import rpy2.robjects as robjects
import HTSeq
from bipy import gtf


def _load_htseq_count_file(filename):
    return pd.read_csv(filename, sep="\t", index_col=0, header=None)


def _make_length_dict(gtf_file):
    HTSeq.GFF_Reader(gtf_file)
    lengths = {}
    for feature in gtf_file:
        lengths[feature.name] = 0


def calculate_rpkm(count_file, gtf_file):
    """ calculates RPKM for each column in the count_file using
    transcript lengths calculated from the gtf_file and the sum
    of the counts in the count_file. uses the average length of the
    transcripts for a gene as the gene length """
    # calculate the lengths of the genes in gtf_file
    gtflines = gtf.GTFtoDict(gtf_file)
    genes = gtf.aggregateFeaturesByGene(gtflines)
    #genes = gtf.mergeOverlappedExons(genes)
    lengths = {}
    for gene, features in genes.items():
        lengths[gene] = 0
        transcripts = set()
        for feature in features:
            if feature["feature"] == "exon":
                transcripts.add(feature["transcript_id"])
                lengths[gene] += abs(int(feature["start"]) -
                                     int(feature["end"]))
        lengths[gene] = lengths[gene] / float(len(transcripts))

    counts_df = pd.read_table(count_file, header=0, index_col=0)
    length_list = [lengths.get(x, 1000) / 1000.0 for x in counts_df.index]
    length_normalized = (counts_df.transpose() / length_list).transpose()
    rpkm = length_normalized.apply(lambda x: x /
                                   (float(x.sum()) / 1000000))
    return rpkm


def combine_counts(in_files, column_names=None, out_file=None):
    if column_names is None:
        column_names = in_files
    if out_file is None:
        out_file = os.path.join(os.path.dirname(in_files[0]),
                                                "combined.counts")
    r = robjects.r
    r.assign('in_files', robjects.StrVector(in_files))
    r.assign('out_file', out_file)
    r.assign('col_names', robjects.StrVector(column_names))
    r('''
    df = read.table(in_files[1], sep="\t", row.names=1)
    for (file in in_files[-1]) {
    df = cbind(df, read.table(file, sep="\t", row.names=1))
    }
    colnames(df) = col_names
    write.table(df, file=out_file, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
    ''')
    return out_file

def _get_outfilename(input_file):
    out_file = replace_suffix(os.path.basename(input_file), "counts")
    return out_file


def run(input_file, gtf_file, options=None, out_file=None):
    if options is None:
        options = []
    if out_file is None:
        out_file = _get_outfilename(input_file)

    safe_makedir(os.path.dirname(out_file))

    if os.path.exists(out_file):
        return out_file

    cmd = map(str, flatten(["htseq-count", options, input_file, gtf_file]))
    with open(out_file, "w") as out_handle:
        subprocess.check_call(cmd, stdout=out_handle)

    return out_file


def run_with_config(input_file, config, stage, out_file=None):
    stage_config = config["stage"][stage]
    options = stage_config.get("options", [])

    if out_file is None:
        out_dir = os.path.join(config["dir"].get("results", None), stage)
        out_file = os.path.join(out_dir, _get_outfilename(input_file))

    safe_makedir(out_dir)
    if "annotation" not in config:
        logger.error("annotation must appear in the config file, see example "
                     "configuration files.")
        exit(1)
    ref = prepare_ref_file(config["annotation"], config)
    out_file = run(input_file, ref, options, out_file)
    return out_file

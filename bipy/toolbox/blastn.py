""" This code is snagged from Brad Chapman almost wholesale """
import bcbio.utils as utils
from bipy.utils import (replace_suffix, append_stem, build_results_dir)
import os
import subprocess
import csv
from bipy.log import logger


HEADER_FIELDS = ('qseqid sseqid pident length mismatch gapopen qstart qend '
                 'sstart send evalue bitscore qlen slen')


def get_id_of_hits(outfile):
    """ returns a set of all of the ids that had hits in an outfile """
    s = set()
    with open(outfile) as in_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        for line in reader:
            linedict = dict(zip(HEADER_FIELDS.split(" "), line))
            s.add(linedict["qseqid"])
    return s


def is_long_enough(linedict, cutoff):
    """ a nicer version of filter_results_by_length, not tested yet though """
    prefixes = ("s", "q")

    def percent_match(prefix):
        length = abs(float(linedict[prefix + "start"]) - float(
            linedict[prefix + "end"]))
        if length / float(linedict[prefix + "len"]) > cutoff:
            return True
        else:
            return False

    return all(map(percent_match, prefixes))


def filter_results_by_length(filename, cutoff):
    """ filters the tsv results by the metric that both the overlap
    of the query sequence and the subject sequence must both be
    > cutoff of their length. This might be a little too restrictive though
    """
    def query_match(linedict):
        length = abs(float(linedict["qstart"]) - float(linedict["qend"]))
        if length / float(linedict["qlen"]) > (cutoff / float(100)):
            return True
        else:
            return False

    def subject_match(linedict):
        length = abs(float(linedict["sstart"]) - float(linedict["send"]))
        if length / float(linedict["slen"]) > (cutoff / float(100)):
            return True
        else:
            return False

    out_fname = append_stem(filename, str(cutoff) + "_filt")
    # skip if it already exists
    if os.path.exists(out_fname):
        return out_fname

    with open(filename) as in_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        with open(out_fname, "w") as out_handle:
            writer = csv.writer(out_handle, delimiter="\t")
            writer.writerow(HEADER_FIELDS.split(" "))
            for line in reader:
                linedict = dict(zip(HEADER_FIELDS.split(" "), line))
                if query_match(linedict) & subject_match(linedict):
                    writer.writerow(line)
    return out_fname


def run(in_file, ref, blastn_config, config):
    logger.info("Preparing the reference file for %s." % (ref.get("name")))
    ref_file = prepare_ref_file(ref, config)
    logger.info("Preparing the blast database for %s." % (ref.get("name")))
    blast_db = prepare_blast_db(ref_file, "nucl")
    logger.info("Blasting %s against %s." % (in_file, ref.get("name")))

    results_dir = build_results_dir(blastn_config, config)
    utils.safe_makedir(results_dir)

    out_file = os.path.join(results_dir,
                            replace_suffix(os.path.basename(in_file),
                                           ref.get("name") + "hits.tsv"))
    tmp_out = out_file + ".tmp"

    blast_results = blast_search(in_file, blast_db, tmp_out)
    #logger.info("Filtering results for at least %f percent of the "
    #            "sequences covered." %(0.5*100))
    #filtered_results = filter_results_by_length(blast_results, 0.5)
    #logger.info("Filtered output file here: %s" %(filtered_results))
    with open(blast_results) as in_handle:
        reader = csv.reader(in_handle, delimiter="\t")
        with open(out_file, "w") as out_handle:
            writer = csv.writer(out_handle, delimiter="\t")
            writer.writerow(HEADER_FIELDS.split(" "))
            for line in reader:
                writer.writerow(line)

    return out_file


# de-parallelized for now
def blast_search(in_file, blast_db, out_file):
    if not os.path.exists(out_file):
        _do_blast(in_file, blast_db, out_file)
    return out_file


def prepare_blast_db(db_file, dbtype):
    exts = {"prot": "pin", "nucl": "nin"}
    base, _ = os.path.splitext(db_file)
    if not os.path.exists("%s.%s" % (base, exts[dbtype])):
        cl = ["makeblastdb", "-in", db_file, "-out", base, "-dbtype", dbtype]
        subprocess.check_call(cl)
    return base


def _do_blast(in_file, blast_db, out_file):
    cl = ["blastn", "-query", in_file, "-outfmt", '6 ' + HEADER_FIELDS, "-out",
          out_file, "-db", blast_db, "-num_alignments", "1",
          "-num_descriptions", "1", "-evalue", "0.1"]
    subprocess.check_call(cl)


def prepare_ref_file(ref, config):
    """Get a reference file, either by URL or locally.
    """
    url = ref.get("url", None)
    if url:
        ref_file = _download_ref(url, config["dir"]["ref"])
    else:
        ref_file = ref.get("file", None)
    assert ref_file is not None and os.path.exists(ref_file), ref_file
    return ref_file


def _download_ref(url, ref_dir):
    dl_file = os.path.basename(url)
    ref_file = None
    for supported_ext, extract_cmd in [(".gz", "gunzip")]:
        if dl_file.endswith(supported_ext):
            ref_file = os.path.join(ref_dir, dl_file[:-len(supported_ext)])
            break
    assert ref_file is not None, url
    if not os.path.exists(ref_file):
        with utils.chdir(ref_dir):
            cl = ["wget", url]
            subprocess.check_call(cl)
            cl = [extract_cmd, dl_file]
            subprocess.check_call(cl)
    return ref_file

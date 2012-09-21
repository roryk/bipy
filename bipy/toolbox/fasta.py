""" tools for handling fasta files """
from Bio import SeqIO
from contextlib import nested
import os


def count_seqio(in_file, predicate, kind):
    """ counts the records in a file that pass a predicate """
    with open(in_file) as in_handle:
        return reduce(lambda x, y: x + 1,
                      filter(predicate, SeqIO.parse(in_handle, kind)), 0)


def passed_seqio(in_file, predicate, kind):
    """ returns the number of records that pass a predicate
     and the total  number of records as a tuple (passed, total) """
    with open(in_file) as in_handle:
        passed = 0
        total = 0
        for record in SeqIO.parse(in_handle, kind):
            if predicate(record):
                passed += 1
            total += 1
    return (passed, total)


def apply_seqio(in_file, f, kind):
    """ apply a function f to every record in in_file """
    with open(in_file) as in_handle:
        return map(f, SeqIO.parse(in_handle, kind))


def filter_fasta(in_file, predicate, out_file):
    """ filters a fasta by a predicate """
    filter_seqio(in_file, predicate, out_file, kind="fasta")
    return out_file


def filter_seqio(in_file, predicate, out_file, kind="fasta"):
    # skip if the output file already exists
    if os.path.exists(out_file):
        return out_file

    with nested(open(in_file, "rU"),
                open(out_file, "w")) as (in_handle, out_handle):
        def output_writer(x):
            return(SeqIO.write(x, out_handle, kind))
            #map(SeqIO.write(x, out_handle, kind),
        map(output_writer, filter(predicate, SeqIO.parse(in_handle, kind)))
    return out_file

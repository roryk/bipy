=====
bipy
=====
This is an ultra alpha set of wrappers around commonly-used bioinformatics tools to
run bioinformatics pipelines using ipython to handle parallelization across
cores and nodes. Borrows heavily from https://github.com/chapmanb/bcbb.


===========
input files
===========
input files should have a consistent naming scheme consisting of a unique
identifier followed by conditions of the sample and ending in _1 (or _2).
for example:

id1_control_batch1_rep1_1.fq

This specifies the id of the file is id1, it is the control condition,
in batch 1 and it is the first replicate of the condition and is the
forward read.

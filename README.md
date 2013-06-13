## bipy
This is an alpha set of wrappers around commonly-used bioinformatics tools to
run bioinformatics pipelines using IPython to handle parallelization across
cores and nodes. Uses several wrappers and ideas from https://github.com/chapmanb/bcbb and will eventually be merged into it.

This documentation is not complete yet.

### installation
    mkdir ~/src; cd ~/src
    pip install numpy
    git clone git@github.com:roryk/bipy.git
    cd bipy
    python setup.py install

### cluster configuration
If you have a standard scheduler/queue setup, setting up bipy to run on your cluster is
straightforward:

```
cluster:
  cores: 4 # number of ipython engines to spin up
  scheduler: lsf
  queue: hsph
```

This will run a total of 4 IPython engines on the Platform LSF scheduler on the queue 'hsph'.
Other values for scheduler which will work are: torque and sge.

If you have a more complicated setup you can still use bipy. You will need to set up
an IPython parallel profile that describes your cluster setup and then use that
profile like this:

```
cluster:
  cores: 4 # number of ipython engines to spin up
  profile: your_profile
  scheduler: lsf
  queue: hsph
```

### quickstart configuration
It will take a small bit of fiddling around to get everything working as there are a couple of pieces that need
to interact together properly to work. You can get going testing everything locally by setting up a default IPython
profile to run on your local machine.

    ipython profile create --parallel --profile=bipy_test

### troubleshooting
If the cluster unit test is failing, something might not be set up properly with IPython.

First, test to see that IPython is functioning correctly:
    ipcluster start --n=4 --profile=bipy_test

You should see something that looks like this:
```
2013-01-23 17:40:56,874 [IPClusterStart] Using existing profile dir: u'/Users/rory/.ipython/profile_bipy_test'
2013-01-23 17:40:56.877 [IPClusterStart] Starting ipcluster with [daemon=False]
2013-01-23 17:40:56.877 [IPClusterStart] Creating pid file: /Users/rory/.ipython/profile_bipy_test/pid/ipcluster.pid
2013-01-23 17:40:56.878 [IPClusterStart] Starting Controller with LocalControllerLauncher
2013-01-23 17:40:57.878 [IPClusterStart] Starting 4 Engines with LocalEngineSetLauncher
2013-01-23 17:41:28.202 [IPClusterStart] Engines appear to have started successfully
```

If you open another terminal on your machine you can manually test to see if you can talk to the engines:
```
rory@clotho:~/cache/bipy (master)$ python
Python 2.7.3 (default, Apr 17 2012, 10:42:40)
[GCC 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2335.15.00)] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> from IPython.parallel import Client
>>> c = Client(profile="bipy_test")
>>> view = c[:]
>>> c.ids
[0, 1, 2, 3]
>>> view.map_sync(lambda x: x**10, range(1, 10))
[1, 1024, 59049, 1048576, 9765625, 60466176, 282475249, 1073741824, 3486784401]
```


### requirements
You should have Bowtie, Tophat, fastqc and samtools installed along with annotation genomes for your organism of interest.
To run the rseqc tools you should also have fetchChromSizes and wigToBigWig
somewhere in your path. You can get those two utilities here in either
the Linux or OSX directory:
http://hgdownload.cse.ucsc.edu/admin/exe/

### input files
FASTQ input files should have a naming scheme where mate pairs have the same
name except for one containing _1 and the other containing _2 and ending
in either .fastq or .fq. For example:
    sample1_rep1_1.fastq
    sample1_rep1_2.fastq
    sample1_rep2_1.fastq
    sample1_rep2_2.fastq

These files should be in the directory specified in the data field of
the YAML file. bipy will find all of the fastq files in all subdirectories
of input_dir as well.

### tool configuration

#### fastqc
FastQC runs some standard quality metrics on your FASTQ files. This stage should be run
before and after adapter trimming so you can see the effect of trimming the adapters on your
sequence composition.

#### cutadapt
The cutadapt tool is what trims adapters and other contaminating sequences from the
ends of reads. For most libraries, this configuration should work fine:
```
  cutadapt:
    program: cutadapt
    chemistry: [truseq]
    trim_polya: True
    options:
      error-rate: 0.1
      quality-cutoff: 20
```
That will search for the first 13 base pairs of the Truseq adapters and trim it
off right right end of the read. It will also search for the reverse complement
of the Truseq adapters. Finally after trimming those off, it will trim off the
polyA tail.

For older libraries, you may need to use "illumina" instead of "truseq". This will
trim off the pre Truseq Illumina adapters. For libraries made with the NextEra chemistry,
use "nextera" to trim the adapters off.

If you need to supply your own adapters you can use a configuration that looks like this:
```
  cutadapt:
    program: cutadapt
    chemistry: [truseq]
    trim_polya: True
    adapters: ["your_adapter_1", "your_adapter_2"]
    options:
      error-rate: 0.1
      quality-cutoff: 20
```
For those adapters if you want to trim the reverse complement, you need to supply that as well.
The above example will trim off both the Truseq adapters as well as your supplied sequence from
the right side of the read.

The *._trimmed.fixed.fastq files are the results of running cutadapt. These have been
trimmed with cutadapt and reads < 20 bases are removed. If the input lanes were paired,
and one of the pair is removed, the other one is as well, and placed in the .singles.fastq
file.

#### tophat
You can pass arbitraty options to Tophat by editing the Tophat portion of the YAML file.
For example if you wanted to run Tophat, using a custom transcriptome located in
/my/custom/transcriptome and run a quick Bowtie alignment you could pass
--b2-very-fast option and the custom transcriptome option to Tophat like this:

```
  tophat:
    name: tophat
    program: tophat
    options:
      b2-very-fast: True
      transcriptome-index: /my/custom/transcriptome
    quality_format: sanger
```
If your reads are of one of the old non-standard Illumina format you can set the
quality_format to "illumina" to handle older-style illumina reads. If you aren't sure,
stick to the sanger format as it is probably the correct format.

### output files
#### quality control
There are several places to look for quality control. In the results/fastqc directory
there are two sets of results, one for the untrimmed lanes and one for the trimmed lanes.
You should compare the results for each lane and see if they make sense; it is possible
there is a contaminating adapter that was missed; if that is true you can add it to cutadapt
and rerun to trim it off.

The second place to look for quality control is to look at the alignments themselves. You can
load the Tophat mapped, sorted and indexed BAM file in
results/tophat/your_file_name_tophat.sorted.sorted.bam into IGV and visually inspect the
reads. You should see relatively even coverage across your genes, if that is not the
case then there may be a problem with your library.

There are two places to look for alignment metrics. The first is under results/rnaseq_metrics;
this is the output from Picard's CollectRnaSeqMetrics tool. In that file are two tab
delimited rows which are important, under ## METRICS. In there you can see total bases that
map to various features of the genome, including coding regions, UTRs and rRNA. It is also
import to look at the 5'-3' bias, you should not see much of a bias in your sample.

The second place to look for mapping metrics is in the rseqc directory. In the
subdirectories under rseqc there are the following quality information:

bam_stat: this is a high level summary of the overall number of records mapped.
junction: these are two pie charts showing the overall class of junctions covered

RPKM_count: this is a rough RPKM calculation. RPKM is known to introduce biases so use this
with caution; this is mostly to get an overall look at your data. If you want to have
normalized data, I suggest loading the count files from htseq-count (see below) into
edgeR and extracting the normalized counts there. The RPKM file you want to look at is
the fixed file, this aggregates the counts by gene_id instead of at the exon level.

RPKM_saturation: the PDF file in this directory is a plot showing what would happen to your
RPKM calculation if you downsampled your data. This is not all that useful, it has been
shown in several places that if you are doing a standard gene-level differential expression
experiment, you can do well with 10-15 million usable reads. If you have many more reads than
that, it is a good idea to run more replicates rather than sequence the same sample more.

saturation: the PDF file in this directory shows what happens to your ability to detect
known and novel splice junctions as you drop the read depth, another way of determining
what you would lose if you sequenced less.

#### downstream analysis
The results/tophat/your_file/your_file.sorted.sorted.bam file can be used to
visualize the reads. For performing differential expression, you can find a list
of gene-level counts in results/htseq-count. You can use these gene-level counts
to perform a differential expression analysis using edgeR or DESeq.  Example R
markdown scripts to do this analysis are in scripts/deseq.Rmd and
scripts/edgeR_paired.Rmd.

### recommendations
Set test_pipeline: True in the YAML configuration file to run the whole pipeline
on a small subset of your data first, to make sure it works correctly.

### contributors
Thank you to Brad Chapman, Oliver Hoffman, John Hutchinson, Sara Dempster, Giles
Hall and Georgios Marnellos for their contributions.

### bipy
This is an alpha set of wrappers around commonly-used bioinformatics tools to
run bioinformatics pipelines using IPython to handle parallelization across
cores and nodes. Uses several wrappers and ideas from https://github.com/chapmanb/bcbb and will eventually be merged into it.

This documentation is not complete yet.

#### installation
    mkdir ~/src; cd ~/src
    pip install numpy
    git clone git@github.com:roryk/bipy.git
    cd bipy
    python setup.py install

#### cluster configuration
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

#### quickstart configuration
It will take a small bit of fiddling around to get everything working as there are a couple of pieces that need
to interact together properly to work. You can get going testing everything locally by setting up a default IPython
profile to run on your local machine.

    ipython profile create --parallel --profile=bipy_test

#### troubleshooting
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


#### requirements
You should have Bowtie, Tophat, fastqc and samtools installed along with annotation genomes for your organism of interest.
To run the rseqc tools you should also have fetchChromSizes and wigToBigWig
somewhere in your path. You can get those two utilities here in either
the Linux or OSX directory:
http://hgdownload.cse.ucsc.edu/admin/exe/

#### input files
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

#### tool information
##### cutadapt
For most libraries, this configuration should work fine:
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


#### recommendations
Set test_pipeline: True in the YAML configuration file to run the whole pipeline
on a small subset of your data first, to make sure it works correctly.

#### Contributors
Thank you to Brad Chapman, Oliver Hoffman, John Hutchinson, Sara Dempster, Giles
Hall and Georgios Marnellos for their contributions.

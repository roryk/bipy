# bipy
This is an alpha set of wrappers around commonly-used bioinformatics tools to
run bioinformatics pipelines using IPython to handle parallelization across
cores and nodes. Uses several wrappers and ideas from https://github.com/chapmanb/bcbb and
will eventually be merged into it.

This documentation is not complete yet.

## installation
    mkdir ~/src; cd ~/src
    pip install numpy
    git clone git@github.com:roryk/bipy.git
    cd bipy
    python setup.py install
    
## quickstart configuration
It will take a small bit of fiddling around to get everything working as there are a couple of pieces that need
to interact together properly to work. You can get going testing everything locally by setting up a default IPython
profile to run on your local machine.

    ipython profile create --parallel --profile=bipy_test

## troubleshooting
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


## requirements
You should have Bowtie, Tophat, fastqc and samtools installed along with annotation genomes for your
organism of interest.

## input files
FASTQ input files should have a naming scheme where mate pairs have the same
name except for one containing _1 and the other containing _2 and ending
in either .fastq or .fq. For example:
    sample1_rep1_1.fastq
    sample1_rep1_2.fastq
    sample1_rep2_1.fastq
    sample1_rep2_2.fastq

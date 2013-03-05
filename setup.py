#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name="bipy",
      version="0.1.0",
      author="Rory Kirchner",
      author_email="rory.kirchner@gmail.com",
      description="Simple analysis pipelines",
      url="https://github.com/roryk/bipy",
      license="MIT",
      namespace_packages=["bipy"],
      packages=find_packages(),
      package_data={'bipy': ['toolbox/data/*']},
      install_requires=[
          "numpy >= 1.6.2",
          "ipython >= 0.13.1",
          "bcbio-nextgen",
          "biopython >= 1.60",
          "cutadapt >= 1.2.1",
          "pandas >= 0.1.0",
          "HTSeq >= 0.5.3p9",
          "pysam == 0.7",
          "sh >= 1.0.8",
          "rseqc == 2.3.5"])

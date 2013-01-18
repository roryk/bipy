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
      dependency_links = ['https://github.com/roryk/sh/tarball/master#egg=sh-1.07'],
      package_data={'bipy': ['toolbox/data/*']},
      install_requires=[
          "cython",
          "ipython >= 0.11",
          "bcbio-nextgen",
          "biopython == 1.60",
          "cutadapt == 1.2.1",
          "pandas == 0.1.0",
          "HTSeq == 0.5.3p9",
          "sh"])

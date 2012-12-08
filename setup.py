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
          "ipython >= 0.11",
          "bcbio",
          "sh == 1.07",
          "biopython == 1.60",
          "cutadapt == 1.2.1"])

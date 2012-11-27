#!/usr/bin/env python
from setuptools import setup, find_packages

setup(name="bipy",
      version="0.0.1",
      author="Rory Kirchner",
      author_email="rory.kirchner@gmail.com",
      description="Simple analysis pipelines",
      url="https://github.com/roryk/bipy",
      license="MIT",
      namespace_packages=["bipy"],
      packages=find_packages(),
      install_requires=[
          "ipython >= 0.11",
          "bcbio",
          "sh"])

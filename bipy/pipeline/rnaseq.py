"""
example pipelines for running RNA-seq experiments

"""
import yaml
from bipy.utils import nested_lookup
from functools import partial
from bipy.pipeline import AbstractStage
from bipy.toolbox.fastqc import FastQCStage



class AbstractPipeline(object):


    def __init__(self, config_file):
        with open(config_file) as in_handle:
            self.config = yaml.load(in_handle)

    def validate_config(self):
        """
        make sure the necessary keys are available in the config file
        """
        WANT = [("dir")]
        have = map(partial(nested_lookup, self.config), WANT)
        return all(have)


class PipelineBuilder(AbstractPipeline):

    def __init__(self, config):
        self.config = config
        self.stages = self.config["stage"]
        self.to_run = self.config["run"]



class RNASeqDEPipeline(AbstractPipeline):
    STAGES = ["fastqc", "cutadapt", "fastqc", "tophat", "htseq-count",
             "deseq"]

    def __init__(self, config_file, in_files=None):
        super(RNASeqDEPipeline, self).__init__(config_file)

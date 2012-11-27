"""
general functions for handling common pipeline steps
"""
import networkx as nx
import pylab
from bipy.pipeline.stages import STAGE_LOOKUP
from bipy.cluster import start_cluster, stop_cluster
from bcbio.utils import safe_makedir, file_exists
from bipy.log import logger, setup_logging
import datetime
import os
import yaml
import abc
from Queue import Queue


def setup_pipeline(config):
    """
    creates output directories
    and performs some minor validation on the configuration file
    """
    # make initial directories
    if "dir" not in config:
        logger.error("'dir' must be in config file, see example "
                     " configurations.")
        exit(-1)
    config = _setup_config(config)
    map(safe_makedir, config["dir"].values())

    _write_config(config)

    return config


def _write_config(config):
    """
    output a yaml file of the configuration for this run in the results
    directory
    """
    yaml_file = os.path.join(config["dir"].get("results", "results"),
                             "config.yaml")
    with open(yaml_file, "w") as out_handle:
        out_handle.write(yaml.dump(config))


def _setup_config(config):
    """
    configuration file setup
    """
    if "pipeline" in config:
        # add timestamp for the run
        if config["pipeline"].get("timestamp", None):
            time_fmt = '%Y-%m-%d-%H-%M-%S'
            now = datetime.datetime.now().strftime(time_fmt)
            stamped_dir = os.path.join(config["dir"]["results"], now)
            config["dir"]["results"] = stamped_dir

    return config


class AbstractRunner(object):

    @abc.abstractmethod
    def run(self):
        """Runs the command"""
        return

    @abc.abstractmethod
    def output(self):
        """Returns summary output of running the command in a dictionary"""
        return

    @abc.abstractmethod
    def command(self):
        """Returns the actual command line for running the command as a string"""
        return

    @abc.abstractmethod
    def accepts(self):
        """Returns a list of extensions the command accepts for input files"""
        return

    @abc.abstractmethod
    def next(self):
        """Returns the files that would be used post this stage of processing"""
        return


class GenericPipeline(object):

    def __init__(self, in_files, config):
        self.in_files = in_files
        self.config = config
        setup_logging(self.config)
        start_cluster(self.config)
        from bipy.cluster import view
        self.view = view
        self.curr_files = [in_files]
        self.to_run = Queue()
        map(self.to_run.put, self.setup_stages(config["run"]))

    def setup_stages(self, stages):
        stage_classes = map(STAGE_LOOKUP.get, stages)
        stage_objects = [stage(self.config) for stage in stage_classes]
        return stage_objects

    def validate_config(self):
        """
        generic config file validator for things that should be in
        all pipeline configuration files
        """
        WANT = ["stage", "run", "dir", "cluster"]
        for w in WANT:
            if w not in self.config:
                raise ValueError('Could not find %s in the config file')

    def get_next_files(self):
        return self.curr_files[-1]

    def set_next_files(self, files):
        self.curr_files.append(files)

    def process(self):
        """
        simple pipeline processor. can't handle tasks where stages
        are forked off but can ignore the output for any stage.

        """
        # keep pulling items off of the queue until it is empty
        while not self.to_run.empty():
            stage = self.to_run.get()
            next_files = self.get_next_files()

        # run in a non-blocked manner since we are not waiting on them
        if not stage.keep:
            self.view.map(stage, next_files, block=False)
            yield None
        else:
            processed = self.view.map(stage, self.get_next_files())
            yield processed
            self.set_next_files(processed)

    def process_next_stage(self):
        stage = self.to_run.get()
        if not stage:
            return None
        processed = self.view.map(stage, self.get_next_files)
        self.set_next_files(processed)
        return processed


"""
parse a configuration dict to set up the pipeline
"""
import networkx as nx
import pylab
import yaml
import os

class Node(object):

    def __init__(self, name, block):
        self.name = name
        if block:
            self.block = "block"
        else:
            self.block = "unblock"

    def __str__(self):
        return "%s - %s" % (self.name, self.block)


class PipelineGraph(object):

    def __init__(self, config):
        self.config = config
        #self.dag = pgv.AGraph(strict=True, directed=True)
        self.dag = nx.DiGraph()
        self.root = Node("start", True)
        self.dag.add_node(self.root)
        self._parse(config.get("run", {}), self.root)

    def _parse(self, run, parent):
        if parent not in self.dag:
            self.dag.add_node(parent)
        for item in run:
            if isinstance(item, str):
                child = Node(item, False)
                self.dag.add_node(child)
                self.dag.add_edge(parent, child)
            elif isinstance(item, dict):
                child = Node(item.keys()[0], True)
                self.dag.add_node(child)
                self.dag.add_edge(parent, child)
                self._parse(item.values()[0], child)

    def graph(self, out_file=None):
        if out_file is None:
            doc_dir = self.config["dir"].get("doc")
            out_file = os.path.join(doc_dir, "pipeline_viz.png")
        if file_exists(out_file):
            return out_file
        pd = nx.to_pydot(self.dag)
        pd.write_png(out_file, prog="dot")

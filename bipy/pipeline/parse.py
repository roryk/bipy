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







config_file = "/Users/rory/cache/bipy/test/pipeviz/test_pipeviz.yaml"
with open(config_file) as in_handle:
    config = yaml.load(in_handle)

p = Parser(config)
p.parse(p.config["run"], p.root)
pd = nx.to_pydot(p.dag)
out_file = "/Users/rory/tmp/lol.png"
pd.write_png(out_file, prog="dot")
im = pylab.imread(out_file)
pylab.imshow(im)
pylab.axis('off')
#pylab.show()

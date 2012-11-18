"""
parse a configuration dict to set up the pipeline
"""
import networkx as nx


class Parser(object):

    def __init__(self, config):
        self.config = config
        self.dag = nx.Graph()

        # XXX NOT COMPLETE
    def create_layer(self, parent, nodes):
        self.dag.add_nodes_from(nodes)

    def parse(self):
        pass

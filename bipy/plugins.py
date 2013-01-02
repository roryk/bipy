"""
plugin handler

automatically finds all classes inheriting from AbstractStage.

inherit from bipy.pipeline.stage.AbstractStage and implement
__init__ and __call__ to create a Stage which can participate
in the pipeline process

place the .py file into the plugin directory specified in your
pipelines config.yaml file

mostly lifted from Giles Hall https://github.com/vishnubob
"""

import sys
import os
import traceback
import imp
from bipy.pipeline.stages import AbstractStage
from bipy.utils import get_in
from bipy.log import logger

PluginDirectory = os.path.join(os.path.split(__file__)[0], "toolbox")
#PluginDirectory = (os.path.split(__file__)[0] or '.')
PackagePrefix = str.join('.', __name__.split('.')[:1])
sys.path.append(PluginDirectory)


def import_(modname, prefix='', globals=None, loals=None, fromlist=None):
    # Fast path: see if the module has already been imported.
    if prefix:
        name = prefix + "." + modname
    else:
        name = modname
    try:
        return sys.modules[name]
    except KeyError:
        pass

    # If any of the following calls raises an exception,
    # there's a problem we can't handle -- let the caller handle it.
    fp, pathname, description = imp.find_module(modname)

    try:
        return imp.load_module(name, fp, pathname, description)
    finally:
        # Since we may exit via an exception, close fp explicitly.
        if fp:
            fp.close()


class StageRepository(object):
    """
    builds a repository of all of the stage objects. loads all of the
    AbstractStage objects in the toolbox directory and adds to it anything
    in the plugins directory.

    will overwrite the main toolbox stages with the ones in the plugins
    directory

    """

    def __init__(self, config):
        self.config = config
        self.plugins = {}
        self.scan(get_in(config, "dir", "plugins"))
        #        if get_in(config, ("dir", "plugins")):
        #    self.scan(get_in(config, ("dir", "plugins")))

    def scan(self, plugin_dir=None):

        files = os.listdir(PluginDirectory)
        if plugin_dir:
            sys.path.append(plugin_dir)
            files += os.listdir(plugin_dir)

        plugins = []
        for fn in files:
            if fn.endswith('.py'):
                plugins.append(fn)
        mods = [fn.split('.')[0] for fn in plugins]
        # build the map of plugins
        for modname in mods:
            try:
                mod = import_(modname, PackagePrefix)
            except:
                logger.error("Error loading plugin: %s" % modname)
                traceback.print_exc()
                continue
            self.scan_module(mod)

    def scan_module(self, mod):
        ns = dir(mod)
        for name in ns:
            # skip over the base class
            thing = getattr(mod, name)
            if thing == AbstractStage:
                continue
            if type(thing) is not type:
                continue
            if issubclass(thing, AbstractStage):
                # found a plugin!
                self.plugins[thing.stage] = thing

    def __getitem__(self, key):
        return self.plugins.get(key, None)

"""
Repository = None


def get_plugins():
    global Repository
    if Repository is None:
        Repository = StageRepository()
    return Repository
    """

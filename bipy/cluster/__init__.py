from IPython.parallel import Client
import subprocess
from time import sleep
from bipy.log import logger
import atexit
import yaml
import uuid

DEFAULT_CLUSTER_TIMEOUT = 300
DEFAULT_DELAY = 10

# XXX come up with a saner set of defaults, especially for the cluster.
DEFAULT_CONFIG = {"dir": {"results": "results",
                           "tmp": "tmp",
                           "log": "log"},
                  "cluster": {"profile": "default",
                              "cores": 4}}

# these should be modified once, at the beginning of the pipeline and
# be read only after that
config = DEFAULT_CONFIG
cluster = None
client = None
view = None
direct_view = None


def load_config(config_file):
    """ load the config file, overwriting the defaults if it has them """
    global config
    with open(config_file) as in_handle:
        cur_config = yaml.load(in_handle)
        # this code lifted from bcbio.utils
        # https://github.com/chapmanb/bcbb/blob/master/nextgen/bcbio/utils.py
        for k, v in cur_config.iteriterms():
            if k in config and isinstance(config[k], dict):
                config[k].update(v)
            else:
                config[k] = v


@atexit.register
def stop_cluster():
    if cluster:
        cluster.stop()


def start_cluster(cluster_config):
    global cluster, view, client, direct_view
    cluster = Cluster(**cluster_config["cluster"])
    logger.info("Starting the cluster with %d nodes." % (cluster.n))
    cluster.start()
    sleep(cluster.delay)

    # only continue when the cluster is completely up
    slept = 0
    while(not cluster.is_up()):
        sleep(cluster.delay)
        slept = slept + cluster.delay
        if(slept > cluster_config["cluster"].get("timeout",
                                                 DEFAULT_CLUSTER_TIMEOUT)):
            logger.error("Cluster startup timed out.")
            cluster.stop()
            exit(-1)
    # only continue if at least one engine is up

    logger.info("Cluster up.")
    client = cluster.client()
    view = cluster.view()
    direct_view = cluster.direct_view()
    engine_config = cluster_config.copy()
    engine_config["engine_log"] = True
    direct_view['config'] = engine_config
    direct_view.execute('from bipy.log import setup_logging')
    direct_view.execute('setup_logging(config)')


class Cluster(object):
    def __init__(self, **kwargs):
        self.profile = kwargs.get("profile", "default")
        self.n = kwargs.get("cores", 1)
        self.delay = kwargs.get("delay", DEFAULT_DELAY)
        self._client = None
        self._view = None
        self._direct_view = None
        self._work = kwargs.get("work", ".")
        self._log_level = kwargs.get("log_level", 30)
        self._cluster_id = str(uuid.uuid1())

    def start(self):
        """starts the cluster and connects the client to the controller"""
        narg = "--n=%d" % (self.n)
        parg = "--profile=%s" % (self.profile)
        """
         XXX: in the future, add "--cluster-id=" + self._cluster_id to
         this, to run each new cluster with a different ID, so we
         can reuse the same profile. right now there is a bug in
         ipython that doesn't support this"""
        return_code = subprocess.call(["ipcluster", "start",
                                       "--daemonize=True",
                                       "--delay=" + str(self.delay),
                                       "--log-level=" + str(self._log_level),
                                       "--debug",
                                       narg, parg])

    def client(self):
        """ returns a handle to the client """

        # add cluster_id=self._cluster_id to this call when the bug
        # is fixed in iPython
        if not self._client:
            self._client = Client(profile=self.profile)
            return self._client
        return self._client

    def new_client(self):
        if self._client:
            self._client.close()
        self._client = Client(profile=self.profile)

    def view(self):
        """ returns a blocking, load balanced view to the cluster engines """
        if self._view:
            return self._view

        if not self._client:
            self._client = Client(profile = self.profile)

        self._view = self._client.load_balanced_view()
        self._view.block = True
        return self._view

    def direct_view(self):
        if self._direct_view:
            return self._direct_view
        if not self.client:
            self._client = Client(profile = self.profile)
        self._direct_view = self._client[:]
        return self._direct_view

    def stop(self):
        parg = "--profile=%s" % (self.profile)
        # add carg = "--cluster-id=%s" % (self._cluster_id) when
        # this gets fixed in iPython
        return_code = subprocess.call(["ipcluster", "stop", parg])

    def is_up(self):
        """ returns True if the cluster is completely up and false otherwise """
        try:
            up = len(self.client().ids)
        except IOError:
            logger.info("Waiting for the controller to come up.")
            return False
        else:
            not_up = self.n - up
            if not_up > 0:
                logger.info("Waiting for %d engines to come up." %(not_up))
                return False
            else:
                return True


def mappable_function(x):
    logger.error("This is an error.")
    logger.info("This is info.")
    return x ** 10

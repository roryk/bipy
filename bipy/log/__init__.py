"""Utility functionality for logging.
"""
import os
import sys
import logging
from bcbio import utils
from bipy.utils import get_in
import zmq
from zmq.log.handlers import PUBHandler
import socket
from IPython.utils.path import get_security_file
import zmq
from IPython.zmq.session import Session
from IPython.parallel.util import disambiguate_url
from IPython.utils.py3compat import str_to_bytes
from IPython.utils.path import get_security_file
import sh
import time

LOG_NAME = "bipy"

logger = logging.getLogger(LOG_NAME)

def setup_logging(config):
    logger.setLevel(logging.INFO)
    if not logger.handlers:
        formatter = logging.Formatter('[%(asctime)s] %(message)s')
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        log_dir = config.get("log_dir", None)
        if log_dir:
            logfile = os.path.join(utils.safe_makedir(log_dir),
                                   "{0}.log".format(LOG_NAME))
            handler = logging.FileHandler(logfile)
            handler.setFormatter(formatter)
            logger.addHandler(handler)

import logbook
logger2 = logbook.Logger(LOG_NAME)

def create_log_handler(config):
    log_dir = config.get("log_dir", None)
    email = config.get("email", None)

    if log_dir:
        utils.safe_makedir(log_dir)
        handler = logbook.FileHandler(os.path.join(log_dir, "%s.log" % LOG_NAME))
    else:
        handler = logbook.StreamHandler(sys.stdout)

    if email:
        handler = logbook.MailHandler(email, [email],
                                      format_string=u'''Subject: [BCBB pipeline] {record.extra[run]} \n\n {record.message}''',
                                      level='INFO', bubble = True)
    return handler


def setup_ipython_logging(connect_to=None):
    """
    set up a zeromq publisher bound to a port.
    returns the address + port to talk to the publisher

    """
    # find the ip address of this machine
    if not connect_to:
        try:
            address = socket.gethostbyname(socket.gethostname())
        except:
            address = "127.0.0.1"

            connect_to = "tcp://" + address
            context = zmq.Context()
            a = context.socket(zmq.PUB)
            port = a.bind_to_random_port(connect_to)
            handler = BipyHandler(a)
            logger.addHandler(handler)
            a.send(__file__)
            return connect_to + ":" + str(port)

    context = zmq.Context()
    a = context.socket(zmq.PUB)
    a.connect(connect_to)
    a.send(__file__)
    handler = BipyHandler(a)
    logger.addHandler(handler)
    return connect_to


class BipyHandler(PUBHandler):
    """
    simple handler to format log messages in the same manner as bcbio
    and publish them to a socket

    """

    def __init__(self, sock, context=None):
        PUBHandler.__init__(self, sock)
        self.socket = sock
        self.ctx = context
        self.formatter = logging.Formatter('[%(asctime)s] %(message)s')

    def format(self, record):
        return self.formatter.format(record)

    def emit(self, record):
        self.socket.send(self.format(record))


def launch_iowatcher(config):
    """
    launch a script to watch traffic on the IOPub channel and
    log to a file

    """
    log_dir = config.get("log_dir", None)
    cluster_id = get_in(config, ("cluster", "cluster_id"), None)
    if cluster_id:
        log_name = LOG_NAME + "-" + cluster_id + ".log"
    else:
        log_name = LOG_NAME + ".log"

    log_file = os.path.join(log_dir, log_name)

    sec_file = _security_file_name(config)
    logger.info("Starting ZMQ IOPub watcher from %s and writing to %s"
                % (sec_file, log_file))
    sh.python("/Users/rory/cache/bipy/scripts/iopubwatcher.py",
              sec_file, log_file, _bg=True)
    time.sleep(5)

def _security_file_name(config):
    """
    Returns the absolute path of a security file, including the cluster-id
    in the name if it has been specified

    """
    profile = get_in(config, ("cluster", "profile"), "default")
    cluster_id = get_in(config, ("cluster", "cluster_id"), None)
    sec_file = "ipcontroller-client.json"
    if cluster_id:
        sec_file = "ipcontroller-" + "cluster_id" + "-client.json"
    else:
        sec_file = "ipcontroller-client.json"

    return get_security_file(sec_file, profile)

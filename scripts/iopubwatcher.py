"""A script for watching all traffic on the IOPub channel (stdout/stderr/pyerr)
of engines.

Cribbed mostly from here, from MinRK:
https://raw.github.com/ipython/ipython/master/docs/examples/parallel/iopubwatcher.py

Authors
-------
* MinRK
"""

import sys
import json
import zmq
from IPython.zmq.session import Session
from IPython.parallel.util import disambiguate_url
from IPython.utils.py3compat import str_to_bytes


def main(connection_file, log_file):
    """watch iopub channel, and print messages"""

    ctx = zmq.Context.instance()

    with open(connection_file) as f:
        cfg = json.loads(f.read())

    log = open(log_file, "w")
    print log

    location = cfg['location']
    reg_url = cfg['interface']
    session = Session(key=str_to_bytes(cfg['exec_key']))

    #query = ctx.socket(zmq.DEALER)
    #query.connect(disambiguate_url(reg_url, location))
    #session.send(query, "connection_request")
    #idents, msg = session.recv(query, mode=0)
    #c = msg['content']
    #iopub_url = disambiguate_url(c['iopub'], location)
    sub = ctx.socket(zmq.SUB)
    # This will subscribe to all messages:
    sub.setsockopt(zmq.SUBSCRIBE, b'')
    iopub_url = reg_url + ":" + str(cfg['iopub'])
    # replace with b'' with b'engine.1.stdout' to subscribe only to
    # engine 1's stdout
    # 0MQ subscriptions are simple 'foo*' matches, so 'engine.1.' subscribes
    # to everything from engine 1, but there is no way to subscribe to
    # just stdout from everyone.
    # multiple calls to subscribe will add subscriptions, e.g. to subscribe to
    # engine 1's stderr and engine 2's stdout:
    # sub.setsockopt(zmq.SUBSCRIBE, b'engine.1.stderr')
    # sub.setsockopt(zmq.SUBSCRIBE, b'engine.2.stdout')
    sub.connect(iopub_url)
    while True:
        try:
            idents, msg = session.recv(sub, mode=0)
        except KeyboardInterrupt:
            return
        # ident always length 1 here
        topic = idents[0]
        if msg['msg_type'] == 'stream':
            # stdout/stderr
            # stream names are in msg['content']['name'], if you want to handle
            # them differently
            log.write("%s: %s" % (topic, msg['content']['data']))
            log.flush()
        elif msg['msg_type'] == 'pyerr':
            # Python traceback
            c = msg['content']
            log.write(topic + ':')
            log.flush()
            for line in c['traceback']:
                # indent lines
                log.write('    ' + line)
                log.flush()

if __name__ == '__main__':
    cf = sys.argv[1]
    if len(sys.argv) == 2:
        lf = sys.argv[2]
    else:
        lf = "bipy.log"

    main(cf, lf)

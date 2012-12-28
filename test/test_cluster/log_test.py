from IPython.config.application import Application
import logging

logger = logging.getLogger('bipy')

def logging_test(x):
    print("From inside engine: %s" % (str(x**10)))
    return x**10

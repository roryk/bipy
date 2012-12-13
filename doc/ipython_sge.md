============================================================
Guide for getting ipython going with platform lsf on odyssey
============================================================

Create a parallel profile named sge:
ipython profile create --parallel --profile=sge

Edit ~/.config/.ipython/profile_lsf/ipcluster_config.py
 c.IPClusterStart.controller_launcher_class = 'SGE'
 c.IPClusterStart.engine_launcher_class = 'SGE'
 c.SGELauncher.queue = u'your_default_queue_name'

Edit ~/.ipython/profile_odyssey_lsf/ipcontroller_config.py
 c.HubFactory.ip = '*'

Finally launch a mini cluster from the command line.
 ipcluster start -n 2 --profile=sge

- Now you should be able to connect to it and go from inside python:

 from IPython.parallel import Client

 c = Client(profile="sge")
 c.ids
 c[:].apply_sync(lambda: "Hello, World")

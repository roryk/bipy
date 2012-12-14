#Guide for getting ipython talking to SGE 

## create a parallel profile named sge
ipython profile create --parallel --profile=sge

## edit ~/.config/ipython/profile_sge/ipcluster_config.py
```python
    c.IPClusterStart.controller_launcher_class = 'SGE'
    c.IPClusterStart.engine_launcher_class = 'SGE'
    c.SGELauncher.queue = u'your_default_queue_name'
```

## edit ~/.config/ipython/profile_sge/ipcontroller_config.py
```python
    c.HubFactory.ip = '*'
```

## finally launch a mini cluster from the command line.
    ipcluster start -n 2 --profile=sge

## now you should be able to connect to it and do this from inside python:
```python
    from IPython.parallel import Client

    c = Client(profile="sge")
    c.ids
    c[:].apply_sync(lambda: "Hello, World")
```

You should be able to look at the jobs and see three running, one ipcontroller
and two engines. When you run c[:].apply_sync(lambda: "Hello, World") that executes
on the two engines.

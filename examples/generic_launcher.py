import argparse
from bipy.cluster import start_cluster, stop_cluster
from bipy.log import setup_logging
import sys

def main():
    view.map(lambda x: "Hello World", xrange(10))

def _args_valid_for_scheduler(args):
    return args.scheduler and args.queue

def _args_valid_for_local(args):
    return not args.scheduler and not args.queue

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='generic launcher')
    parser.add_argument('--profile', required=True,
                        help="IPython profile name to use")
    parser.add_argument('--cores', required=True,
                        help="Number of IPython engines to start.")
    parser.add_argument('--queue',
                        help="Name of queue to use.")
    parser.add_argument('--scheduler', default="",
                        help="Name of scheduler to use (LSF or SGE)")
    args = parser.parse_args()

    cluster_config = {"cluster":
                      {"profile": args.profile,
                       "cores": int(args.cores),
                       "queue": args.queue,
                       "scheduler": args.scheduler}}
    setup_logging(cluster_config)
    if _args_valid_for_scheduler(args) or _args_valid_for_local(args):
        start_cluster(cluster_config)
        from bipy.cluster import view
        main()
        stop_cluster()
    else:
        print parser.usage()
        sys.exit(1)

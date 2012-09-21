from bipy import cluster
import yaml
import sys

def main():
    f = open(sys.argv[1])
    config = yaml.load(f)
    f.close()
    cluster.run(config)


if __name__ == "__main__":
    main()

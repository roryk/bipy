import yaml
import os
import collections
from bcbio import utils
import subprocess
import functools
import stat
import difflib
from itertools import repeat
from functools import reduce


def in2out(in_file, word, transform=True, out_dir=None):
    """
    creates an out_file name from an in_file name by adding word to the
    in_file name and writing it to the output_directory if specified.
    if transform is True, the word replaces the extension of in_file.
    if False it adds it to the stem
    transform: in2out("test.bam", "sam", transform=True) -> "test.sam"
    non transform: in2out("test.bam", "sam", transform=False) -> "test_bam.sam"

    """
    (base, ext) = os.path.splitext(in_file)
    if out_dir:
        base = os.path.join(out_dir, os.path.basename(base))

    if transform:
        return "%s.%s" % (base, word)

    else:
        return "%s_%s.%s" % (base, ext, word)



def get_stem(filename):
    return os.path.basename(remove_suffix(filename))


def combine_pairs(input_files):
    """ calls files pairs if they are completely the same except
    for one has _1 and the other has _2 returns a list of tuples
    of pairs or singles """
    PAIR_FILE_IDENTIFIERS = ["1", "2"]

    pairs = []
    used = []
    for in_file in input_files:
        if in_file in used:
            continue
        for comp_file in input_files:
            if comp_file in used:
                continue
            s = difflib.SequenceMatcher(a=in_file, b=comp_file)
            blocks = s.get_matching_blocks()
            # length 3 means on match in the middle of the string
            if len(s.get_matching_blocks()) is not 3:
                continue
            if comp_file[blocks[0][2]] in PAIR_FILE_IDENTIFIERS:
                used.append(in_file)
                used.append(comp_file)
                pairs.append([in_file, comp_file])
                break
        if in_file not in used:
            pairs.append([in_file])
            used.append(in_file)

    return pairs


def freeze_files(files, directory):
    """ freezes a list of files, copying them to directory and
    setting them to be read only """
    new_files = [os.path.join(directory, x) for x in
                 map(os.path.basename, files)]
    [os.chmod(x, stat.S_IREAD | stat.S_IRGRP) for x in new_files]
    return new_files


# XXX I'm not 100% sure if this actually works, I forget use with caution
def memoize_outfile_to_dir(res_dir, ext):
    """Creates outfile from input file and ext, running if outfile not present.

    This requires a standard function usage. The first arg, or kwarg 'in_file', needs
    to be the input file that is being processed. The output name is created with the
    provided ext relative to the input. The function is only run if the created
    out_file is not present.
    Memoize_outfile was lifted from bcbio.utils (Brad Chapman), added the
    writing to a directory portion
    """
    utils.safe_makedir(res_dir)
    def decor(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            if len(args) > 0:
                in_file = args[0]
            else:
                in_file = kwargs['in_file']
            out_fname = "%s%s" % (os.path.splitext(in_file)[0], ext)
            out_file = os.path.join(res_dir, out_fname)
            if not os.path.exists(out_file) or os.path.getsize(out_file) == 0:
                kwargs['out_file'] = out_file
                f(*args, **kwargs)
            return out_file
        return wrapper
    return decor


def prepare_ref_file(ref, config):
    """Get a reference file, either by URL or locally.
    Lifted from Brad Chapman
    """
    url = ref.get("url", None)
    if url:
        ref_file = _download_ref(url, config["dir"]["ref"])
    else:
        ref_file = ref.get("file", None)
    assert ref_file is not None and os.path.exists(ref_file), ref_file
    return ref_file


def _download_ref(url, ref_dir):
    #Lifted from Brad Chapman
    dl_file = os.path.basename(url)
    ref_file = None
    for supported_ext, extract_cmd in [(".gz", "gunzip"),
                                       (".tgz", ("tar", "zxvf"))]:
        if dl_file.endswith(supported_ext):
            ref_file = os.path.join(ref_dir, dl_file[:-len(supported_ext)])
            break
    assert ref_file is not None, url
    if not os.path.exists(ref_file):
        with utils.chdir(ref_dir):
            cl = ["wget", url]
            subprocess.check_call(cl)
            cl = list(flatten([extract_cmd, dl_file]))
            subprocess.check_call(cl)
    return ref_file


def build_results_dir(stage_config, config):
    outdir = os.path.join(config["dir"]["results"],
                          stage_config["name"])
    return outdir


def transform_infile(filename, stage, delim="."):
    return replace_suffix(filename, stage.get("name", "transformed"), delim)


def filter_infile(filename, stage, delim="."):
    return append_stem(filename, stage.get("name", "filtered"), delim)


def remove_suffix(filename):
    filename, extension = os.path.splitext(filename)
    return filename


def append_stem(filename, word, delim="."):
    """ returns a filename with word appended to the stem
    example: /path/to/test.run.sam -> /path/to/test.run.sorted.sam """
    dirname = os.path.dirname(filename)
    filename = os.path.basename(filename)
    fsplit = filename.split(delim)
    fsplit.insert(len(fsplit) - 1, word)
    return os.path.join(dirname, delim.join(fsplit))


def replace_suffix(filename, suffix, delim="."):
    """ returns a filename with the suffix replaced
    with a new suffix
    exampoe: /path/to/test.run.sam -> /path/to/test.run.fa"""
    dirname = os.path.dirname(filename)
    filename = os.path.basename(filename)
    fsplit = filename.split(delim)[:-1]
    fsplit.append(suffix)
    return os.path.join(dirname, delim.join(fsplit))


def flatten_options(config):
    """ returns a list of pairs of options from a config file; returns
    an empty list if no options are found """
    return(flatten(config.get("options", {}.items())))


def flatten(l):
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def rfind_key(d, item):
    for k,v in d.items():
        if k == item:
            return (k, v)
        elif isinstance(v, dict):
            rfind_key(v, item)


def rfind_value(d, item):
    for k, v in d.items():
        if v == item:
            return (k, v)
        elif isinstance(v, dict):
            rfind_value(v, item)


def validate_config(config, data=True):
    # add the full path to any executables to avoid confusion later
    for k, v in config["program"].items():
        config["program"][k] = _add_full_path(v)

    if not all([os.path.exists(x) for x in config["program"].values()]):
        print("Cannot find all of the files or directories listed "
              "in 'program'.")
        exit(-1)

    if not all(x in get_stages(config) for x in config["run"]):
        print("Cannot find all of the stages listed to run.")
        exit(-1)

    if data and not all([os.path.exists(x) for x in config["input"]]):
        print("All of the input files do not exist.")
        exit(-1)

    return(config)


def get_stages(config):
    return config["stage"].keys()


def _add_full_path(program):
    """ add the full path to any executables """
    new_path = which(program)
    if new_path:
        program_config = new_path
    return program_config


def find_one_if(predicate, coll):
    """ returns a single entry in a collection if it contains exactly
    one true instance of a predicate and None otherwise """
    results = filter(predicate, coll)
    if len(results) != 1:
        return None
    else:
        print results
        return results[0]


def load_yaml(filename):
    with open(filename) as f:
        body = yaml.load(f)
        return body


def which(program):
    """ returns the path to an executable"""

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


"""
simple file tracking to keep track of the transformations a file has made
just playing with it for now
"""

class FileWithHistory(object):
    """ this tracks transformations of an input file so we can get the
    metadata associated with it at any point """

    def __init__(self, in_file, **kwargs):
        self._name_list = [in_file]
        [setattr(self, key, value) for key, value in kwargs.items()]

    def curr_file(self):
        return self._name_list[len(self._name_list)]

    def add_file(self, in_file):
        self._name_list.append(in_file)

    def was_file(self, in_file):
        return in_file in self._name_list


def update_tracker(tracker, curr_files):
    new_tracker = [track.add_file(new_file) for (track, new_file)
                   in zip(tracker, curr_files)]
    return new_tracker


def curr_files(tracker):
    return [x.curr_file() for x in tracker]


def dict_to_vectors(d):
    """ convert a dictionary to two vectors. returns a tuple of
    (key, value). For example: {"a": [1, 2, 3], "b": [4, 5]} returns
    (["a" "a" "a" "b" "b"], [1 2 3 4 5])

    """

    ks = []
    vs = []
    for k, v in d.items():
        ks += list(repeat(k, len(v)))
        vs += v
    return (ks, vs)


class dotdict(dict):
    """
    access dictioanry items via dot syntax
    via: http://parand.com/say/index.php/2008/10/24/python-dot-notation-dictionary-access/:
    """

    def __getattr__(self, attr):
        return self.get(attr, None)
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__


def nested_lookup(d, t):
    """
    look up if you can get a tuple of values from a nested dictionary,
    each item in the tuple a deeper layer
    """
    return reduce(lambda d, t: d.get(t, {}), t, d)

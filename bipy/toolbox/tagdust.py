"""
tagdust -fdr fdr_rate -o clean_tags_file -a dirty_tags_file
-singleline -quiet mode
"""
from bipy.utils import flatten, flatten_options, append_stem
import subprocess
import os
from bcbio.utils import safe_makedir


def _build_output_file(input_file, suffix, config):
    base = os.path.basename(input_file)
    return os.path.join(config["dir"]["results"], "tagdust",
                        append_stem(base, suffix))


def _build_output_files(input_file, tagdust_config, config):
    return [_build_output_file(input_file, x, config)
            for x in tagdust_config["keep"]]


def _build_command(input_file, tagdust_config, config):
    cl = [tagdust_config["program"], flatten_options(tagdust_config)]
    if "clean" in tagdust_config["keep"]:
        cl += ["-o", _build_output_file(input_file, "clean", config)]
    if "dirty" in tagdust_config["keep"]:
        cl += ["-a", _build_output_file(input_file, "dirty", config)]
    cl += [tagdust_config["contaminants"], input_file]
    return list(map(str, flatten(cl)))


def run(input_file, tagdust_config, config):
    cl = _build_command(input_file, tagdust_config, config)
    safe_makedir(os.path.join(config["dir"]["results"],
                              tagdust_config["name"]))
    output_files = list(_build_output_files(input_file, tagdust_config,
                                            config))
    # if all of the output files already exist, skip and exit
    if all([os.path.exists(x) for x in output_files]):
        return output_files

    subprocess.check_call(cl)
    return list(_build_output_files(input_file, tagdust_config, config))

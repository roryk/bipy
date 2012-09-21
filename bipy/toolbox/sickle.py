import subprocess
from bipy.utils import append_stem
import os


def run(in_file, end="se", qual="sanger", l="20", out_file=None):
    if not out_file:
        out_file = append_stem(in_file, "trimmed")

    if os.path.exists(out_file):
        return out_file

    cmd = ["sickle", end, "-f", in_file, "-o", out_file,
           "-t", qual, "-l", l]

    subprocess.check_call(cmd)
    return out_file

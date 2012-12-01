"""
wrapper around Trim galore! for trimming off common adapter
sequences used in NGS

http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

"""

from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, safe_makedir
from bipy.pipeline.stages import AbstractStage
from bipy.utils import flatten, append_stem, get_in
import sh
import os
from pkg_resources import resource_stream
import yaml
from Bio.Seq import Seq

with resource_stream(__name__, 'data/adapters.yaml') as in_handle:
    ADAPTERS = yaml.load(in_handle)


class TrimGalore(AbstractStage):
    """
    runs trim_galore on the data to trim off adapters, polyA tails and
    other contaminating sequences
    http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

    """

    stage = "trim_galore"

    def __init__(self, config):
        self.config = config
        self.stage_config = get_in(config, ("stage", self.stage), {})
        self.chemistry = self.stage_config.get("chemistry", "truseq")
        if not isinstance(self.chemistry, list):
            self.chemistry = [self.chemistry]

        self.options = self.stage_config.get("options", "")
        self.trim_galore = sh.Command(self.stage_config.get("program",
                                                            "trim_galore"))
        self.out_dir_prefix = get_in(config, ("dir", "results"),
                                     "results") + self.stage
        self.trim_polya = self.stage_config.get("trim_polya", True)
        safe_makedir(self.out_dir_prefix)

    def get_adapters(self, chemistry):
        return list(flatten([["-a", x] for x in ADAPTERS.get(chemistry, [])]))

    def _in2out(self, in_file):
        base, _ = os.path.splitext(in_file)
        return base + "_trimmed.fq"

    def __call__(self, in_file):
        raise NotImplementedError("Waiting to hear back from maintainer to "
                                  "handle multiple adapters before finishing.")
        adapters = list(flatten(map(self.get_adapters, self.chemistry)))
        # if it is a list assume these are pairs
        if isinstance(in_file, list):
            out_files = map(self._in2out, in_file)
            if all(map(file_exists, out_files)):
                return out_files
            self.trim_galore(in_file, self.options, adapters, paired=True)
            return out_files
        # if it is only one file just run it
        else:
            out_file = self._in2out(in_file)
            if file_exists(out_file):
                return out_file
            self.trim_galore(in_file, self.options, adapters)
            return out_file


class Cutadapt(AbstractStage):

    stage = "cutadapt"

    def __init__(self, config):
        self.config = config
        self.stage_config = get_in(config, ("stage", self.stage), {})
        self.chemistry = self.stage_config.get("chemistry", "truseq")
        if not isinstance(self.chemistry, list):
            self.chemistry = [self.chemistry]

        self.options = self.stage_config.get("options", "")
        self.cutadapt = sh.Command(self.stage_config.get("program",
                                                         "cutadapt"))
        self.out_dir = os.path.join(get_in(config, ("dir", "results"),
                                           "results"), self.stage)
        self.user_adapters = self.stage_config.get("adapters", [])
        self.trim_polya = self.stage_config.get("trim_polya", True)
        safe_makedir(self.out_dir)

    def in2out(self, in_file):
        """
        return expected out_file name if cutadapt is run on
        in_file

        """
        basename = os.path.basename(in_file)
        base, _ = os.path.splitext(basename)
        return os.path.join(self.out_dir, base + "_trimmed.fq")

    def _get_adapters(self, chemistry):
        adapters = [ADAPTERS.get(x, []) for x in chemistry]
        adapters += self.user_adapters
        adapters = list(flatten(adapters))
        adapters += self._rc_adapters(adapters)
        adapter_args = [["-a", adapter] for adapter in adapters]
        return list(flatten(adapter_args))

    def _rc_adapters(self, adapters):
        rc = [str(Seq(x).reverse_complement()) for x in adapters]
        return rc

    def __call__(self, in_file):
        adapters = self._get_adapters(self.chemistry)
        out_file = self.in2out(in_file)
        if file_exists(out_file):
            return out_file
        if self.trim_polya:
            cut_file = append_stem(out_file, "cut")
            # first cut off the adapters
            with file_transaction(cut_file) as tmp_cut_file:
                self.cutadapt(in_file, self.options, adapters,
                              _out=tmp_cut_file)
            # then trim the poly-a sequences
            with file_transaction(out_file) as tmp_out_file:
                polya = ADAPTERS.get("polya")
                self.cutadapt(cut_file, self.options, "-a",
                              polya, "-a", self._rc_adapters(polya),
                              out=tmp_out_file)
            # clean up the cut file
            os.remove(cut_file)

        else:
            # if not trimming polya, just trim the adadpters off
            with file_transaction(out_file) as tmp_out_file:
                self.cutadapt(in_file, self.options, adapters,
                              _out=cut_file)

        return out_file

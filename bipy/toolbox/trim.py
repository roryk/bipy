"""
wrappers around Trim galore! and cutadapt for trimming off common adapter
sequences used in NGS

trim galore!: http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
cutadapt: https://github.com/marcelm/cutadapt
sickle: https://github.com/najoshi/sickle
"""

from bcbio.distributed.transaction import file_transaction
from bcbio.utils import file_exists, safe_makedir
from bipy.pipeline.stages import AbstractStage
from bipy.utils import flatten, append_stem, get_in, is_pair
import sh
import os
from pkg_resources import resource_stream
import yaml
from Bio.Seq import Seq
import tempfile
from bipy.toolbox.fastq import DetectFastqFormat
from bipy.toolbox import fastq
from bcbio.provenance import do
from bcbio.log import logger, setup_local_logging

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
        return base + "_trimmed.fastq"

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

        self.options = self.stage_config.get("options", {})
        # it is important that we don't drop any reads out during cutadapt
        self.options["minimum-length"] = 0
        self.user_adapters = self.stage_config.get("adapters", [])
        self.out_dir = os.path.join(get_in(self.config, ("dir", "results"),
                                           "results"), self.stage)
        self.length_cutoff = self.stage_config.get("length_cutoff", 30)
        self.quality_format = self.stage_config.get("quality_format", None)

    def _detect_fastq_format(self, in_file):
        formats = DetectFastqFormat.run(in_file)
        sanger = ["sanger", "illumina_1.8+"]

        # if it is a sanger variant, 33 if not return the base 64
        if any([x in formats for x in sanger]):
            return "sanger"
        else:
            return "illumina"

    def in2trimmed(self, in_file):
        """
        return expected out_file name if cutadapt is run on
        in_file

        """
        basename = os.path.basename(in_file)
        base, _ = os.path.splitext(basename)
        safe_makedir(self.out_dir)
        return os.path.join(self.out_dir, base + "_trimmed.fastq")


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

    def _cut_file(self, in_file):
        """
        run cutadapt on a single file

        """
        adapters = self._get_adapters(self.chemistry)
        out_file = self.in2trimmed(in_file)
        if file_exists(out_file):
            return out_file
        cutadapt = sh.Command(self.stage_config.get("program",
                                                    "cutadapt"))

        quality_format = self.quality_format
        if not quality_format:
            quality_format = self._detect_fastq_format(in_file)
        if quality_format == "sanger":
            logger.info("Quality format detected as sanger.")
            quality_base = 33
        elif quality_format == "illumina":
            logger.info("Quality format set to illumina 1.5/1.3")
            quality_base = 64
        else:
            logger.error("Quality format could not be detected. Quality "
                         "Detected or set as %s. It should be illumina "
                         "or sanger.")
            exit(1)

        # if we want to trim the polya tails we have to first remove
        # the adapters and then trim the tail
        if self.stage_config.get("trim_polya", True):
            temp_cut = tempfile.NamedTemporaryFile(suffix=".fastq",
                                                   dir=self.out_dir)
            # trim off adapters
            cmd = str(cutadapt.bake(in_file, self.options, adapters,
                                    quality_base=quality_base, out=temp_cut.name))
            do.run(cmd, "Cutadapt trim of adapters of %s." % (in_file), None)
            with file_transaction(out_file) as temp_out:
                polya = ADAPTERS.get("polya")
                # trim off polya
                cmd = str(cutadapt.bake(temp_cut.name, self.options, "-a",
                                        polya, "-a", self._rc_adapters(polya),
                                        quality_base=quality_base, out=temp_out))
                do.run(cmd, "Cutadapt trim of polyA tail of %s." % (temp_cut.name),
                       None)
            return out_file
        else:
            with file_transaction(out_file) as temp_out:
                cmd = str(cutadapt.bake(in_file, self.options, adapters,
                                    out=temp_out))
                do.run(cmd, "Cutadapt trim of %s." % (in_file))
            return out_file

    def _get_lf_file(self, in_file):
        base, ext = os.path.splitext(in_file)
        out_file = base + ".fixed" + ext
        return out_file

    def _run_se(self, in_file):
        # cut polyA tails and adapters off
        logger.info("Running cutadapt in single end mode on %s." % (in_file))
        trimmed_file = self._cut_file(in_file)
        out_file = self._get_lf_file(trimmed_file)
        if file_exists(out_file):
            return out_file
        fastq.filter_single_reads_by_length(trimmed_file,
                                            self.length_cutoff)

        return out_file

    def _run_pe(self, in_files):
        logger.info("Running cutadapt in paired end mode on %s." % (in_files))
        trimmed_files = map(self._cut_file, in_files)
        out_files = map(self._get_lf_file, trimmed_files)
        if all(map(file_exists, out_files)):
            return out_files
        fastq.filter_reads_by_length(trimmed_files[0], trimmed_files[1],
                                     self.length_cutoff)

        return out_files

    def __call__(self, in_file):
        if is_pair(in_file):
            out_files = self._run_pe(in_file)
            return out_files
        elif isinstance(in_file, str):
            out_file = self._run_se(in_file)
            return out_file
        elif len(in_file) == 1:
            out_file = self._run_se(in_file[0])
            return out_file
        else:
            raise ValueError("Cutadapt can only run on either a single "
                             "file as a string or a pair of files to "
                             "handle paired end data.")

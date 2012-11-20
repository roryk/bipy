from bipy.toolbox import fastqc
from bipy.log import logger
from bcbio.utils import file_exists
import sh
import os
import csv
from bcbio.variation import effects


class AbstractStage(object):
    """
    Wrapper around a stage for use with the pipeline manager. Subclasses should
    implement their own __init__ in addition to calling the super
    __init__. __init__ should parse the appropriate information from the config
    dict and save it in the objects state. __call__ should be implemented as a
    pure function with no closures.

    """

    stage = "abstract"

    def __init__(self, config):
        self.config = config
        self._validate_config()

    def _start_message(self, in_file, **kwargs):
        if kwargs:
            logger.info("Starting %s on %s with arguments %s." % (self.stage,
                                                                  in_file,
                                                                  kwargs))
        else:
            logger.info("Starting %s on %s." % (self.stage, in_file))

    def _end_message(self, in_file):
        logger.info("%s complete on %s." % (self.stage, in_file))

    def _check_run(self, in_file):
        if not file_exists(in_file):
            raise IOError("%s not found." % (in_file))

    def __call__(self, in_file):
        pass

    def _validate_config(self):
        if "stage" not in self.config:
            raise ValueError('Could not find "stage" in the config file.')


class FastQCStage(AbstractStage):

    stage = "fastqc"

    def __init__(self, config):
        self.config = config
        super(FastQCStage, self).__init__(self.config)
        self.stage_config = config["stage"][self.stage]

    def _start_message(self, in_file):
        logger.info("Starting %s on %s" % (self.stage, in_file))

    def _end_message(self, in_file):
        logger.info("%s complete on %s." % (self.stage, in_file))

    def _check_run(self, in_file):
        if not file_exists(in_file):
            raise IOError('%s not found.' % (in_file))

    def __call__(self, in_file):
        self._start_message(in_file)
        self._check_run(in_file)
        out_file = fastqc.run(in_file, self.stage_config, self.config)
        self._end_message(in_file)
        return out_file


class IlluminaVCFFixer(AbstractStage):

    stage = "illumina_fixer"
    variation_version = "bcbio.variation-0.0.6-SNAPSHOT-standalone.jar"

    def __init__(self, config):
        self.config = config
        super(IlluminaVCFFixer, self).__init__(self.config)
        self.grc_file = config["ref"]["grc_file"]
        self.ucsc_file = config["ref"]["ucsc_file"]
        self.fixer = os.path.join(config["program"]["bcbio.variation"],
                                  self.variation_version)
        self.sample_lookup = self._lane2sample(self.config)

    def _lane2sample(self, config):
        """
        return dict to translate illumina sample ids to harvard ids
        """
        sample_field = "sample_file"
        sample_file = config.get(sample_field, None)

        if not sample_file or not file_exists(sample_file):
            raise ValueError("%s in annotate.yaml needs to point to "
                             "a CSV file of the samples to look up."
                             % (sample_field))

        with open(sample_file, 'rb') as in_handle:
            reader = csv.reader(in_handle, delimiter=",")
            d = {}
            for line in reader:
                d[line[0]] = line[1]

        return d

    def _run_fixer(self, vcf_dir):
        sample = self.sample_lookup.get(os.path.basename(vcf_dir), None)
        if not sample:
            raise ValueError("Could not find lane %s in the lookup table "
                             "%s" % (self.sample_lookup))

        out_file = os.path.join(vcf_dir, sample + ".vcf")
        if file_exists(out_file):
            return out_file

        sh.java("-jar", self.fixer, "variant-utils", "illumina", vcf_dir,
                sample, self.grc_file, self.ucsc_file)

        return out_file

    def __call__(self, vcf_dir):
        self._start_message(vcf_dir)
        out_file = self._run_fixer(vcf_dir)
        self._end_message(vcf_dir)
        return out_file


class SnpEff(AbstractStage):

    stage = "snpeff"

    def __init__(self, config):
        self.config = config
        super(SnpEff, self).__init__(self.config)
        self.genome = config["ref"]["name"]

    def __call__(self, in_file):
        self._start_message(in_file)
        out_file = effects.snpeff_effects(in_file, self.genome, self.config)
        self._end_message(in_file)
        return out_file


class GeminiLoader(AbstractStage):

    stage = "geminiloader"

    def __init__(self, config):
        self.config = config
        gemini_stage = self.config["stage"]["gemini"]
        super(GeminiLoader, self).__init__(self.config)
        self.gemini = config["program"].get("gemini", "gemini")
        self.type = gemini_stage.get("type", "snpEff")
        self.db = os.path.join(config["dir"].get("results", "results"),
                               self.stage,
                               gemini_stage.get("db", "combined.db"))

    def _load_gemini(self, in_file):
        sh.gemini.load(self.db, v=in_file, t=self.type)

    def __call__(self, in_file):
        self._start_message(in_file)
        self._load_gemini(in_file)
        self._end_message(in_file)
        return self.db


STAGE_LOOKUP = {"fastqc": FastQCStage,
                "geminiloader": GeminiLoader,
                "snpeff": SnpEff,
                "illumina_fixer": IlluminaVCFFixer}

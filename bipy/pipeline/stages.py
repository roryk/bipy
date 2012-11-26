from bipy.toolbox import fastqc
from bipy.log import logger
from bcbio.utils import file_exists, safe_makedir
from bipy.utils import append_stem, replace_suffix
import sh
import os
import csv
from bcbio.variation import effects
from bcbio.distributed.transaction import file_transaction
import uuid


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
        self.java_memory = "-Xmx" + self.config["algorithm"].get("java_memory",
                                                                 "1g")

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

        out_file = os.path.join(vcf_dir, "Variations", sample + ".vcf")
        if file_exists(out_file):
            return out_file

        sh.java(self.java_memory,
                "-jar", self.fixer, "variant-utils", "illumina", vcf_dir,
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


class Vep(AbstractStage):

    stage = "vep"

    def __init__(self, config):
        self.config = config
        super(Vep, self).__init__(self.config)
        self.vep_config = config["stage"].get("vep", {})
        self.species = self.vep_config.get("species", "human")
        self.options = self.vep_config.get("options", None)
        self.vep = self.config["program"].get("vep",
                                              "variant_effect_predictor.pl")

    def __call__(self, in_file):
        self._start_message(in_file)
        out_file = self._run_vep(in_file)
        self._end_message(in_file)
        return out_file

    def _run_vep(self, in_file):
        out_file = append_stem(in_file, "vep")
        if file_exists(out_file):
            return out_file

        with file_transaction(out_file) as tmp_out_file:
            sh.perl(self.vep, "-i", in_file, "-o", tmp_out_file,
                    species=self.species, _convert_underscore=False,
                    **self.options)

        return out_file


class GeminiLoader(AbstractStage):

    stage = "geminiloader"

    def __init__(self, config):
        self.config = config
        gemini_stage = self.config["stage"]["gemini"]
        super(GeminiLoader, self).__init__(self.config)
        self.gemini = config["program"].get("gemini", "gemini")
        self.type = gemini_stage.get("type", "snpEff")
        out_dir = os.path.join(config["dir"].get("results", "results"),
                               self.stage)
        safe_makedir(out_dir)
        self.db = os.path.join(out_dir,
                               gemini_stage.get("db", "combined.db"))
        self.log = config["dir"].get("log", "log")

    def _load_gemini(self, in_file):
        log_id = os.path.join(self.log,
                              "gemini" + "_" + uuid.uuid4() + ".log")
        sh.gemini.load(self.db, v=in_file, t=self.type,
                       _out=append_stem(log_id, "out"),
                       _err=append_stem(log_id, "err"))


    def __call__(self, in_file):
        self._start_message(in_file)
        self._load_gemini(in_file)
        self._end_message(in_file)
        return self.db


class BreakVcfByChromosome(AbstractStage):

    stage = "break_vcf_by_chromosome"

    def __init__(self, config):
        self.config = config
        self.fasta_file = self.config["ref"]["fasta"]
        self.fasta_index = self.fasta_file + ".fai"
        self.tabix = config["program"].get("tabix", "tabix")
        self.samtools = config["program"].get("samtools", "samtools")

    def _break_vcf(self, in_file):
        if not file_exists(self.fasta_index):
            sh.samtools.faidx(self.fasta_file)

        # if file is not compressed, compress it
        (_, ext) = os.path.splitext(in_file)
        if ext is not ".gz":
            gzip_file = in_file + ".gz"
            sh.bgzip("-c", in_file, _out=gzip_file)
            in_file = gzip_file

        # create tabix index if it does not exist already
        if not file_exists(in_file + ".tbi"):
            sh.tabix("-p", "vcf", in_file)

        # find the chromosome names from the fasta index file
        chroms = str(sh.cut("-f1", self.fasta_index)).split()
        break_dir = os.path.join(os.path.dirname(in_file), "break")
        safe_makedir(break_dir)

        def chr_out(chrom):
            out_file = os.path.join(break_dir, append_stem(in_file, chrom))
            out_file = replace_suffix(out_file, "vcf")
            return out_file

        def tabix(chrom):
            out_file = chr_out(chrom)
            if file_exists(out_file):
                return out_file
            with file_transaction(out_file) as tmp_out_file:
                sh.tabix("-h", in_file, chrom, _out=tmp_out_file)
            return out_file

        # use tabix to separate out the variants based on chromosome
        out_files = map(tabix, chroms)

        return out_files

    def __call__(self, in_file):
        self._start_message(in_file)
        out_files = self._break_vcf(in_file)
        self._end_message(in_file)
        return out_files


STAGE_LOOKUP = {"fastqc": FastQCStage,
                "geminiloader": GeminiLoader,
                "snpeff": SnpEff,
                "illumina_fixer": IlluminaVCFFixer,
                "vep": Vep,
                "breakvcf": BreakVcfByChromosome}

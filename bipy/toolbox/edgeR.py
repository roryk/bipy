"""
wrapper to run EdgeR analysis in R:
http://www.bioconductor.org/packages/release/bioc/html/edgeR.html

"""
from rpy2.robjects.vectors import DataFrame, FactorVector, StrVector
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from bipy.pipeline.stages import AbstractStage
from bipy.utils import remove_suffix
from bcbio.utils import file_exists


class EdgeR(AbstractStage):

    def __init__(self, config):
        self.config = config
        self.r = robjects.r

    def _load_csv_file(self, in_file, r):
        r.assign('in_file', in_file)
        r('''
        counts = read.table(in_file, header=TRUE, row.names=1)
        ''')
        return r

    def _make_count_set(self, in_file, conds, r):
        self._load_csv_file(in_file, r)
        r.assign('group', StrVector.factor(StrVector(conds)))
        r('''
        require('edgeR')
        cds = edgeR.DGEList(counts=counts, group=group)
        ''')
        return r

    def _make_mds_plot(self, in_file, r):
        out_file = remove_suffix(in_file) + "_MDSplot.pdf"
        if file_exists(out_file):
            return out_file
        r.assign("plot_name", out_file)
        r('''
        pdf(plot_name, width=7, height=7)
        plotMDS.dge(cds, main="MDS Plot for Count Data",
        labels = colnames(cds$counts))
        dev.off()
        ''')
        return out_file

    def _estimate_dispersions(self, r, prior_n):
        r.assign("prior.n", prior_n)
        r('''
        cds = estimateCommonDisp(cds)
        cds = estimateTagwiseDisp(cds, prior.n)
        ''')
        return r

    def _make_mva_plot(self, in_file, r):
        out_file = remove_suffix(in_file) + "_MvAplot.pdf"
        if file_exists(out_file):
            return out_file
        r.assign("mva_plot_name", out_file)
        r('''
        pdf(mva_plot_name, width=7, height=7)
        plotMeanVar(cds, show.raw.vars=TRUE, show.tagwise.vars=TRUE,
        show.binned.common.disp.vars=FALSE, show.ave.raw.vars=FALSE,
        dispersion.method="qcml", NBline=TRUE,
        nbins=100,
        pch=16,
        xlab="Mean Expression (log10 scale)",
        ylab="Variance (log10 scale)",
        main="Mean-Variance plot")
        ''')
        return out_file

    def __call__(self, in_file, conds):
        r = robjects.r
        r = self._make_count_set(in_file, conds, r)
        mds_plot = self._make_mds_plot(in_file, r)
        prior_n = 50 / (len(conds) - len(Set(conds)))
        r = self._estimate_dispersions(r, prior_n)
        mva_plot = self._make_mva_plot(in_file, r)

"""
targets is a pandas dataframe
convert it to a robject dataframe
subject treatment
design = model.matrix(~subject + treat, targets)
fit = glmFit(y, design)
lrt = glmLRT(fit)
topTags(fit)



# simple example sequence
count_table_file = "/Users/rory/cache/bipy/test/data/pasilla_gene_counts.tsv"
count_table_file = "xxx"
edgeR = importr("edgeR")
base = importr("base")
utils = importr("utils")
groups = ["u", "u", "u", "u", "t", "t", "t"]
grp = FactorVector(groups)
#counts = utils.read_table(count_table_file, sep="\t", header=True,
#                          row_names = 1)
counts = DataFrame.from_csvfile(count_table_file, sep="\t", header=True,
                                row_names = 1)

cds = edgeR.DGEList(counts=counts, group=grp)"""

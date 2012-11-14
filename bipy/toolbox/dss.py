"""
wrapper to perform a two condition DE testing of RNA-seq data using
the DSS package. requires in_file to be a table of counts and conds to
be a vector describing which condition each column of counts refers to
"""
import rpy2.robjects as robjects
import rpy2.robjects.vectors as vectors
from bcbio.utils import safe_makedir
import os


def _plot_disp_ests(r, dispersion_plot_out):
    """
    make a plot of the dispersion estimation
    XXX: this needs to get ported to work with DSS
    """
    r.assign("dispersion_plot_out", dispersion_plot_out)
    r('''
        plotDispEsts <- function(cds) {
          plot(rowMeans(counts(cds, normalized=TRUE)),
          fitInfo(cds)$perGeneDispEsts,
          pch = '.', log="xy")
          xg <- 10^seq(-.5,5,length.out=300)
          lines(xg,fitInfo(cds)$dispFun(xg),col="red")}
        pdf(dispersion_plot_out)
        plotDispEsts(cds)
        dev.off()
          ''')
    return r


def _plot_MvA(r, mva_plot_out):
    """
    make a MvA plot of the differential expression
    XXX: calculate basemean for dss
    """
    r.assign("mva_plot_out", mva_plot_out)
    r('''
    plotMvA <- function(res) {
      plot(res$baseMean,res$lfc,log="x",pch=20,cex=.3,
      col=ifelse(res$padj < .1, "red","black")) }
    pdf(mva_plot_out)
    plotMvA(res)
    dev.off()
    ''')
    return r


def load_count_file_as_matrix(in_file, r):
    """
    load in counts from a table and convert into matrix form for
    DSS
    """
    r.assign('in_file', in_file)
    r('''
    count_table = read.table(in_file, header=TRUE, row.names=1)
    count_matrix = as.matrix(count_table)
    dimnames(count_matrix) = NULL
    ''')
    return r


def make_count_set(conds, r):
    """
    returns an r session with a new count data set loaded as cds
    """
    #r.assign('conds', vectors.StrVector.factor(vectors.StrVector(conds)))
    r.assign('conds', vectors.StrVector(conds))
    r('''
    require('DSS')
    cds = newSeqCountSet(count_matrix, as.character(conds))
    ''')
    return r


def run(in_file, conds, out_prefix):
    dss_table_out = out_prefix + ".dss.txt"

    safe_makedir(os.path.dirname(out_prefix))
    r = robjects.r
    r.assign('dss_table_out', dss_table_out)

    r = load_count_file_as_matrix(in_file, r)
    r = make_count_set(conds, r)
    r('''
    cds = estNormFactors(cds)
    cds = estDispersion(cds)
    conditions = levels(factor(conds))
    res = waldTest(cds, conditions[0], conditions[1])
    write.table(res, file=dss_table_out, quote=FALSE,
    row.names=FALSE, sep="\t")
    ''')

    return dss_table_out

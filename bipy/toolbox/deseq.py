import os
import rpy2.robjects as robjects
import rpy2.robjects.vectors as vectors
import pandas as pd
from bcbio.utils import safe_makedir

this_dir, this_filename = os.path.split(__file__)
RFILE = os.path.join(this_dir, "R", "diffexp.R")


def run(in_file, conds, out_prefix):
    deseq_table_out = out_prefix + ".deseq.txt"
    dispersion_plot_out = out_prefix + ".dispersions.pdf"
    mva_plot_out = out_prefix + ".MvA.pdf"

    safe_makedir(os.path.dirname(out_prefix))

    r = robjects.r
    r.assign('in_file', in_file)
    r.assign('deseq_table_out', deseq_table_out)
    r.assign('mva_plot_out', mva_plot_out)
    r.assign('conds',
             vectors.StrVector.factor(vectors.StrVector(conds)))
    r('''
    require('DESeq')
    counts = read.table(in_file, header=TRUE, row.names=1)
    print(head(counts))
    cds = newCountDataSet(counts, conds)
    print(head(cds))
    cds = estimateSizeFactors(cds)
    ''')

    # if there are no replicates, use the replicate cheating mode
    if len(set(conds)) == len(conds):
        r('''
        cds = estimateDispersions(cds, method="blind", sharingMode="fit-only")
        ''')
    else:
        r('''
        cds = estimateDispersions(cds)
        ''')

    _plot_disp_ests(r, dispersion_plot_out)

    # if there are two conditions use the standard deseq diffexpression
    sconds = set(conds)
    if len(sconds) == 2:
        r.assign('cond1', str(list(sconds)[0]))
        r.assign('cond2', str(list(sconds)[1]))
        r('''
        res = nbinomTest(cds, cond1, cond2)
        ''')
        _plot_MvA(r, mva_plot_out)
    r('''write.table(res, file=deseq_table_out, quote=FALSE,
    row.names=FALSE, sep="\t")''')
    return deseq_table_out


def _plot_disp_ests(r, dispersion_plot_out):
    """
    make a plot of the dispersion estimation
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
    """
    r.assign("mva_plot_out", mva_plot_out)
    r('''
    plotMvA <- function(res) {
      plot(res$baseMean,res$log2FoldChange,log="x",pch=20,cex=.3,
      col=ifelse(res$padj < .1, "red","black")) }
    pdf(mva_plot_out)
    plotMvA(res)
    dev.off()
    ''')
    return r


def run_with_config(in_file, gtf, config):
    pass

import os
import rpy2.robjects as robjects
import rpy2.robjects.vectors as vectors
import pandas as pd
from bcbio.utils import safe_makedir

this_dir, this_filename = os.path.split(__file__)
RFILE = os.path.join(this_dir, "R", "diffexp.R")

""" what does deseq need. a counts table of gene name and condition. i
think i have that via the htseq-count script, just not combined. maybe
need to take those in and then merge them
step 1: merge the htseq-count files
step 2: load in the table
step 3: load in the experimental design table
row names are the conditions in the count table
need condition the experiment condition
plus other conditions as columns

how to figure out the conditions
filename: condition i guess makes the most sense
or just a list corresponding to the filenames

then we take those counts and merge them into one file
load the conds
call cds <- newCountDataSet with and that and the conds
call estimateSizeFactors(cds)
call estimateDispersions(cds)
# do those diagnostic plots after this works

now need to determine which type of analysis to do
it can either take a parameter that says which type
or figure it out based on the number of conditions

## handling DE for two treatments
if only two and each has more than one:
res <- nbinomTest( cds, "untreated", "treated" )
# make plotDe function
# make pval hist function
# write table out

# if you only have replicates for one condition
this is the same as above

# if you have no replicates
estimateDispersions(cds, method="blind", sharingMode="fit-only")
res <- nbinomTest( cds, "untreated", "treated" )

## handling DE for multiple treatments
this is a little more involved, keep this as an option for now
for now output a not implemented message and fail gracefully
> fit1 <- fitNbinomGLMs( cdsFull, count ~ libType + condition )
> fit0 <- fitNbinomGLMs( cdsFull, count ~ libType  )

 cds <- newCountDataSet( countsTable, conds )
 cds <- estimateSizeFactors( cds )
 cds <- estimateDispersions(cds)
"""


def run(in_file, conds, out_prefix):
    deseq_table_out = out_prefix + ".deseq.txt"
    dispersion_plot_out = out_prefix + ".dispersions.pdf"
    mva_plot_out = out_prefix + ".MvA.pdf"

    safe_makedir(os.path.dirname(out_prefix))

    r = robjects.r
    r.assign('in_file', in_file)
    r.assign('deseq_table_out', deseq_table_out)
    r.assign('dispersion_plot_out', dispersion_plot_out)
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

    # if there are two conditions use the standard deseq diffexpression
    sconds = set(conds)
    if len(sconds) == 2:
        r.assign('cond1', str(list(sconds)[0]))
        r.assign('cond2', str(list(sconds)[1]))
        r('''
        res = nbinomTest(cds, cond1, cond2)
        ''')
        r('''
        plotDE <- function(res) {
          plot(res$baseMean,res$log2FoldChange,log="x",pch=20,cex=.3,
          col=ifelse(res$padj < .1, "red","black")) }
        pdf(mva_plot_out)
        plotDE(res)
        dev.off()
        ''')
    r('''write.table(res, file=deseq_table_out, quote=FALSE,
    row.names=FALSE, sep="\t")''')
    return deseq_table_out


def run_with_config(in_file, gtf, config):
    pass

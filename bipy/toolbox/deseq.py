"""
wrapper to run DESeq analyses. at the moment it will only run with a simple
two-condition comparison. more complicated experimental setups could be
added by creating a dataframe of the conditions and passing it in
"""
import os
import rpy2.robjects as robjects
import rpy2.robjects.vectors as vectors
import pandas as pd
from bcbio.utils import safe_makedir, file_exists
from bipy.toolbox.reporting import LatexReport, panda_to_latex
from mako.template import Template



def load_count_file(in_file, r):
    """
    returns an r session with a read count file loaded as 'counts'
    """
    r.assign('in_file', in_file)
    r('''
    counts = read.table(in_file, header=TRUE, row.names=1)
    ''')
    return r


def make_count_set(conds, r):
    """
    returns an r session with a new count data set loaded as cds
    """
    r.assign('conds', vectors.StrVector.factor(vectors.StrVector(conds)))
    r('''
    require('DESeq')
    cds = newCountDataSet(counts, conds)
    ''')
    return r


def run(in_file, conds, out_prefix):
    deseq_table_out = out_prefix + ".deseq.txt"
    dispersion_plot_out = out_prefix + ".dispersions.pdf"
    mva_plot_out = out_prefix + ".MvA.pdf"

    safe_makedir(os.path.dirname(out_prefix))

    r = robjects.r
    r.assign('deseq_table_out', deseq_table_out)
    r.assign('mva_plot_out', mva_plot_out)

    r = load_count_file(in_file, r)
    r = make_count_set(conds, r)
    r('''
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


class DeseqParser(object):
    """
    parse a directory of a deseq experiment to generate a summary report

    """
    GRAPH_SUFFIXES = ((".dispersions.pdf", "", 1.0),
                      (".MvA.pdf", "", 1.0))

    INFO_SUFFIXES = {"annotated": ".annotated.deseq.txt",
                     "deseq": ".deseq.txt"}

    def __init__(self, base_dir):
        self._dir = base_dir
        self._comparison = os.path.basename(base_dir)
        self._top_max = 25

    def get_deseq_graphs(self):
        final_graphs = []
        for suffix, caption, size in self.GRAPH_SUFFIXES:
            final_f = os.path.join(self._dir, self._comparison + suffix)
            if file_exists(final_f):
                final_graphs.append((final_f, caption, size))
        return final_graphs

    def get_top_genes(self):
        final_info = {}
        for info, suffix in self.INFO_SUFFIXES.items():
            final_f = os.path.join(self._dir, self._comparison + suffix)
            if file_exists(final_f):
                final_info[info] = final_f

        # prefer to use the annotated file
        if "annotated" in final_info:
            top_table = self._get_top_by(final_info["annotated"])
        else:
            top_table = self._get_top_by(final_info["deseq"])

        return (top_table, "")

    def _get_top_by(self, in_file, by="padj"):
        if "annotated" in in_file:
            cols = ["id", "log2FoldChange", "padj", "symbol"]
        else:
            cols = ["id", "log2FoldChange", "padj"]
        if by not in cols:
            by = "padj"
        table = pd.read_table(in_file, header=0)
        # grab the top genes up to top_max genes
        sub_table = table[cols].sort(columns=by)[0:self._top_max]

        return sub_table


class DeseqReport(LatexReport):

    def __init__(self):
        pass

    def template(self):
        return self._template

    def generate_report(self, name, figures=None, top_hits=None):
        template = Template(self._template)
        section = template.render(name=name, figures=figures,
                                  top_hits=top_hits)
        return section

    _template = r"""
\subsection*{DESeq report for ${name}}

% if figures:
    % for i, (figure, caption, size) in enumerate(figures):
        \begin{figure}[htbp]
            \centering
            \includegraphics[width=${size}\linewidth] {${figure}}
            \caption{${caption}}
        \end{figure}
    % endfor
% endif

% if top_hits:
    <%
        from bipy.toolbox.reporting import panda_to_latex
        table, caption = top_hits
        table_out = panda_to_latex(table, caption)
    %>
    ${table_out}
% endif
"""


#    ${panda_to_latex(table, caption)}
#    ${panda_to_latex(top_hits[0], top_hits[1])}

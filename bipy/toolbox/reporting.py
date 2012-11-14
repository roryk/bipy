"""reporting.py provides classes to inherit from to provide output
reports for stages of pipelines.

"""

from mako.template import Template
import sh
import os
import abc
from itertools import repeat
from bipy.utils import remove_suffix
from bcbio.distributed.transaction import file_transaction
from bcbio.utils import curdir_tmpdir
import shutil


def safe_latex(to_fix):
    """Escape characters that make LaTeX unhappy.
    Lifted from Brad Chapman (bcbio)
    """
    chars = ["%", "_", "&", "#"]
    for char in chars:
        to_fix = to_fix.replace(char, "\\%s" % char)
    return to_fix


class LatexPdf(object):
    """Generates a PDF document from a list of sections.
    usage: generate_pdf(sections)"""

    @classmethod
    def generate_pdf(self, sections=None, out_file=None):
        out_tmpl = Template(self._base_template)
        if not out_file:
            latex_file = "latex.tex"
            out_file = "latex.pdf"
        else:
            latex_file = remove_suffix(out_file) + ".tex"

        with open(latex_file, "w") as latex_handle:
            latex_handle.write(out_tmpl.render(sections=sections))
        sh.pdflatex(latex_file)

        return out_file

        #        sh.pdflatex(latex_file,
        #with file_transaction(latex_file) as tmp_latex_file:
        #     with open(tmp_latex_file, "w") as latex_handle:
        #        latex_handle.write(out_tmpl.render(sections=sections))
        #        sh.pdflatex(tmp_
        #return out_file

    _base_template = r"""
\documentclass{article}
\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{placeins}
\usepackage[multidot]{grffile}

\begin{document}
% for section in sections:
    ${section}
% endfor
\end{document}
"""


class LatexReport(object):
    """Abstract class for generating Latex reports. Inherit from this
    and implement template and generate_report"""
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def template(self):
        return

    @abc.abstractmethod
    def generate_report(self, *args, **kwargs):
        """Generate the Latex output using the template"""
        return


def make_latex_table_header(length):
    tabular = "|l|" + "".join(repeat("r|", length - 1))
    table_header = r"\begin{tabular}{%s}" % (tabular)
    return table_header


def panda_to_latex(table, caption=""):
    header = make_latex_table_header(len(table.keys()))

    def _panda_row_to_latex_row(row):
        vals = map(str, [e for e in row[1]])
        return " & ".join(vals) + " \\\\"

    rows = [_panda_row_to_latex_row(x) for x in table.iterrows()]

    _table_template = r"""
\begin{table}[h]
    \centering
    ${header}

    % for row in rows:
        ${row}
    % endfor
    \hline
    \end{tabular}
    \caption{${caption}}
\end{table}
"""

    template = Template(_table_template)
    latex_table = template.render(header=header, rows=rows, caption=caption)
    return latex_table

"""reporting.py provides classes to inherit from to provide output
reports for stages of pipelines.

"""

from mako.template import Template
import sh
import os
import abc
from itertools import repeat


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
            out_file = "latex.tex"
        with open(out_file, "w") as out_handle:
            out_handle.write(out_tmpl.render(sections=sections))
        sh.pdflatex(out_file)
        return "%s.pdf" % os.path.splitext(out_file)[0]

    _base_template = r"""
\documentclass{article}
\usepackage{fullpage}
\usepackage{graphicx}
\usepackage{placeins}

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

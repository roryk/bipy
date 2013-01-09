from bipy.utils import append_stem
from rpy2 import robjects
import os
from bipy.log import logger

ORG_TO_ENSEMBL = {"opossum": {"gene_ensembl": "mdomestica_gene_ensembl",
                              "gene_symbol": "hgnc_symbol"},
                  "mouse": {"gene_ensembl": "mmusculus_gene_ensembl",
                            "gene_symbol": "mgi_symbol"},
                  "human": {"gene_ensembl": "hsapiens_gene_ensembl",
                            "gene_symbol": "hgnc_symbol"},
                  "taz": {"gene_ensembl": "sharrisii_gene_ensembl",
                          "gene_symbol": "hgnc_symbol"},
                  "zebrafish": {"gene_ensembl": "drerio_gene_ensembl",
                                "gene_symbol": "zfin_symbol"},
                    "c.elegans": {"gene_ensembl": "celegans_gene_ensembl",
                                  "filter_type": "wormbase_locus",
                                  "gene_symbol": "external_gene_id"}
                                }


def annotate_table_with_biomart(in_file, join_column,
                                filter_type, organism, out_file=None):
    """
    join_column is the column to combine the perform the lookups on
    filter_type describes the type of the join_column (see the getBM
    documentation in R for details), organism is the english name of
    the organism

    example:
    annotate_table_with_biomart(in_file, "id", "ensembl_gene_id",
                                "human")

    """

    if organism not in ORG_TO_ENSEMBL:
        logger.error("organism not supported")
        exit(1)

    logger.info("Annotating %s." % (organism))
    if not out_file:
        out_file = append_stem(in_file, "annotated")
    if os.path.exists(out_file):
        return out_file
    # use biomaRt to annotate the data file
    r = robjects.r
    r.assign('join_column', join_column)
    r.assign('in_file', in_file)
    r.assign('out_file', out_file)
    r.assign('ensembl_gene', ORG_TO_ENSEMBL[organism]["gene_ensembl"])
    r.assign('gene_symbol', ORG_TO_ENSEMBL[organism]["gene_symbol"])
    r.assign('filter_type', filter_type)
    r('''
    library(biomaRt)
    ensembl = useMart("ensembl", dataset = ensembl_gene)
    d = read.table(in_file, header=TRUE)
    a = getBM(attributes=c(filter_type,
                gene_symbol, "description"),
                filters=c(filter_type), values=d[,join_column],
                mart=ensembl)
    m = merge(d, a, by.x=join_column, by.y=filter_type)
    write.table(m, out_file, quote=FALSE, row.names=FALSE, sep="\t")
    ''')

    return out_file

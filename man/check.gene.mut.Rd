% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r-functions.selection.R
\name{check.gene.mut}
\alias{check.gene.mut}
\title{Retreve somatic mutations of a gene from one or more tumor types.}
\usage{
check.gene.mut(gene.symbol, folder, prefix)
}
\arguments{
\item{gene.symbol:}{A gene symbol.}

\item{folder:}{A folder containing a "prefix.mut.summary.txt" file and a "prefix.symbol_2_cds_id.txt" file.}

\item{prefix:}{A prefix string used to label the files. Usually, it contains TCGA.tumor_abbreviation.}
}
\value{
A list with three elements: the first element is a data frame with mutations counts aggregated on protein level; the second element is a data frame with mutations counts aggregated on nucleotide level; the third element is the Ensembl transcript ID.
}
\description{
Retreve somatic mutations of a gene from one or more tumor types.
}

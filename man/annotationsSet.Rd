% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotationSet.R
\name{annotationsSet}
\alias{annotationsSet}
\title{Create a geneAnnotation object.}
\usage{
annotationsSet(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "scerevisiae_gene_ensembl", host = "www.ensembl.org")
}
\arguments{
\item{biomart}{BioMart database name you want to connect to.}

\item{dataset}{Dataset you want to use.}

\item{host}{Host to connect to.}
}
\value{
An object of class geneAnnotation with the following components:
       \item{chromosome}{Gene chromosome.}
       \item{start}{Gene transcription start site.}
       \item{end}{Gene transcription termination site.}
       \item{name}{Gene name.}
       \item{score}{Gene score.}
       \item{strand}{Gene strand.}
}
\description{
Create a geneAnnotation object containing the informations on
chromosomes, start and end positions, strand, score and name of genes.
}
\details{
Please see \code{\link[biomaRt]{useMart}} for further information.
}
\examples{
annotationsSet()
}


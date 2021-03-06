% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{coverage}
\alias{coverage}
\title{Report coverage.}
\usage{
coverage(reads, annotation, sf = NULL, outfile = NULL)
}
\arguments{
\item{reads}{a \code{\link[data.table]{data.table}} of values
or a character vector indicating the input file name.
Note: the data must be in a BED-like format.}

\item{annotation}{a object of type \code{\link{annotationsSet}}
or a character vector indicating the input BED file name.}

\item{sf}{a number indicating the scaling factor.
Each coverage value is multiplied by this factor before being reported.
Useful for normalizing coverage by, e.g., reads per million (RPM).}

\item{outfile}{a character vector indicating the output file name.
If not provided, the result will be internalized in R.}
}
\value{
An object of type \code{\link[data.table]{data.table}} with 8 columns,
        corresponding to:
\enumerate{
   \item Chromosome
   \item Starting position
   \item Ending position
   \item Gene name
   \item Gene strand (+ or -)
   \item Positions. One based positions of the gene.
   \item Depth at each position of the gene.
}
}
\description{
Report the coverage using \href{http://bedtools.readthedocs.org/en/latest/content/tools/coverage.html}{bedtools coverage} algorithm.
}
\examples{
data(annotations)
bam  <- system.file("extdata", "wt_rep1.bam", package = "yCrypticRNAs")
fragments <- bam_to_reads(bam, annotations)
coverage(fragments, annotations, sf = 0.069872847)
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initiation_points_methods.R
\name{initiation_sites}
\alias{initiation_sites}
\title{Cryptic transcription start sites (cTSS)}
\usage{
initiation_sites(name, bamfiles, types, annotations, introns = NULL,
  sf = replicate(length(bamfiles), 1), replicates = 200, percentage = 0.1,
  method = c("methodC_gaussian", "methodA", "methodB", "methodC", "methodD"),
  paired_end = TRUE, as_fragments = TRUE)
}
\arguments{
\item{name}{The name of the gene for which to calculate the cryptic
initiation site(s).}

\item{bamfiles}{a vector of characters indicating the BAM file paths.}

\item{types}{a vector of the same length as bamfiles indicating the type
of data in each file. Example: WT vs MUT or untreated vs treated".}

\item{annotations}{An object of type \code{\link{annotationsSet}}
containing information on genes.}

\item{introns}{an objet of type \code{\link{annotationsSet}}
containing the annotations of the intronic regions.
Note: The introns must have same name as the gene they
are associated with.}

\item{sf}{A vector of the scaling factors to apply to each
sample. Must be the same length as the \code{fragments_file}}

\item{replicates}{The number of time to sample the data. Default = 200.}

\item{percentage}{Th fraction of data to be removed at each simulation.
Default = 0.1.}

\item{method}{A character string or a vector specifying the method to use to
calculate the cryptic initiations sites. Must be one of "methodC_gaussian"
(default), "methodA", "methodB", "methodC" or "methodD".}

\item{paired_end}{logical indicating whether the \code{bamfile}
contains paired_end data.}

\item{as_fragments}{logical indicating if paired_end data must paired
and merged to form fragments.}
}
\value{
A list with the following components:
      \item{methodC_gaussian}{cTSS mean and sd values using the
                              method C (gaussian)}
       \item{methodA}{cryptic zones start and end
                       position using the method A}
       \item{methodB}{cryptic zones start and end
                       position using the method B}
       \item{methodC}{cryptic zones start and end
                       position using the method C}
       \item{methodD}{cryptic zones start and end
                       position using the method D}
       \item{gene_information}{An object of class \code{\link{annotationsSet}}
                             containing the information on the gene.}
}
\description{
For genes expresing a cryptic transcript, the method allows the
  identification of cTSS. Five differents methods are available.
}
\details{
By definition, the observed f value for a gene is the perpendicular
 distance between the differential cumulative RNA-seq values (type1 - type2)
 and a diagonal linking the first and last data points. The simulated
 f max is the maximum f value for a gene after re-sampling the data.

Method A identify a cryptic zone by calculating positions for which
the observed f value is in the distribution of simulated f max.

Method B identify a cryptic zone by calculating positions for which the mean
simulated f value is in the distribution of simulated f max.

Method C identify a cryptic zone by calculating positions for which the
simulated f value is in the distribution of simulated f max.

Method D identify a cryptic zone by calculating the positions for each
simulated f max.

Method C gaussian determine the mean and standard deviation of all the
positions for which the simulated f value is in the distribution of simulated
f max.
}
\examples{
data("annotations")
samples <- c("wt_rep1", "wt_rep2", "mut_rep1", "mut_rep2")
bamfiles <- system.file("extdata", paste0(samples, ".bam"),
                        package = "yCrypticRNAs")
sf <- c(0.069872847, 0.081113079, 0.088520251, 0.069911116)
data(introns)
types = c("wt", "wt", "mut", "mut")


initiation_sites("YER109C", bamfiles, types, annotations,
                 introns, sf = sf, replicates = 5)
#initiation_sites("YER109C", bamfiles, types, annotations,
                # introns, sf = sf, percentage = 0.2, method = "methodA")
}


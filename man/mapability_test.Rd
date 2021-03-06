% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dottests.R
\name{mapability_test}
\alias{mapability_test}
\title{Perform the Low Mapability Test}
\usage{
mapability_test(ase, mapbed, cl = NULL, verbose = TRUE)
}
\arguments{
\item{ase}{An ASE data frame from \code{\link{read_ase}}}

\item{mapbed}{Path to a BED file with low mapability regions}

\item{cl}{
A cluster object created by \code{\link{makeCluster}},
or an integer to indicate number of child-processes
(integer values are ignored on Windows) for parallel evaluations
(see Details on performance).
}

\item{verbose}{Show progress updates}
}
\value{
\code{ase} with the following column added:
\itemize{
 \item \dQuote{\code{MAP_WARNING}}: Flag indicating if the gene  is
  associated with a variant that is in a region of low mapability
}
}
\description{
Perform the Low Mapability Test
}
\seealso{
Other ASE filtering tests:
\code{\link{bias_test}()},
\code{\link{homozygous_test}()},
\code{\link{monoallelic_test}()}
}
\concept{filter_test}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dottests.R
\name{bias_test}
\alias{bias_test}
\title{Perform the Mapping Bias Test}
\usage{
bias_test(ase, biasvcf, cl = NULL, verbose = TRUE)
}
\arguments{
\item{ase}{An ASE data frame from \code{\link{read_ase}}}

\item{biasvcf}{Path to VCF with bias sites}

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
 \item \dQuote{\code{BIAS_WARNING}}: Flag indicating if the gene is
  associated with a variant prone to mapping bias
}
}
\description{
Perform the Mapping Bias Test
}
\seealso{
Other ASE filtering tests:
\code{\link{homozygous_test}()},
\code{\link{mapability_test}()},
\code{\link{monoallelic_test}()}
}
\concept{filter_test}

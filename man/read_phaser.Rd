% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{read_phaser}
\alias{read_phaser}
\title{Read in phASER Gene AE Input}
\usage{
read_phaser(infile, verbose = TRUE)
}
\arguments{
\item{infile}{Path to input ASE data}

\item{verbose}{Show progress updates}
}
\value{
A \code{\link[data.table:data.table-class]{data.table}} with the
following columns:
\itemize{
 \item ...
}
followed by remaining columns from the input ASE table
}
\description{
Read in phASER Gene AE Input
}
\seealso{
Other functions for reading input ASE data:
\code{\link{read_ase}()},
\code{\link{read_gatk}()}
}
\concept{ase_input}

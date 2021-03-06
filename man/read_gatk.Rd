% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{read_gatk}
\alias{read_gatk}
\title{Read in GATK \code{ASEReadCounter} Input}
\usage{
read_gatk(infile, type = c("annotated"), verbose = TRUE)
}
\arguments{
\item{infile}{Path to input ASE data}

\item{type}{Type of ASE data being read in; choose from:
\itemize{
 \item \dQuote{\code{annotated}}: ASE annotated by \code{annotate_ase.py}
}}

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
Read in GATK \code{ASEReadCounter} Input
}
\seealso{
Other functions for reading input ASE data:
\code{\link{read_ase}()},
\code{\link{read_phaser}()}
}
\concept{ase_input}

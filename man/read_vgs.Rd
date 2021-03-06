% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.R
\name{read_vgs}
\alias{read_vgs}
\title{Read in Vg Estimates}
\usage{
read_vgs(vgfile, verbose = TRUE)
}
\arguments{
\item{vgfile}{Path to the Vg estimates}

\item{verbose}{Show progress updates}
}
\value{
A \code{\link[data.table:data.table-class]{data.table}} with two
columns:
\itemize{
 \item \dQuote{\code{GENE_ID}}: Gene IDs
 \item \dQuote{\code{Vg}}: Vg estimates
}
}
\description{
Read in Vg Estimates
}
\concept{data_input}

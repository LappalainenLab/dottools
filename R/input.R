#' @include zzz.R
#'
NULL

#' Read in ASE Input
#'
#' @param infile Path to input ASE data
#' @param type Type of ASE data being read in; choose from:
#' \itemize{
#'  \item \dQuote{\code{phaser}}: ASE from phASER Gene AE
#'  \item \dQuote{\code{annotated}}: ASE annotated by \code{annotate_ase.py}
#  \item \dQuote{\code{gatk}}: Raw ASE from \code{ASEReadCounter}
#' }
#' @param verbose Show progress updates
#'
#' @return A \code{\link[data.table:data.table-class]{data.table}} with the
#' following columns:
#' \itemize{
#'  \item ...
#' }
#' followed by remaining columns from the input ASE table
#'
#' @export
#'
#' @family ase_input
#'
read_ase <- function(
  infile,
  type = c('phaser', 'annotated'),
  verbose = TRUE
) {
  type <- type[1L]
  type <- match.arg(arg = type)
  if (!file.exists(infile)) {
    stop("Cannot find input ASE data ", infile)
  } else if (isTRUE(x = verbose)) {
    message(.msg("Reading ASE data", infile, "in", type, "format"))
  }
  ase <- fread(
    file = infile,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE,
    verbose = FALSE,
    showProgress = verbose
  )
  reqd.cols <- .reqd_cols(type = type)
  if (!all(reqd.cols %in% colnames(x = ase))) {
    stop("Input ASE missing some required columns")
  }
  ase <- ase[, reqd.cols, with = FALSE]
  colnames(x = ase) <- names(x = reqd.cols)
  return(switch(
    EXPR = type,
    'phaser' = {
      ase$gw_phased <- as.logical(x = ase$gw_phased)
      if (isTRUE(x = verbose)) {
        message(.msg("Adding 'OTHER_COUNT' from alternate variants"))
      }
      ase$OTHER_COUNT <- vapply(
        X = ase$n_variants,
        FUN = function(n) {
          return(max(n - 1L, 0L))
        },
        FUN.VALUE = integer(length = 1L)
      )
      ase
    },
    'gatk' = {
      ase
    },
    'annotated' = {
      ase
    }
  ))
}

#' Read in a Gene Blacklist
#'
#' @inheritParams read_ase
#' @param blfile Path to the gene blacklist
#'
#' @return A character vector of genes to exclude from analysis
#'
#' @export
#'
#' @concept data_input
#'
read_blacklist <- function(blfile, verbose = TRUE) {
  if (!file.exists(blfile)) {
    stop(.msg("Cannot find gene blacklist", blfile))
  } else if (isTRUE(x = verbose)) {
    message(.msg("Reading gene blacklist", blfile))
  }
  bl <- readLines(con = blfile)
  bl <- bl[!grepl(pattern = '^#', x = bl)]
  return(bl)
}

#' Read in GATK \code{ASEReadCounter} Input
#'
#' @inheritParams read_ase
#' @param type Type of ASE data being read in; choose from:
#' \itemize{
#'  \item \dQuote{\code{annotated}}: ASE annotated by \code{annotate_ase.py}
#  \item \dQuote{\code{gatk}}: Raw ASE from \code{ASEReadCounter}
#' }
#'
#' @inherit read_ase return
#'
#' @export
#'
#' @family ase_input
#'
read_gatk <- function(infile, type = c('annotated'), verbose = TRUE) {
  type <- type[1L]
  type <- match.arg(arg = type)
  return(read_ase(infile = infile, type = type, verbose = verbose))
}

#' Read in phASER Gene AE Input
#'
#' @inheritParams read_ase
#'
#' @inherit read_ase return
#'
#' @export
#'
#' @family ase_input
#'
read_phaser <- function(infile, verbose = TRUE) {
  return(read_ase(infile = infile, type = 'phaser', verbose = verbose))
}

#' Read in a Gene Symbol Map
#'
#' @inheritParams read_ase
#' @param mapfile Path to the gene ID/symbol map file; should be a
#' tab-separated file with ...
#'
#' @return A \code{\link[data.table:data.table-class]{data.table}} with two
#' columns:
#' \itemize{
#'  \item \dQuote{\code{GENE_ID}} ...
#'  \item \dQuote{\code{GENE_SYMBOL}} ...
#' }
#'
#' @export
#'
#' @concept data_input
#'
read_symbol_map <- function(mapfile, verbose = TRUE) {
  if (!file.exists(mapfile)) {
    stop("Cannot find gene ID/symbol map ", mapfile)
  } else if (isTRUE(x = verbose)) {
    message(.msg("Reading gene ID/symbol map", mapfile))
  }
  map <- fread(
    file = mapfile,
    sep = '\t',
    header = TRUE,
    stringsAsFactors = FALSE,
    verbose = FALSE,
    showProgress = verbose
  )
  colnames(x = map) <- c('GENE_ID', 'GENE_SYMBOL')
  return(map)
}

#' Read in Vg Estimates
#'
#' @inheritParams read_ase
#' @param vgfile Path to the Vg estimates
#'
#' @return A \code{\link[data.table:data.table-class]{data.table}} with two
#' columns:
#' \itemize{
#'  \item \dQuote{\code{GENE_ID}}: Gene IDs
#'  \item \dQuote{\code{Vg}}: Vg estimates
#' }
#'
#' @export
#'
#' @concept data_input
#'
read_vgs <- function(vgfile, verbose = TRUE) {
  if (!file.exists(vgfile)) {
    stop("Cannot find Vg estimates ", vgfile)
  } else if (isTRUE(x = verbose)) {
    message(.msg("Reading in Vg estimates", vgfile))
  }
  vgs <- fread(
    file = vgfile,
    sep = '\t',
    header = TRUE,
    select = c(1L, 2L),
    stringsAsFactors = FALSE,
    verbose = FALSE,
    showProgress = verbose
  )
  colnames(x = vgs) <- c('GENE_ID', 'Vg')
  return(vgs)
}

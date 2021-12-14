#' @importFrom data.table fread
#' @importFrom pbapply pblapply pboptions pbsapply
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Package options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @section Package Options:
#' \describe{
#'  \item{\code{dottools.filter.coverage}}{Minimum coverage for considering a
#'  gene in ANEVA-DOT}
#'  \item{\code{dottools.test.mono_cutoff}}{Other allele ratio for monoallelic
#'  test; see \code{\link{monoallelic_test}} for more details}
#'  \item{\code{dottools.test.mono_threshold}}{\emph{p}-value threshold for
#'  monoallelic test; see \code{\link{monoallelic_test}} for more details}
#'  \item{\code{dottools.test.homozygous_method}}{Method for running homozygous
#'  test; defaults to \dQuote{auto}}
#' }
#' @docType package
#'
"_PACKAGE"

dottools.options <- list(
  dottools.filter.coverage = 8L,
  dottools.test.mono_cutoff = 0.05,
  dottools.test.mono_threshold = 0.01,
  dottools.test.homozygous_method = 'auto'
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Reexports
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom ANEVADOT ANEVADOT_test
#' @export
#'
ANEVADOT::ANEVADOT_test

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create a Message
#'
#' @inheritParams base::paste
#'
#' @return A character vector with the message, wrapped to the terminal width
#'
#' @keywords internal
#'
#' @noRd
#'
.msg <- function(..., sep = ' ', collapse = NULL, recycle0 = FALSE) {
  return(paste(
    strwrap(x = paste(..., sep = sep, collapse = collapse, recycle0 = recycle0)),
    collapse = '\n'
  ))
}

#' Get Required Column Information
#'
#' @inheritParams read_ase
#'
#' @return A named vector with column names; values are the column names that
#' should be present in the input data, names are the names and order that the
#' ASE data frame can be renamed and reordered to
#'
#' @keywords internal
#'
#' @noRd
#'
.reqd_cols <- function(type = c('phaser', 'annotated', 'gatk')) {
  type <- type[1L]
  type <- match.arg(arg = type)
  reqd <- switch(
    EXPR = type,
    'phaser' = c(
      CHROM = 'contig',
      GENE_ID = 'name',
      REF_COUNT = 'aCount',
      ALT_COUNT = 'bCount',
      TOTAL_COUNT = 'totalCount',
      'start',
      'stop',
      'log2_aFC',
      'n_variants',
      VARIANTS = 'variants',
      'gw_phased',
      'bam'
    ),
    'annotated' = c(
      CHROM = 'CHR',
      GENE_ID = 'GENE_ID',
      REF_COUNT = 'REF_COUNT',
      ALT_COUNT = 'ALT_COUNT',
      TOTAL_COUNT = 'TOTAL_COUNT',
      'POS',
      'rsID',
      'REF_ALLELE',
      'ALT_ALLELE',
      'LOW_MAPQ_COUNT',
      'LOW_BASEQ_COUNT',
      'RAW_TOTAL_COUNT',
      OTHER_COUNT = 'OTHER_ALLELE_COUNT',
      'IMPROPER_PAIR_COUNT',
      VARIANTS = 'VARIANT_ID',
      'SAMPLE_ID',
      'REF_RATIO',
      'GENOTYPE',
      'LOW_MAPABILITY',
      'MAPPING_BIAS_SIM',
      'GENOTYPE_WARNING',
      'NULL_RATIO',
      'BINOM_P',
      'BINOM_P_ADJUSTED'
    ),
    'gatk' = c(
      CHROM = 'contig',
      REF_COUNT = 'refCount',
      ALT_COUNT = 'altCount',
      TOTAL_COUNT = 'totalCount',
      'position',
      'variantID',
      'refAllele',
      'altAllele',
      'lowMAPQDepth',
      'lowBaseQDepth',
      'rawDepth',
      OTHER_COUNT = 'otherBases',
      'improperPairs'
    )
  )
  unnamed <- which(x = !nchar(x = names(x = reqd)))
  names(x = reqd)[unnamed] <- reqd[unnamed]
  return(reqd)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Exported functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' ANEVA-DOT Script
#'
#' Get the path to the ANEVA-DOT script
#'
#' @return The path to the ANEVA-DOT script
#'
#' @keywords internal
#'
#' @export
#'
dot_script <- function() {
  return(system.file(
    'exec',
    'aneva-dot.R',
    package = 'dottools',
    mustWork = TRUE
  ))
}

#' Prepare a Vgs for ANEVA-DOT
#'
#' Generate a vector of Vgs for input into ANEVA-DOT
#'
#' @param ase A data frame with ASE from \code{\link{read_ase}}
#' @param vgs A data frame from \code{\link{read_vgs}}
#' @param verbose Show progress updates
#'
#' @return A named vector with the Vgs to use; names are the gene IDs, values
#' are the Vgs
#'
#' @export
#'
prep_vgs <- function(ase, vgs, verbose = TRUE) {
  if (isTRUE(x = verbose)) {
    message(.msg("Preparing Vgs for ANEVA-DOT"))
  }
  if (!'GENE_ID' %in% colnames(x = ase)) {
    stop("Malformed ASE data")
  } else if (!all(c('GENE_ID', 'Vg') %in% colnames(x = vgs))) {
    stop("Malformed Vg estimates")
  }
  vgs.all <- vgs$Vg
  names(x = vgs.all) <- vgs$GENE_ID
  genes.use <- intersect(x = names(x = vgs.all), y = ase$GENE_ID)
  if (!length(x = genes.use)) {
    stop("No genes matching between Vgs and ASE")
  }
  vgs.use <- vgs.all[genes.use]
  return(sqrt(x = vgs.use))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Load hooks
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  toset <- setdiff(x = names(x = dottools.options), y = names(x = options()))
  if (length(x = toset)) {
    options(dottools.options[toset])
  }
}

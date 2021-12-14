#' @include zzz.R
#'
NULL

#' Perform the Mapping Bias Test
#'
#' @inheritParams read_ase
#' @inheritParams pbapply::pbapply
#' @param ase An ASE data frame from \code{\link{read_ase}}
#' @param biasvcf Path to VCF with bias sites
#'
#' @return \code{ase} with the following column added:
#' \itemize{
#'  \item \dQuote{\code{BIAS_WARNING}}: Flag indicating if the gene is
#'   associated with a variant prone to mapping bias
#' }
#'
#' @export
#'
#' @family filter_test
#'
bias_test <- function(ase, biasvcf, cl = NULL, verbose = TRUE) {
  if (!isTRUE(x = verbose)) {
    x <- pboptions(type = NULL)
    on.exit(expr = pboptions(x))
  }
  if (!file.exists(biasvcf)) {
    stop("Cannot find mapping VCF file ", biasvcf)
  } else if (!'VARIANTS' %in% colnames(x = ase)) {
    stop("No variant information present in ase")
  }
  vcf <- fread(
    file = biasvcf,
    sep = '\t',
    header = TRUE,
    skip = '##',
    verbose = FALSE,
    showProgress = verbose
  )
  bvars <- unique(x = paste(vcf$`#CHROM`, vcf$POS, sep = '_'))
  ase$BIAS_WARNING <- pbsapply(
    X = ase$VARIANTS,
    FUN = function(vars) {
      vars <- unlist(x = strsplit(x = vars, split = ','), use.names = FALSE)
      pos <- vector(mode = 'character', length = length(x = vars))
      for (i in seq_along(along.with = vars)) {
        vi <- unlist(x = strsplit(x = vars[i], split = '_'), use.names = FALSE)
        pos[i] <- paste(vi[c(1L, 2L)], collapse = '_')
      }
      return(any(pos %in% bvars))
    },
    USE.NAMES = FALSE,
    simplify = TRUE,
    cl = cl
  )
  rm(vcf)
  return(ase)
}

#' Perform the Homozygous Test
#'
#' @inheritParams bias_test
#' @param method Method for running the homozygous test; choose from
#' \itemize{
#'  \item \dQuote{\code{genotype}}: Use genotype information for determining
#'   heterozygous/homozygous status
#'  \item \dQuote{\code{counts}}: Use ref/alt counts to determine
#'   heterozygous/homozygous status
#'  \item \dQuote{\code{auto}}: If a \dQuote{\code{GENOTYPE}} column present,
#'   use the \dQuote{\code{genotype}} method; otherwise, use the
#'   \dQuote{\code{counts}} method
#' }
#' @param missing For \code{method == "genotype"}, behavior for missing
#' genotypes; choose from:
#' \itemize{
#'  \item \dQuote{\code{homozygous}}: Treat missing genotypes as homozygous
#'  \item \dQuote{\code{na}}: Mark missing genotypes with \code{NA}
#' }
#' @param other For \code{method == "counts"}, include other counts with alt
#' counts in heterozygous/homozygous test
#'
#' @return \code{ase} with the following column added:
#' \itemize{
#'  \item \dQuote{\code{HOMOZYGOUS_WARNING}}: Flag to indicate if the variant
#'   associated with the gene is homozygous
#' }
#'
#' @export
#'
#' @family filter_test
#'
homozygous_test <- function(
  ase,
  method = getOption(x = 'dottools.test.homozygous_method'),
  missing = c('homozygous', 'na'),
  other = TRUE,
  cl = NULL,
  verbose = TRUE
) {
  if (!isTRUE(x = verbose)) {
    x <- pboptions(type = NULL)
    on.exit(expr = pboptions(x))
  }
  method <- method[1L]
  method <- match.arg(arg = method, choices = c('auto', 'genotype', 'counts'))
  if (method == 'auto') {
    method <- ifelse(
      test = 'GENOTYPE' %in% colnames(x = ase),
      yes = 'genotype',
      no = 'counts'
    )
  }
  ase$HOMOZYGOUS_WARNING <- switch(
    EXPR = method,
    'genotype' = {
      if (!'GENOTYPE' %in% colnames(x = ase)) {
        stop("No genotype information present in ase")
      }
      missing <- missing[1L]
      missing <- match.arg(arg = missing)
      missing.value <- switch(EXPR = missing, 'homozygous' = TRUE, 'na' = NA)
      pbsapply(
        X = ase$GENOTYPE,
        FUN = function(x) {
          x <- gsub(pattern = '/|\\|', replacement = ',', x = x)
          x <- unlist(x = strsplit(x = x, split = ','))
          if (any(grepl(pattern = '\\.', x = x))) {
            return(missing.value)
          }
          return(length(x = unique(x = x)) == 1L)
        },
        USE.NAMES = FALSE,
        cl = cl
      )
    },
    'counts' = {
      if (!all(c('REF_COUNT', 'ALT_COUNT') %in% colnames(x = ase))) {
        stop("No counts information in ase")
      }
      if (isTRUE(x = other) && 'OTHER_COUNT' %in% colnames(x = ase)) {
        ase$ALT_OTHER <- ase$ALT_COUNT + ase$OTHER_COUNT
        alt <- 'ALT_OTHER'
      } else {
        alt <- 'ALT_COUNT'
      }
      !(ase$REF_COUNT != 0 & ase[[alt]] != 0)
    }
  )
  return(ase)
}

#' Perform the Low Mapability Test
#'
#' @inheritParams bias_test
#' @param mapbed Path to a BED file with low mapability regions
#'
#' @return \code{ase} with the following column added:
#' \itemize{
#'  \item \dQuote{\code{MAP_WARNING}}: Flag indicating if the gene  is
#'   associated with a variant that is in a region of low mapability
#' }
#'
#' @export
#'
#' @family filter_test
#'
mapability_test <- function(ase, mapbed, cl = NULL, verbose = TRUE) {
  if (!isTRUE(x = verbose)) {
    x <- pboptions(type = NULL)
    on.exit(expr = pboptions(x))
  }
  if (!file.exists(mapbed)) {
    stop("Cannot find mapability BED file ", mapbed)
  } else if (!'VARIANTS' %in% colnames(x = ase)) {
    stop("No variant information present in ase")
  }
  bed <- fread(
    file = mapbed,
    sep = '\t',
    header = FALSE,
    select = seq.int(from = 1L, to = 3L),
    stringsAsFactors = FALSE,
    verbose = FALSE,
    showProgress = verbose
  )
  colnames(x = bed) <- c('CHROM', 'START', 'END')
  ase$MAP_WARNING <- pbsapply(
    X = ase$VARIANTS,
    FUN = function(vars) {
      vars <- unlist(x = strsplit(x = vars, split = ','), use.names = FALSE)
      vchroms <- unique(x = vapply(
        X = strsplit(x = vars, split = '_'),
        FUN = '[[',
        FUN.VALUE = character(length = 1L),
        1L
      ))
      map.chrom <- map.chrom[map.chrom$CHROM %in% vchroms, , drop = FALSE]
      low.site <- FALSE
      for (v in vars) {
        pos <- as.integer(x = unlist(x = strsplit(x = v, split = '_'))[2L])
        if (any(pos >= map.chrom$START & map.chrom$END >= pos)) {
          low.site <- TRUE
          break
        }
      }
      return(low.site)
    },
    cl = cl
  )
  rm(bed)
  return(ase)
}

#' Perform the Monoallelic Test
#'
#' @inheritParams bias_test
#' @param coverage ...
#' @param cutoff Maximum proportion of reads arising from non ref/alt read for
#' variant to be counted
#' @param threshold \emph{p}-value threshold
#'
#' @return \code{ase} with the following columns added:
#' \itemize{
#'  \item \dQuote{\code{mono_pvals}}: \emph{p}-values for the monoallelic test
#'  \item \dQuote{\code{mono_qvals}}: FDR-adjusted \emph{p}-values
#'   (\emph{q}-values) for the monoallelic test
#'  \item \dQuote{\code{GENOTYPE_WARNING}}: Flag to indicate if the gene fails
#'   the monoallelic test
#' }
#'
#' @importFrom stats pbinom p.adjust
#'
#' @export
#'
#' @family filter_test
#'
monoallelic_test <- function(
  ase,
  coverage = getOption(x = 'dottools.filter.coverage'),
  cutoff = getOption(x = 'dottools.test.mono_cutoff'),
  threshold = getOption(x = 'dottools.test.mono_threshold'),
  verbose = TRUE
) {
  if (!'OTHER_COUNT' %in% colnames(x = ase)) {
    stop("No count of other alleles in ASE data")
  }
  # TODO: move this to aneva_dot
  ase <- ase[ase$TOTAL_COUNT >= coverage, , drop = FALSE]
  other.ratio <- ase$OTHER_COUNT / (ase$TOTAL_COUNT + ase$OTHER_COUNT)
  all.count <- ase$TOTAL_COUNT + ase$OTHER_COUNT
  sum.index <- which(x = other.ratio < cutoff)
  sum.total <- sum(all.count[sum.index])
  sum.other <- sum(ase$OTHER_COUNT[sum.index])
  lamp <- (sum.other / sum.total) / 2
  if (isTRUE(x = verbose)) {
    message(.msg("LAMP =", lamp))
  }
  ase$mono_pvals <- pbinom(
    q = ase$TOTAL_COUNT - ase$REF_COUNT,
    size = ase$TOTAL_COUNT,
    prob = 1 - lamp
  ) + pbinom(
    q = ase$REF_COUNT,
    size = ase$TOTAL_COUNT,
    prob = 1 - lamp
  )
  ase$mono_qvals <- p.adjust(p = ase$mono_pvals, method = 'fdr')
  ase$GENOTYPE_WARNING <- ase$mono_qvals >= threshold
  return(ase)
}

#!/usr/bin/env Rscript

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Dependency Checking
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Check dependencies
deps <- c('dottools', 'argparser', 'tools', 'grDevices')
pkgcheck <- vapply(
  X = deps,
  FUN = requireNamespace,
  FUN.VALUE = logical(length = 1L),
  quietly = TRUE
)

if (!all(pkgcheck)) {
  stop(
    "Cannot find the following packages: ",
    paste(deps[!pkgcheck], collapse = ', ')
  )
}

library(dottools)
library(argparser)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Argument Parsing
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Set arguments
parser <- arg_parser(description = "Run ANVEA-DOT", hide.opts = TRUE)

# Basic I/O
parser <- add_argument( # ASE Input
  parser = parser,
  arg = '--input',
  help = "Input ASE file from phASER Gene AE or annotate_ase.py",
  default = NA_character_,
  short = '-i'
)
parser <- add_argument( # Vg Estimates
  parser = parser,
  arg = '--vg-file',
  help = "Two-column, tab-separated Vg file; requires gene IDs in the first column and Vg estimates in the second",
  default = NA_character_,
  short = '-v'
)
parser <- add_argument( # Output Directory
  parser = parser,
  arg = '--outdir',
  help = "Name of output directory",
  default = file.path(getwd(), 'dot_scores'),
  short = '-o'
)

# Extra inputs/information
parser <- add_argument( # Gene ID/Symbol Map
  parser = parser,
  arg = '--symbol-map',
  help = "Two-column, tab-separated gene ID/symbol map; requires the gene IDs matching the Vg file in the first column and the gene symbols matching the ASE data in the second",
  default = NA_character_,
  short = '-m'
)
parser <- add_argument( # Gene Blacklist
  parser = parser,
  arg = '--blacklist',
  help = "List of blacklisted genes; applied after the ID/symbol conversion",
  default = NA_character_,
  short = '-b'
)
parser <- add_argument( # ASE Input Format
  parser = parser,
  arg = '--ase-type',
  help = "Type of ASE being input; choose from 'phaser' or 'annotated'",
  default = 'phaser'
)

# Filtering
parser <- add_argument( # Coverage
  parser = parser,
  arg = '--coverage',
  help = "Minimum coverage for including a gene in ANEVA-DOT",
  default = 20L,
  short = '-c'
)
parser <- add_argument( # Count Threshold
  parser = parser,
  arg = '--count-threshold',
  help = "Count threshold",
  default = 8L,
  short = '-t'
)
parser <- add_argument( # Filter Homozygous
  parser = parser,
  arg = '--filter-homozygous',
  help = "Filter homozygous sites", # TODO: include details about req'd cols
  flag = TRUE
)
parser <- add_argument( # Run the homozygous test
  parser = parser,
  arg = '--filter-homozygous-test',
  help = "Force run the homozygous test",
  flag = TRUE
)
parser <- add_argument( # Mapping Bias
  parser = parser,
  arg = '--filter-bias',
  help = "Filter sites based on mapping bias",
  flag = TRUE
)
parser <- add_argument( # Mapping Bias VCF File
  parser = parser,
  arg = '--filter-bias-vcf',
  help = "VCF with sites prone to mapping bias; implies '--filter-bias'",
  default = NA_character_
)
parser <- add_argument( # Low Mapability
  parser = parser,
  arg = '--filter-mapability',
  help = "Filter regions of low mapability",
  flag = TRUE
)
parser <- add_argument( # Low Mapability BED File
  parser = parser,
  arg = '--filter-mapability-bed',
  help = "BED file with regions of low mapability; implies '--filter-mapability'",
  default = NA_character_
)
parser <- add_argument( # Monoallelic Test
  parser = parser,
  arg = '--filter-monoallelic',
  help = "Filter sites that are likely monoallelic",
  flag = TRUE
)
parser <- add_argument( # Force Monoallelic Test
  parser = parser,
  arg = '--filter-monoallelic-test',
  help = "Recalculate monoallelic test if not already done; implies '--filter-monoallelic'",
  flag = TRUE
)
parser <- add_argument( # Monoallelic Test Cutoff
  parser = parser,
  arg = '--mono-cutoff',
  help = "mono cutoff",
  default = 0.05
)
parser <- add_argument( # Monoallelic Test Threshold
  parser = parser,
  arg = '--mono-threshold',
  help = "mono threshold",
  default = 0.01
)

# Flags
parser <- add_argument( # Generate plots
  parser = parser,
  arg = '--plot',
  help = "Generate ANEVA-DOT plots",
  flag = TRUE
)
parser <- add_argument( # Show progress updates
  parser = parser,
  arg = '--quiet',
  help = "Suppress progress updates",
  flag = TRUE
)

if (!length(x = commandArgs(trailingOnly = TRUE))) {
  print(x = parser)
  stop("Not enough arguments provided")
}
args <- parse_args(parser = parser)

# Check required arguments
if (is.na(x = args$input)) {
  stop("'-i|--input' is required", call. = FALSE)
} else if (is.na(x = args$vg_file)) {
  stop("'-v|--vg-file' is required", call. = FALSE)
}
verbose <- !args$quiet

# Adjust filtering arguments
if (!is.na(x = args$filter_bias_vcf)) {
  args$filter_bias <- TRUE
}
if (!is.na(x = args$filter_mapability_bed)) {
  args$filter_mapability <- TRUE
}
if (isTRUE(x = args$filter_homozyous_test)) {
  args$filter_homozygous <- TRUE
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# I/O Management
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Read in input data
ase <- read_ase(infile = args$input, type = args$ase_type, verbose = verbose)
vgs <- read_vgs(vgfile = args$vg_file, verbose = verbose)
if (!is.na(x = args$blacklist)) {
  blacklist <- read_blacklist(blfile = args$blacklist, verbose = verbose)
}

# Prepare output directory and names
dir.create(path = args$outdir, showWarnings = FALSE, recursive = TRUE)
bname <- paste0(
  tools::file_path_sans_ext(x = basename(path = args$input)),
  '_dot_scores'
)
bname <- file.path(args$outdir, bname)
ofile <- paste0(bname, '.tsv')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ANEVA-DOT Preparation
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Translate symbols
if (!is.na(x = args$symbol_map)) {
  if (isTRUE(x = verbose)) {
    message("Adding gene IDs from gene ID/symbol map")
  }
  symbols <- read_symbol_map(mapfile = args$symbol_map, verbose = verbose)
  colnames(x = ase) <- gsub(
    pattern = 'GENE_ID',
    replacement = 'GENE_SYMBOL',
    x = colnames(x = ase)
  )
  ase <- merge(x = ase, y = symbols, by = 'GENE_SYMBOL', all.x = TRUE)
}

# TODO: filter blacklists

# Perform filtering
## Filter homozygous sites
if (isTRUE(x = args$filter_homozygous)) {
  if (isTRUE(x = verbose)) {
    message("Filtering sites that are homozygous")
  }
  do.homozygous <- isTRUE(x = args$filter_homozygous_test) ||
    !'HOMOZYGOUS_WARNING' %in% colnames(x = ase)
  if (isTRUE(x = do.homozygous)) {
    ase <- tryCatch(
      expr = {
        ase <- homozygous_test(ase = ase, verbose = verbose)
        ase[!ase$HOMOZYGOUS_WARNING, , drop = FALSE]
      },
      error = function(e) {
        warning(e$message, call. = FALSE, immediate. = TRUE)
        return(ase)
      }
    )
  }
  if ('HOMOZYGOUS_WARNING' %in% colnames(x = ase)) {
    ase <- ase[!ase$HOMOZYGOUS_WARNING, , drop = FALSE]
  } else {
    warning("No homozygous status provided in ASE data", immediate. = TRUE)
  }
}
## Filter by count
if (isTRUE(x = verbose)) {
  message(
    "Filtering sites by count; using a threshold of ",
    args$count_threshold
  )
} else {
  old.pbopts <- pboptions(type = NULL)
}
ase <- pblapply(
  X = unique(x = ase$GENE_ID),
  FUN = function(x) {
    df <- ase[ase$GENE_ID == x, , drop = FALSE]
    df <- df[df$TOTAL_COUNT > args$count_threshold, , drop = FALSE]
    # if (nrow(x = df) > 1L) {
    #   df <- df[df$TOTAL_COUNT < max(df$TOTAL_COUNT), , drop = FALSE]
    # }
    return(df)
  }
)
if (!isTRUE(x = verbose)) {
  pboptions(old.pbopts)
}
ase <- do.call(what = 'rbind', args = ase)
## Filter sites prone to mapping bias
if (isTRUE(x = args$filter_bias)) {
  if (isTRUE(x = verbose)) {
    message("Filtering sites prone to mapping bias")
  }
  if (!is.na(x = args$filter_bias_vcf)) {
    ase <- tryCatch(
      expr = bias_test(
        ase = ase,
        biasvcf = args$filter_bias_vcf,
        verbose = verbose
      ),
      error = function(e) {
        warning(e$message, call. = FALSE, immediate. = TRUE)
        return(ase)
      }
    )
  }
  if ('BIAS_WARNING' %in% colnames(x = ase)) {
    ase <- ase[!ase$BIAS_WARNING, , drop = FALSE]
  } else {
    warning("No mapping bias status provided in ASE data", immediate. = TRUE)
  }
}
### Filter regions of low mapability
if (isTRUE(x = args$filter_mapability)) {
  if (isTRUE(x = verbose)) {
    message("Filtering regions of low mapability")
  }
  if (!is.na(x = args$filter_mapability_bed)) {
    ase <- tryCatch(
      expr = mapability_test(
        ase = ase,
        mapbed = args$filter_mapability_bed,
        verbose = verbose
      ),
      error = function(e) {
        warning(e$message, call. = FALSE, immediate. = TRUE)
        return(ase)
      }
    )
  }
  if ('MAP_WARNING' %in% colnames(x = ase)) {
    ase <- ase[!ase$MAP_WARNING, , drop = FALSE]
  } else {
    warning("No mapability status provided in ASE data", immediate. = TRUE)
  }
}
## Filter monoallelic sites
if (isTRUE(x = args$filter_monoallelic)) {
  if (isTRUE(x = verbose)) {
    message("Filtering monallelic sites")
  }
  if (!'GENOTYPE_WARNING' %in% colnames(x = ase)) {
    args$filter_monoallelic_test <- TRUE
  }
  if (isTRUE(x = args$filter_monoallelic_test)) {
    ase <- monoallelic_test(ase = ase, verbose = verbose)
  }
  ase <- ase[!ase$GENOTYPE_WARNING, , drop = FALSE]
}

# Match up Vgs with genes tested
vgs.use <- prep_vgs(ase = ase, vgs = vgs, verbose = verbose)
ase <- ase[match(x = names(x = vgs.use), table = ase$GENE_ID), , drop = FALSE]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Run ANEVA-DOT
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Prepare plotting
if (isTRUE(x = args$plot)) {
  oplot <- paste0(bname, '.pdf')
  if (isTRUE(x = verbose)) {
    message("Saving plot to ", oplot)
  }
  grDevices::pdf(file = oplot)
}

# Run ANEVA-DOT
dot.scores <- ANEVADOT_test(
  ASEdat = ase,
  output_columns = colnames(x = ase),
  eh1 = 'REF_COUNT',
  eh2 = 'ALT_COUNT',
  Eg_std = vgs.use,
  r0 = ase$NULL_RATIO,
  coverage = args$coverage,
  plot = args$plot
)

# Close the plotting device
if (isTRUE(x = args$plot)) {
  invisible(x = grDevices::dev.off())
}

# Write out the scores
if (isTRUE(x = verbose)) {
  message("Writing DOT scores to ", ofile)
}
data.table::fwrite(
  x = dot.scores,
  file = ofile,
  quote = FALSE,
  sep = '\t',
  row.names = FALSE,
  col.names = TRUE
)

#!/usr/bin/env Rscript
# builds UPGMA (phangorn) or NJ (ape) trees
# writes the resulting newick to --out-tree, summary stats to --out-stats

suppressPackageStartupMessages({
  library(optparse)
  library(phangorn)
  library(ape)
  library(jsonlite)
  library(rlang)
})

#### CLI parsing ####
make_parser <- function() {
  option_list <- list(
    make_option("--fasta", type = "character", help = "aligned FASTA"),
    make_option("--datatype", type = "character", help = "dna or aa"),
    make_option("--model", type = "character", help = "substitution model"),
    make_option("--method", type = "character", help = "UPGMA or NJ"),
    make_option("--out-tree", type = "character", help = ".nwk path"),
    make_option("--out-stats", type = "character", help = "stats.json path")
  )
  OptionParser(option_list = option_list)
}

#### utilities ####
ensure_dir <- function(path) {
  d <- dirname(path)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
} # ensures dir exists or creates it

safe_write_json <- function(obj, path) {
  ensure_dir(path)
  writeLines(toJSON(obj, auto_unbox = TRUE, pretty = TRUE),
             con = path, sep = "\n")
} # writes json using jsonlite

elapsed_timer <- function() {
  t0 <- proc.time()[["elapsed"]]
  function() round(proc.time()[["elapsed"]] - t0, 3)
} # measure runtime

normalize_model <- function(datatype, model) {
  model <- toupper(trimws(model %||% ""))
  if (datatype == "dna") {
    # DNA models
    dna_ok <- c("JC69", "F81", "K80", "HKY", "TN93", "GTR")
    if (model %in% dna_ok) return(model)
    "JC69" # use JC69 if invalid model
  } else {
    # AA models
    aa_ok <- c("JTT", "LG", "WAG", "DAYHOFF")
    if (model %in% aa_ok) return(model)
    "JTT" # use JTT if invalid model
  }
} # selects supported base substitution model

#### core method functions ####

read_alignment <- function(fasta, datatype) {
  type <- if (datatype == "dna") "DNA" else "AA"
  phy <- read.phyDat(file = fasta, format = "fasta", type = type)
  if (is.null(phy) || length(phy) < 3) {
    stop("alignment could not be read into phyDat or has <3 sequences.")
  }
  phy
} # loads aligned FASTA into phyDat object

compute_distance <- function(phy, datatype, model) {
  mdl <- normalize_model(datatype, model)
  phangorn::dist.ml(phy, model = mdl)
} # computes ML distance matrix with model

build_tree_from_dist <- function(d, method) {
  method <- toupper(trimws(method))
  if (method == "UPGMA") {
    phangorn::upgma(d)
  } else if (method == "NJ") {
    ape::nj(d)
  } else {
    stop(sprintf("unsupported method for distance tree: %s", method))
  }
} # runs phangorn::upgma() or ape::nj() using distance matrix

write_tree <- function(tree, out_path) {
  ensure_dir(out_path)
  write.tree(tree, file = out_path)
} # creates .nwk file

#### final output functions ####
gather_stats <- function(phy, d, tree, args, elapsed_sec) {
  nseq <- length(phy)
  # for phyDat, alignment length can be derived from attr
  w <- attr(phy, "weight")
  nsites <- if (!is.null(w)) sum(w) else (attr(phy, "nr") %||% NA_integer_)

  dist_vals <- as.numeric(d)
  has_vals <- length(dist_vals) > 0
  list(
    method = toupper(args$method),
    model = toupper(args$model %||% ""),
    datatype = tolower(args$datatype),
    n_sequences = nseq,
    n_sites = nsites,
    dist_mean = if (has_vals) mean(dist_vals, na.rm = TRUE) else NA_real_,
    dist_min = if (has_vals) min(dist_vals,  na.rm = TRUE) else NA_real_,
    dist_max = if (has_vals) max(dist_vals,  na.rm = TRUE) else NA_real_,
    is_ultrametric = is.ultrametric(tree),
    is_rooted = is.rooted(tree),
    tip_count = length(tree$tip.label),
    elapsed_sec = elapsed_sec,
    versions = list(
      R = as.character(getRversion()),
      ape = as.character(packageVersion("ape")),
      phangorn = as.character(packageVersion("phangorn"))
    )
  )
}
# summarises method, model, datatype, fasta stats, distances distribution,
# ultrametric or rooted, tip count, elapsed time, package versions.
# results written by safe_write_json() and will go to to stats.json

main <- function() {
  parser <- make_parser()
  opts <- parse_args(parser)

  # minimal sanity, heavy validation already done in validate_config.py
  required <- c("fasta", "datatype", "model", "method", "out-tree", "out-stats")
  ok <- vapply(required, function(k) {
    v <- opts[[k]]
    !is.null(v) && length(v) == 1 && nzchar(as.character(v))
  }, logical(1))
  if (!all(ok)) {
    missing <- paste(required[!ok], collapse = ", ")
    stop(sprintf("Missing required arguments: %s", missing))
  }

  datatype <- tolower(opts$datatype)
  if (!(datatype %in% c("dna", "aa"))) {
    stop("datatype must be 'dna' or 'aa'")
  }

  tictoc <- elapsed_timer()

  # read alignment
  phy <- read_alignment(opts$fasta, datatype)
  # distance matrix
  d <- compute_distance(phy, datatype, opts$model)
  # make tree
  tree <- build_tree_from_dist(d, opts$method)
  # outputs
  write_tree(tree, opts[["out-tree"]])
  stats <- gather_stats(phy, d, tree, opts, tictoc())
  safe_write_json(stats, opts[["out-stats"]])
}

# wrap main() in trycatch so snakemake detects failure
tryCatch(
  { main() },
  error = function(e) {
    message(sprintf("[build_dist_tree.R] error: %s", conditionMessage(e)))
    quit(status = 1)
  }
)

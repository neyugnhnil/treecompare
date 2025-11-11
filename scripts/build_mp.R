#!/usr/bin/env Rscript
# builds MP tree using phangorn's ratchet search
# optionally starts from previously generated tree
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
    make_option("--starter", type = "character", help = "'random' or path"),
    make_option("--out-tree", type = "character", help = ".nwk path"),
    make_option("--out-stats", type = "character", help = "stats.json path"),
    make_option("--seed", type = "integer", help = "optional RNG seed")
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

##### starter tree handling ####
# loads starter .nwk, validates it, reorders tips
# falls back to random (NULL) with warnings if anything's off
read_starter_tree <- function(starter, expected_labels) {
  # starter is "random" or a path to .nwk (dispatcher resolved names to paths)
  if (is.null(starter) || tolower(trimws(starter)) == "random") {
    return(NULL)  # random
  }
  if (!file.exists(starter)) {
    warning(sprintf("starter not found: %s; defaulting to random", starter))
    return(NULL)
  }
  tr <- tryCatch(ape::read.tree(starter), error = function(e) NULL)
  if (is.null(tr)) {
    warning(sprintf("starter invalid (%s); defaulting to random", starter))
    return(NULL)
  }
  if (length(tr$tip.label) < 3) {
    warning("starter has <3 tips; defaulting to random")
    return(NULL)
  }
  if (!setequal(sort(tr$tip.label), sort(expected_labels))) {
    warning("starter labels do not match FASTA headers; defaulting to random")
    return(NULL)
  }
  # reorder to alignment tip order
  tr <- keep.tip(tr, expected_labels)
  tr
}

#### core method functions ####
read_alignment <- function(fasta, datatype) {
  # get phyDat object
  type <- if (datatype == "dna") "DNA" else "AA"
  phy <- read.phyDat(file = fasta, format = "fasta", type = type)
  if (is.null(phy) || length(phy) < 3) {
    stop("Alignment could not be read into phyDat or has <3 sequences.")
  }
  phy
} # loads aligned FASTA into phyDat object

build_mp_tree <- function(phy, starter_tree = NULL) {
  # if starter provided, use it; otherwise random starting point
  if (is.null(starter_tree)) {
    best <- phangorn::pratchet(phy, trace = 0)
    # pratchet generates random tree (random addition sequence) by default
  } else {
    best <- phangorn::pratchet(phy, start = starter_tree, trace = 0)
  }
  best # branch length improvements unnecessary for topology
} # heuristic parsimony search (Ratchet)

write_tree <- function(tree, out_path) {
  ensure_dir(out_path)
  write.tree(tree, file = out_path)
} # creates .nwk file

#### final output functions ####
gather_stats <- function(phy, tree, args, starter_kind, elapsed_sec) {
  nseq <- length(phy)

  # alignment length from phyDat weights
  w <- attr(phy, "weight")
  nsites <- if (!is.null(w)) sum(w) else (attr(phy, "nr") %||% NA_integer_)
  # parsimony score (number of steps)
  score <- as.numeric(parsimony(tree, phy))
  # Consistency Index (CI), Retention Index (RI)
  ci <- tryCatch(CI(tree, phy), error = function(e) NA_real_)
  ri <- tryCatch(RI(tree, phy), error = function(e) NA_real_)

  list(
    method = "MP",
    datatype = tolower(args$datatype),
    n_sequences = nseq,
    n_sites = nsites,
    parsimony_score = score,
    CI = ci,
    RI = ri,
    tip_count = length(tree$tip.label),
    starter = starter_kind,     # "random" or "provided"
    elapsed_sec = elapsed_sec,
    versions = list(
      R = as.character(getRversion()),
      ape = as.character(packageVersion("ape")),
      phangorn = as.character(packageVersion("phangorn"))
    )
  )
}
# summarises method, datatype, fasta stats, parsimony score, CI/RI indices,
# ultrametric or rooted, starter, elapsed time, package versions.
# results written by safe_write_json() and will go to to stats.json

main <- function() {
  parser <- make_parser()
  opts <- parse_args(parser)

  # minimal sanity, heavy validation already done in validate_config.py
  required <- c("fasta", "datatype", "starter", "out-tree", "out-stats")
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

  # seed will apply to both ratchet search and starter if random
  if (!is.null(opts$seed) && !is.na(opts$seed) && opts$seed > 0) {
    set.seed(as.integer(opts$seed))
  }

  # read alignment
  phy <- read_alignment(opts$fasta, datatype)
  # starter handling
  labels <- names(phy)
  starter_tree <- read_starter_tree(opts$starter, labels)
  starter_kind <- if (is.null(starter_tree)) "random" else "provided"
  # build MP tree
  tree <- build_mp_tree(phy, starter_tree)
  # outputs
  write_tree(tree, opts[["out-tree"]])
  stats <- gather_stats(phy, tree, opts, starter_kind, tictoc())
  safe_write_json(stats, opts[["out-stats"]])
}

# wrap main() in trycatch so snakemake detects failure
tryCatch(
  { main() },
  error = function(e) {
    message(sprintf("[build_mp.R] error: %s", conditionMessage(e)))
    quit(status = 1)
  }
)

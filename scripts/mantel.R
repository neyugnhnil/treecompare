#!/usr/bin/env Rscript
# given --tree and --meta, correlate topo. distance matrix against feature distance matrices (numeric euclidean or categorical mismatch)
# outputs tsv with metrics for each feature
# assumes inputs are validated upstream

suppressPackageStartupMessages({
  library(optparse)
  library(ape)
  library(vegan)
})

#### CLI parsing #######
make_parser <- function() {
  option_list <- list(
    make_option("--tree", type = "character", help = ".nwk"),
    make_option("--meta", type = "character", help = ".tsv metafile"),
    make_option("--method", type = "character", default = "spearman",
                help = "spearman or pearson"),
    make_option("--permutations", type = "integer", default = 1000,
                help = "number of permutations"),
    make_option("--out", type = "character", help = "output TSV")
  )
  OptionParser(option_list = option_list)
}
#### utilities #########
ensure_dir <- function(path) {
  d <- dirname(path)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
} # ensures dir exists or creates it

cat_dist <- function(v) {
  v <- as.character(v)
  as.dist(outer(v, v, FUN = function(a, b) as.integer(a != b)))
} # builds pairwise 0/1 mismatch distance matrix for a single vector

read_any_tree <- function(path) {
  header <- tryCatch(readLines(path, n = 1), error = function(e) character())
  reader <- if (length(header) > 0 && grepl("^#?\\s*NEXUS", toupper(header[1]))) {
    ape::read.nexus
  } else {
    ape::read.tree
  }
  tree <- reader(path)
  if (inherits(tree, "multiPhylo")) tree <- tree[[1]]
  tree
} # can read either nexus or newick trees

###### main logic #######
main <- function() {
  opt <- parse_args(make_parser())
  tr <- read_any_tree(opt$tree)

  tr$edge.length <- rep(1, nrow(tr$edge)) # sets branch lengths to 1
  dtree <- cophenetic.phylo(tr) # topological distance matrix

  # trust validate_config.py
  meta <- read.table(opt$meta, sep = "\t", header = TRUE, quote = "",
                     comment.char = "", check.names = FALSE,
                     stringsAsFactors = FALSE)
  ids <- meta[[1]]

  # keep only metafile rows present in the tree when subsetting per feature
  in_tree <- ids %in% rownames(dtree)
  meta <- meta[in_tree, , drop = FALSE]
  ids <- meta[[1]] # list of headers

  res <- vector("list", ncol(meta) - 1) # a list to hold dfs (1 per feature)
  k <- 0
  for (j in seq.int(2, ncol(meta))) { # iterate over cols (skip col 1)
    feat <- colnames(meta)[j] # feature name
    v <- meta[[j]] # vector of values

    # drop rows with missing vals for this feature
    ok <- !is.na(v)
    ids_ok <- ids[ok]
    v_ok <- v[ok]

    # are there still 3 or more taxa to work with? (NA in metrics if no)
    if (length(ids_ok) < 3) {
      k <- k + 1
      res[[k]] <- data.frame(feature = feat, r = NA_real_, p = NA_real_,
                             method = tolower(opt$method),
                             permutations = opt$permutations,
                             n = length(ids_ok))
      next
    }

    dtree_sub <- # dtree subsetted to kept taxa
      as.dist(dtree[ids_ok, ids_ok, drop = FALSE])

    dfeat <- # feature value distance matrix
      if (is.numeric(v_ok)) dist(v_ok, method = "euclidean") else cat_dist(v_ok)

    # run mantel test
    mt <- vegan::mantel(dtree_sub, dfeat,
                        method = tolower(opt$method),
                        permutations = opt$permutations)

    k <- k + 1

    # store results for this feature
    res[[k]] <- data.frame(feature = feat,
                           r = unname(mt$statistic),
                           p = unname(mt$signif),
                           method = tolower(opt$method),
                           permutations = opt$permutations,
                           n = length(ids_ok))
  }

  # output
  out <- if (k > 0) do.call(rbind, res[seq_len(k)]) else
    data.frame(feature = character(), r = numeric(), p = numeric(),
               method = character(), permutations = integer(), n = integer())

  ensure_dir(opt$out)
  write.table(out, file = opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
}

tryCatch(main(), error = function(e) {
  message(sprintf("[mantel.R] error: %s", conditionMessage(e)))
  quit(status = 1)
})

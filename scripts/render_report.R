#!/usr/bin/env Rscript
# builds html reports (overview + per-tree)

suppressPackageStartupMessages({
  library(optparse)
  library(rmarkdown)
  library(rlang)
  library(yaml)
  library(purrr)
})

error_exit <- function(msg) {
  message(sprintf("[render_report] %s", msg))
  quit(status = 1)
}

#### CLI parsing #######
make_parser <- function() {
  option_list <- list(
    make_option("--project-dir", type = "character", default = "."),
    make_option("--rmd", type = "character"),
    make_option("--treermd", type = "character", default = "scripts/tree_report.Rmd"),
    make_option("--output-dir", type = "character", default = "results/report", dest = "output_dir"),
    make_option("--output-file", type = "character", default = "index.html", dest = "output_file"),
    make_option("--results", type = "character", default = "results"),
    make_option("--config", type = "character", default = "config.yaml"),
    make_option("--mapping", type = "character", default = "resources/header_map.tsv")
  )
  OptionParser(option_list = option_list)
}

opts <- parse_args(make_parser())

# anchor everything in project directory
project_dir <- opts[["project_dir"]] %||% "."
treermd_arg <- opts$treermd
project_root <- if (nzchar(project_dir)) project_dir else "." # normalize blank input

# helpers to convert command line paths into absolute
is_absolute_path <- function(path) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(path)) 
    return(FALSE)
  startsWith(path, "/") ||
    startsWith(path, "~") ||
    grepl("^[A-Za-z]:[/\\\\]", path)
}
resolve_path <- function(path, label) {
  if (is.null(path) || length(path) == 0) {
    error_exit(sprintf("argument for %s is missing.", label))
  }
  path <- path[[1]]
  if (!nzchar(path)) {
    error_exit(sprintf("argument for %s is empty.", label))
  }
  full <- if (is_absolute_path(path)) path else file.path(project_root, path)
  normalizePath(full, mustWork = FALSE)
}

# input checks
input_rmd <- resolve_path(opts$rmd, "--rmd")
if (!file.exists(input_rmd)) {
  error_exit(sprintf("Rmd file not found at: %s", input_rmd))
}
treermd <- resolve_path(treermd_arg, "--treermd")
if (!file.exists(treermd)) {
  error_exit(sprintf("Tree Rmd file not found at: %s", treermd))
}

# create output dir
output_dir <- resolve_path(opts$output_dir, "--output-dir")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# use resolve on path params
render_params <- list(
  results_dir = resolve_path(opts$results, "--results"),
  config_path = resolve_path(opts$config, "--config"),
  mapping_path = resolve_path(opts$mapping, "--mapping")
)

# check missing inputs
missing <- Filter(function(p) !file.exists(p), render_params)
if (length(missing) > 0) {
  msg <- paste(sprintf("%s", missing), collapse = ", ")
  error_exit(sprintf("Required input files missing: %s", msg))
}

# load yaml config
config_data <- tryCatch(yaml::read_yaml(render_params$config_path), error = function(e) NULL)
if (is.null(config_data)) {
  error_exit(sprintf("Failed to read config: %s", render_params$config_path))
}
tree_definitions <- config_data$trees %||% list()

# renderer helper
render_template <- function(input, output_file, params_list) {
  rmarkdown::render(
    input = input,
    output_file = output_file,
    output_dir = output_dir,
    params = params_list,
    quiet = TRUE
  )
} # calls rmarkdown::render() with given template

# overview report
tryCatch(
render_template(input_rmd, opts$output_file, render_params),
  error = function(e) {
    error_exit(sprintf("Rendering failed: %s", conditionMessage(e)))
  }
)

# per-tree reports
purrr::walk(tree_definitions, function(tree) {
  tree_name <- tree$name %||% ""
  if (!nzchar(tree_name)) return(NULL)
  tree_params <- render_params
  tree_params$tree_name <- tree_name
  tree_output <- sprintf("%s.html", tree_name)
  tryCatch(
    render_template(treermd, tree_output, tree_params),
    error = function(e) {
      error_exit(sprintf("Rendering failed for '%s': %s", tree_name, conditionMessage(e)))
    }
  )
})

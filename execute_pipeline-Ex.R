#!/usr/bin/env Rscript
library(rmarkdown)

# Render one or more R Markdown reports, in the order given, against a single
# config file. Each report renders in a fresh environment so later reports do
# not inherit earlier ones' in-memory objects (they read their inputs from disk).
#
# Usage:
#   Rscript execute_pipeline-Ex.R report1.Rmd [report2.Rmd ...]
#
# Config file:
#   Defaults to config.yaml. Override with the CONFIG_FILE environment variable:
#     CONFIG_FILE=config_test.yaml Rscript execute_pipeline-Ex.R part1.Rmd
#
# This is normally invoked through run_templates.sh (which forwards all of its
# arguments here), not called directly.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Supply one or more .Rmd report files to render, in order.", call. = FALSE)
}

config.file <- Sys.getenv("CONFIG_FILE", unset = "config.yaml")
if (!file.exists(config.file)) {
  stop("Config file not found: ", config.file,
       " (set CONFIG_FILE to override the default config.yaml).", call. = FALSE)
}

render_report <- function(template.filename, config.file) {
  if (!file.exists(template.filename)) {
    stop("Report file not found: ", template.filename, call. = FALSE)
  }
  output.name <- gsub("\\.Rmd$", ".pdf", template.filename)
  message("=================================================================")
  message("Rendering ", template.filename, "  (config: ", config.file, ")")
  message("=================================================================")
  rmarkdown::render(
    template.filename,
    params      = list(config.args = config.file),
    envir       = new.env(parent = globalenv()),
    output_file = output.name
  )
}

# Render sequentially. render() stops on error, so if an earlier report fails
# (e.g. part 2a), the later ones (part 2b) are never run.
for (rmd in args) render_report(rmd, config.file)

message("All reports rendered.")

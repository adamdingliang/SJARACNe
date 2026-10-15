#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(NetBID2))

demo_root <- system.file("demo1/network", package = "NetBID2")
network_files <- list.files(
  demo_root,
  pattern = "consensus_network_ncol_[.]txt$",
  recursive = TRUE,
  full.names = TRUE
)
tf_network_file <- network_files[grepl("output_tf_", network_files)]

if (length(tf_network_file) != 1L) {
  stop("Could not identify the bundled NetBID2 TF consensus network")
}

network <- get.SJAracne.network(tf_network_file)
if (!is.list(network) || !inherits(network$igraph_obj, "igraph")) {
  stop("get.SJAracne.network did not produce the expected igraph object")
}

output_directory <- tempfile("netbid2-qc-smoke-")
dir.create(output_directory)

result <- draw.network.QC(
  network$igraph_obj,
  outdir = output_directory,
  prefix = "TF_smoke_",
  directed = TRUE,
  weighted = FALSE,
  generate_html = TRUE,
  html_info_limit = TRUE
)

html_files <- list.files(output_directory, pattern = "[.]html$", full.names = TRUE)
if (!isTRUE(result) || length(html_files) == 0L) {
  stop("draw.network.QC did not produce its HTML report")
}

cat("NetBID2 QC smoke test passed\n")
cat("Input: ", tf_network_file, "\n", sep = "")
cat("Vertices: ", igraph::vcount(network$igraph_obj), "\n", sep = "")
cat("Edges: ", igraph::ecount(network$igraph_obj), "\n", sep = "")
cat("Temporary report: ", html_files[[1L]], "\n", sep = "")

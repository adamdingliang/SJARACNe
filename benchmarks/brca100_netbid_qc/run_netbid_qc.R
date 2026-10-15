#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(
    "Usage: run_netbid_qc.R CONSENSUS_NCOL DRIVER_LIST OUTPUT_DIR PREFIX",
    call. = FALSE
  )
}

network_file <- normalizePath(args[[1L]], mustWork = TRUE)
driver_file <- normalizePath(args[[2L]], mustWork = TRUE)
output_directory <- args[[3L]]
prefix <- args[[4L]]

suppressPackageStartupMessages(library(NetBID2))

if (dir.exists(output_directory)) {
  stop("QC output directory already exists: ", output_directory)
}
dir.create(output_directory, recursive = TRUE)

network <- get.SJAracne.network(network_file)
if (!is.list(network) || !inherits(network$igraph_obj, "igraph")) {
  stop("get.SJAracne.network did not return the expected network object")
}
if (igraph::ecount(network$igraph_obj) != nrow(network$network_dat)) {
  stop("NetBID2 edge count does not match the consensus table")
}

result <- draw.network.QC(
  network$igraph_obj,
  outdir = output_directory,
  prefix = prefix,
  directed = TRUE,
  weighted = FALSE,
  generate_html = TRUE,
  html_info_limit = TRUE
)
if (!isTRUE(result)) {
  stop("draw.network.QC returned a non-success result")
}

html_file <- file.path(output_directory, paste0(prefix, "netQC.html"))
if (!file.exists(html_file) || file.info(html_file)$size == 0) {
  stop("NetBID2 did not create its expected HTML report: ", html_file)
}

drivers <- trimws(readLines(driver_file, warn = FALSE))
drivers <- drivers[nzchar(drivers)]
if (anyDuplicated(drivers)) {
  stop("Driver list contains duplicate accessions")
}

graph <- network$igraph_obj
out_degree <- igraph::degree(graph, mode = "out")
driver_degree <- stats::setNames(rep.int(0L, length(drivers)), drivers)
shared <- intersect(names(out_degree), drivers)
driver_degree[shared] <- as.integer(out_degree[shared])
active_degree <- driver_degree[driver_degree > 0L]

components <- igraph::components(graph, mode = "weak")
degree_distribution <- igraph::degree_distribution(graph)
degree_table <- data.frame(
  k = seq_len(max(igraph::degree(graph))),
  pk = degree_distribution[-1L] + 1 / igraph::vcount(graph)
)
scale_free_model <- stats::lm(log10(pk) ~ log10(k), data = degree_table)
scale_free_adjusted_r2 <- summary(scale_free_model)$adj.r.squared

quantile_value <- function(values, probability) {
  as.numeric(stats::quantile(values, probability, names = FALSE, type = 7))
}

summary_table <- data.frame(
  metric = c(
    "candidate_drivers", "active_drivers", "edges", "incident_nodes",
    "weak_components", "largest_weak_component", "density",
    "target_size_zero_mean", "target_size_zero_median",
    "target_size_zero_q25", "target_size_zero_q75", "target_size_zero_max",
    "target_size_active_mean", "target_size_active_median",
    "target_size_active_q25", "target_size_active_q75",
    "target_size_active_max", "scale_free_adjusted_r2"
  ),
  value = c(
    length(drivers), length(active_degree), igraph::ecount(graph),
    igraph::vcount(graph), components$no, max(components$csize),
    igraph::edge_density(graph, loops = FALSE), mean(driver_degree),
    stats::median(driver_degree), quantile_value(driver_degree, 0.25),
    quantile_value(driver_degree, 0.75), max(driver_degree),
    mean(active_degree), stats::median(active_degree),
    quantile_value(active_degree, 0.25), quantile_value(active_degree, 0.75),
    max(active_degree), scale_free_adjusted_r2
  )
)

utils::write.table(
  summary_table,
  file.path(output_directory, "network_summary.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
utils::write.table(
  data.frame(driver = names(driver_degree), target_count = as.integer(driver_degree)),
  file.path(output_directory, "driver_target_sizes.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("NetBID2 QC complete\n")
cat("Network: ", network_file, "\n", sep = "")
cat("Drivers: ", length(drivers), "\n", sep = "")
cat("Active drivers: ", length(active_degree), "\n", sep = "")
cat("Edges: ", igraph::ecount(graph), "\n", sep = "")
cat("HTML: ", html_file, "\n", sep = "")

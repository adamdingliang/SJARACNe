#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(NetBID2))

package_field <- function(description, field) {
  value <- description[[field]]
  if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) {
    return("")
  }
  as.character(value[[1L]])
}

environment_table <- function() {
  description <- utils::packageDescription("NetBID2")
  data.frame(
    component = c("R", "NetBID2", "NetBID2_remote_sha", "igraph"),
    version = c(
      R.version.string,
      package_field(description, "Version"),
      package_field(description, "RemoteSha"),
      as.character(utils::packageVersion("igraph"))
    ),
    stringsAsFactors = FALSE
  )
}

if (length(args) == 1L && identical(args[[1L]], "--probe")) {
  utils::write.table(
    environment_table(), stdout(), sep = "\t", row.names = FALSE,
    col.names = TRUE, quote = FALSE, na = "NA"
  )
  quit(save = "no", status = 0L)
}

if (length(args) != 5L) {
  stop(
    paste(
      "Usage: run_netbid_qc.R CONSENSUS_NCOL DRIVER_LIST",
      "OUTPUT_DIR PREFIX GENERATE_HTML"
    ),
    call. = FALSE
  )
}

network_file <- normalizePath(args[[1L]], mustWork = TRUE)
driver_file <- normalizePath(args[[2L]], mustWork = TRUE)
output_directory <- args[[3L]]
prefix <- args[[4L]]
html_token <- tolower(args[[5L]])

if (!grepl("^[A-Za-z0-9_.-]+$", prefix)) {
  stop("PREFIX contains unsafe path characters", call. = FALSE)
}
if (!html_token %in% c("true", "false")) {
  stop("GENERATE_HTML must be true or false", call. = FALSE)
}
generate_html <- identical(html_token, "true")

if (dir.exists(output_directory) || file.exists(output_directory)) {
  stop("QC output path already exists: ", output_directory)
}
dir.create(output_directory, recursive = TRUE)

network <- get.SJAracne.network(network_file)
if (!is.list(network) || !inherits(network$igraph_obj, "igraph")) {
  stop("get.SJAracne.network did not return the expected network object")
}
if (!is.data.frame(network$network_dat)) {
  stop("get.SJAracne.network did not return network_dat as a data frame")
}
if (igraph::ecount(network$igraph_obj) != nrow(network$network_dat)) {
  stop("NetBID2 edge count does not match the consensus table")
}

drivers <- trimws(readLines(driver_file, warn = FALSE))
drivers <- drivers[nzchar(drivers)]
if (length(drivers) == 0L) {
  stop("Driver list is empty")
}
if (anyDuplicated(drivers)) {
  stop("Driver list contains duplicate accessions")
}

graph <- network$igraph_obj
out_degree <- igraph::degree(graph, mode = "out")
driver_degree <- stats::setNames(rep.int(0L, length(drivers)), drivers)
shared <- intersect(names(out_degree), drivers)
driver_degree[shared] <- as.integer(out_degree[shared])
active_degree <- driver_degree[driver_degree > 0L]

if (igraph::vcount(graph) > 0L) {
  components <- igraph::components(graph, mode = "weak")
  weak_component_count <- components$no
  largest_weak_component <- max(components$csize)
  largest_weak_component_fraction <- largest_weak_component / igraph::vcount(graph)
} else {
  weak_component_count <- 0L
  largest_weak_component <- 0L
  largest_weak_component_fraction <- NA_real_
}

scale_free_adjusted_r2 <- NA_real_
if (igraph::vcount(graph) > 0L) {
  all_degree <- igraph::degree(graph)
  max_degree <- max(all_degree)
  if (max_degree >= 3L) {
    degree_distribution <- igraph::degree_distribution(graph)
    degree_table <- data.frame(
      k = seq_len(max_degree),
      pk = degree_distribution[-1L] + 1 / igraph::vcount(graph)
    )
    scale_free_model <- stats::lm(log10(pk) ~ log10(k), data = degree_table)
    candidate_r2 <- summary(scale_free_model)$adj.r.squared
    if (is.finite(candidate_r2)) {
      scale_free_adjusted_r2 <- candidate_r2
    }
  }
}

quantile_value <- function(values, probability) {
  if (length(values) == 0L) {
    return(NA_real_)
  }
  as.numeric(stats::quantile(values, probability, names = FALSE, type = 7))
}

mean_value <- function(values) {
  if (length(values) == 0L) NA_real_ else mean(values)
}

median_value <- function(values) {
  if (length(values) == 0L) NA_real_ else stats::median(values)
}

max_value <- function(values) {
  if (length(values) == 0L) NA_real_ else max(values)
}

summary_table <- data.frame(
  metric = c(
    "candidate_drivers", "active_drivers", "active_driver_fraction",
    "edges", "incident_nodes", "weak_components",
    "largest_weak_component", "largest_weak_component_fraction", "density",
    "target_size_zero_mean", "target_size_zero_median",
    "target_size_zero_q25", "target_size_zero_q75", "target_size_zero_max",
    "target_size_active_mean", "target_size_active_median",
    "target_size_active_q25", "target_size_active_q75",
    "target_size_active_max", "scale_free_adjusted_r2"
  ),
  value = c(
    length(drivers), length(active_degree), length(active_degree) / length(drivers),
    igraph::ecount(graph), igraph::vcount(graph), weak_component_count,
    largest_weak_component, largest_weak_component_fraction,
    if (igraph::vcount(graph) > 1L) {
      igraph::edge_density(graph, loops = FALSE)
    } else {
      NA_real_
    },
    mean_value(driver_degree), median_value(driver_degree),
    quantile_value(driver_degree, 0.25), quantile_value(driver_degree, 0.75),
    max_value(driver_degree), mean_value(active_degree),
    median_value(active_degree), quantile_value(active_degree, 0.25),
    quantile_value(active_degree, 0.75), max_value(active_degree),
    scale_free_adjusted_r2
  ),
  stringsAsFactors = FALSE
)

html_file <- ""
if (generate_html) {
  result <- draw.network.QC(
    graph,
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
}

utils::write.table(
  summary_table,
  file.path(output_directory, "network_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
utils::write.table(
  data.frame(
    driver = names(driver_degree),
    target_count = as.integer(driver_degree),
    stringsAsFactors = FALSE
  ),
  file.path(output_directory, "driver_target_sizes.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)
utils::write.table(
  environment_table(),
  file.path(output_directory, "netbid_environment.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE, na = "NA"
)

cat("NetBID2 QC complete\n")
cat("Network: ", network_file, "\n", sep = "")
cat("Drivers: ", length(drivers), "\n", sep = "")
cat("Active drivers: ", length(active_degree), "\n", sep = "")
cat("Edges: ", igraph::ecount(graph), "\n", sep = "")
cat("Generate HTML: ", generate_html, "\n", sep = "")
if (generate_html) {
  cat("HTML: ", html_file, "\n", sep = "")
}

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(NetBID2))

description <- utils::packageDescription("NetBID2")
required_exports <- c("get.SJAracne.network", "draw.network.QC")
missing_exports <- setdiff(required_exports, getNamespaceExports("NetBID2"))

if (length(missing_exports) > 0L) {
  stop("Missing NetBID2 exports: ", paste(missing_exports, collapse = ", "))
}

cat("R version: ", R.version.string, "\n", sep = "")
cat("NetBID2 version: ", as.character(description$Version), "\n", sep = "")
cat("NetBID2 commit: ", as.character(description$RemoteSha), "\n", sep = "")
cat("Verified exports: ", paste(required_exports, collapse = ", "), "\n", sep = "")

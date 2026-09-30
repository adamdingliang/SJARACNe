#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(NetBID2))

description <- utils::packageDescription("NetBID2")
lines <- c(
  paste0("NetBID2 version: ", as.character(description$Version)),
  paste0(
    "NetBID2 remote: ", as.character(description$RemoteUsername), "/",
    as.character(description$RemoteRepo)
  ),
  paste0("NetBID2 commit: ", as.character(description$RemoteSha)),
  "",
  capture.output(utils::sessionInfo())
)
cat(sub("[[:space:]]+$", "", lines), sep = "\n")

#!/usr/bin/env Rscript

netbid_repository <- "jyyulab/NetBID"
netbid_commit <- "5defa454d600b94f5dd6d1f9f4428f99759a6821"
netbid_version <- "2.2.0"

installed_commit <- tryCatch(
  utils::packageDescription("NetBID2", fields = "RemoteSha"),
  error = function(e) NA_character_
)
installed_version <- tryCatch(
  as.character(utils::packageVersion("NetBID2")),
  error = function(e) NA_character_
)

if (!identical(installed_commit, netbid_commit) ||
    !identical(installed_version, netbid_version)) {
  remotes::install_github(
    paste0(netbid_repository, "@", netbid_commit),
    dependencies = FALSE,
    upgrade = "never",
    build_vignettes = FALSE
  )
}

description <- utils::packageDescription("NetBID2")
stopifnot(
  identical(as.character(description$Version), netbid_version),
  identical(as.character(description$RemoteSha), netbid_commit)
)

message(
  "Installed NetBID2 ", description$Version,
  " from commit ", description$RemoteSha
)

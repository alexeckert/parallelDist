#!/usr/bin/env Rscript

if (!requireNamespace("covr", quietly = TRUE)) {
  stop("The covr package must be installed to run coverage", call. = FALSE)
}

wrapper <- normalizePath(file.path("tools", "gcov-wrapper.sh"), winslash = "/", mustWork = TRUE)

old_opts <- options(covr.gcov = wrapper)
on.exit(options(old_opts), add = TRUE)

token <- Sys.getenv("CODECOV_TOKEN", "")

if (nzchar(token)) {
  covr::codecov(token = token)
} else {
  covr::codecov()
}

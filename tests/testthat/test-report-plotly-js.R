library(testthat)

get_report_fn <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    get(name, envir = globalenv())
  } else if (requireNamespace("TALOS", quietly = TRUE)) {
    getFromNamespace(name, "TALOS")
  } else {
    stop("Function ", name, " not found. Load the package with devtools::load_all().")
  }
}

test_that("local Plotly bundle discovery never requires the plotly package", {
  locate <- get_report_fn(".local_plotly_js")
  old <- options(TALOS.plotly_js = NULL)
  on.exit(options(old), add = TRUE)
  path <- locate()
  expect_length(path, 1L)
  expect_true(is.character(path))
  expect_true(is.na(path) || (file.exists(path) && file.access(path, 4L) == 0L))
})

test_that("an explicitly configured Plotly bundle is used", {
  locate <- get_report_fn(".local_plotly_js")
  old <- options(TALOS.plotly_js = NULL)
  on.exit(options(old), add = TRUE)
  f <- tempfile(fileext = ".js")
  writeLines("window.Plotly = window.Plotly || {};", f)
  on.exit(unlink(f), add = TRUE)
  options(TALOS.plotly_js = f)
  expect_identical(normalizePath(locate()), normalizePath(f))
})

test_that("HTML Plotly script selection falls back to CDN instead of erroring", {
  script_tag <- get_report_fn(".plotly_script_tag")
  old <- options(TALOS.plotly_js = NULL)
  on.exit(options(old), add = TRUE)
  options(TALOS.plotly_js = tempfile("missing-plotly-"))
  out <- script_tag("auto")
  expect_false(out$inline)
  expect_match(out$tag, "cdn.jsdelivr.net")
})

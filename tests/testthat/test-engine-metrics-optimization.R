library(testthat)

get_fn_opt <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    get(name, envir = globalenv())
  } else if (requireNamespace("TALOS", quietly = TRUE)) {
    getFromNamespace(name, "TALOS")
  } else {
    stop("Function ", name, " not found. Load the package with devtools::load_all().")
  }
}

test_that("indel-aware vectorized junction screening preserves orientation", {
  find <- get_fn_opt(".find_tandem_junction_candidates")
  template <- "ACGTTGCAACCTGATCGTAGGCTAACGATCGTACCTGAGCTTACGATCGGA"
  flank <- 12L
  junction <- paste0(substr(template, nchar(template) - flank + 1L, nchar(template)),
                     substr(template, 1L, flank))
  forward <- paste0("TTTT", junction, "AAAA")
  reverse <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(forward)))
  got <- find(c(forward, "NNNN", NA_character_, reverse, "ZZZ",
                "ACGTACGTACGTACGTACGTACGT"), template,
              anchor_width = flank, max_mismatch = 0L)
  expect_equal(got, c(1L, 4L))
  expect_identical(attr(got, "orientation"), c("forward", "reverse"))
})

test_that("template-derived indexes and junction windows are cached", {
  cache_get <- get_fn_opt(".talos_cache_get")
  seed_index <- get_fn_opt(".template_unique_seed_map")
  junction <- get_fn_opt(".template_junction_window")
  template <- "ACGTTGCAACCTGATCGTAGGCTAACGATCGTACCTGAGCTTACGATCGGA"
  a <- seed_index(template, 9L)
  expect_identical(seed_index(template, 9L), a)
  expect_false(is.null(cache_get(paste0("seed|9|", nchar(template), "|", template))))
  j <- junction(template, 12L)
  expect_identical(junction(template, 12L), j)
  expect_false(is.null(cache_get(paste0("junction|12|", nchar(template), "|", template))))
})

test_that("ambiguity coverage is reconstructed from per-read positional evidence", {
  align <- get_fn_opt(".align_reads_to_itd_template")
  ambiguity <- get_fn_opt(".align_reads_to_itd_template_ambiguity")
  unit <- "ATCCGTTCCTAATAAGGAATGGTGATTCCC"
  template <- paste0(unit, unit)
  reads <- c(unit, substr(unit, 1L, 22L))
  best <- align(reads, template, min_pid = 95, min_aligned = 20L)
  got <- ambiguity(reads, template, min_pid = 95, min_aligned = 20L,
                   best_result = best, max_alternatives = 4L, score_delta = 0)
  expect_length(got$per_read_certain_positions, length(reads))
  expect_length(got$per_read_possible_positions, length(reads))
  expected_certain <- tabulate(unlist(lapply(got$per_read_certain_positions, unique),
                                      use.names = FALSE), nbins = nchar(template))
  expected_possible <- tabulate(unlist(lapply(got$per_read_possible_positions, unique),
                                       use.names = FALSE), nbins = nchar(template))
  expect_equal(got$certain_cov, expected_certain)
  expect_equal(got$possible_cov, expected_possible)
})

test_that("HTML error diagnostics include context, class and call stack", {
  diagnostic <- get_fn_opt(".html_error_diagnostic")
  err <- tryCatch(stop("synthetic report error"), error = identity)
  msg <- diagnostic(err, context = "HTML report test")
  expect_match(msg, "synthetic report error", fixed = TRUE)
  expect_match(msg, "simpleError", fixed = TRUE)
  expect_match(msg, "stack=", fixed = TRUE)
})


test_that("junction screening accepts small insertions and deletions", {
  find <- get_fn_opt(".find_tandem_junction_candidates")
  template <- "ACGTTGCAACCTGATCGTAGGCTAACGATCGTACCTGAGCTTACGATCGGA"
  flank <- 12L
  junction <- paste0(substr(template, nchar(template) - flank + 1L, nchar(template)),
                     substr(template, 1L, flank))
  # Insert one base within the crossing window and delete one base in a
  # separate read. Both must still align across almost the complete window.
  insertion <- paste0(substr(junction, 1L, 10L), "A", substr(junction, 11L, nchar(junction)))
  deletion <- paste0(substr(junction, 1L, 10L), substr(junction, 12L, nchar(junction)))
  got <- find(c(paste0("TT", insertion, "GG"),
                paste0("TT", deletion, "GG")), template,
              anchor_width = flank, max_mismatch = 2L)
  expect_equal(got, c(1L, 2L))
  expect_identical(attr(got, "orientation"), c("forward", "forward"))
})

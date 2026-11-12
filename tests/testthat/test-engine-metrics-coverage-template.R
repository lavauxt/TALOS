library(testthat)

get_fn <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    get(name, envir = globalenv())
  } else if (requireNamespace("TALOS", quietly = TRUE)) {
    getFromNamespace(name, "TALOS")
  } else {
    stop("Function ", name, " not found. Load the package with devtools::load_all().")
  }
}

test_that("unresolved ITD sequence gets a reference-only coverage template", {
  resolve <- get_fn(".resolve_itd_coverage_template")
  ref <- paste(rep(c("A", "C", "G", "T"), 100L), collapse = "")

  got <- resolve(
    itd_seq = NA_character_,
    best_len = 300L,
    ref_dna = ref,
    local_seg_start = 25L,
    is_ptd_scale = FALSE,
    seg_in_ref = TRUE,
    sequence_imputed = FALSE
  )

  expect_identical(got$sequence, substr(ref, 25L, 324L))
  expect_identical(got$source, "reference_segment_for_coverage_only")
})

test_that("reference-only coverage template is refused for PTD-scale or truncated reference", {
  resolve <- get_fn(".resolve_itd_coverage_template")
  ref <- paste(rep("ACGT", 100L), collapse = "")

  ptd <- resolve(
    itd_seq = NA_character_, best_len = 300L, ref_dna = ref,
    local_seg_start = 25L, is_ptd_scale = TRUE, seg_in_ref = TRUE
  )
  truncated <- resolve(
    itd_seq = NA_character_, best_len = 300L, ref_dna = ref,
    local_seg_start = 250L, is_ptd_scale = FALSE, seg_in_ref = TRUE
  )

  expect_true(is.na(ptd$sequence))
  expect_identical(ptd$source, "unavailable")
  expect_true(is.na(truncated$sequence))
  expect_identical(truncated$source, "unavailable")
})

test_that("read-derived and fully imputed sequences retain their source labels", {
  resolve <- get_fn(".resolve_itd_coverage_template")
  read_seq <- paste(rep("A", 300L), collapse = "")

  observed <- resolve(read_seq, 300L, "CCCC", 1L, FALSE, FALSE, FALSE)
  imputed <- resolve(read_seq, 300L, "CCCC", 1L, FALSE, FALSE, TRUE)

  expect_identical(observed$sequence, read_seq)
  expect_identical(observed$source, "read_derived_sequence")
  expect_identical(imputed$source, "reference_imputed_sequence")
})


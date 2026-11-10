test_that("PTD graph path scoring uses the minimum cost per exon transition", {
  edges <- data.frame(
    from = c("Exon_1", "Exon_1", "Exon_2", "Exon_3"),
    to = c("Exon_2", "Exon_2", "Exon_3", "Exon_1"),
    prob = c(exp(-2), exp(-1), exp(-3), 0.5),
    weight = c(2, 1, 3, 0),
    support_count = c(0L, 0L, 0L, 2L),
    qnames = c("", "", "", "r1;r2"),
    evidence_type = c("canonical", "split_read", "canonical", "split_read"),
    dup_start = c(NA_integer_, NA_integer_, NA_integer_, 101L),
    dup_end = c(NA_integer_, NA_integer_, NA_integer_, 200L),
    n_split = c(NA_integer_, NA_integer_, NA_integer_, 2L),
    n_pairs = c(NA_integer_, NA_integer_, NA_integer_, 0L),
    stringsAsFactors = FALSE
  )
  cycles <- .detect_ptd_cycles(edges, c("Exon_1", "Exon_2", "Exon_3"))
  expect_length(cycles, 1L)
  expect_equal(cycles[[1L]]$nll_score, 4 - log(0.5 + 1e-9), tolerance = 1e-8)
})

test_that("soft-clip filtering keeps MAPQ weights aligned", {
  got <- .filter_weighted_sequences(c("AAA", "", "TTT", NA_character_), c(20, 60, 40, 50))
  expect_equal(got$seqs, c("AAA", "TTT"))
  expect_equal(got$weights, c(20, 40))
  expect_error(.filter_weighted_sequences(c("AAA"), c(20, 30)))
})

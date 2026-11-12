test_that("all-position k-mer matching preserves repeated reference hits", {
  hits <- TALOS:::.kmer_match_positions(
    c("AAA", "CCC", "TTT"),
    c("AAA", "CCC", "AAA", "GGG")
  )
  expect_identical(hits[[1]], c(1L, 3L))
  expect_identical(hits[[2]], 2L)
  expect_identical(hits[[3]], integer(0))
})

test_that("PTD merge gap is inclusive", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr11", "chr11"),
    ranges = IRanges::IRanges(start = c(100L, 121L), end = c(109L, 130L))
  )
  # The number of uncovered bases between 100-109 and 121-130 is 11.
  expect_equal(length(GenomicRanges::reduce(gr, min.gapwidth = .ptd_merge_min_gapwidth(10L))), 2L)
  expect_equal(length(GenomicRanges::reduce(gr, min.gapwidth = .ptd_merge_min_gapwidth(11L))), 1L)
  expect_equal(.ptd_merge_min_gapwidth(0L), 1L)
  expect_error(.ptd_merge_min_gapwidth(-1L))
})

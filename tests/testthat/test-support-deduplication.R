test_that("read-level support metrics use one row per QNAME", {
  rows <- data.frame(
    read_name = c("read1", "read1", "read2", NA_character_, ""),
    breakpoint = c(110, 101, 100, 100, 100),
    mapq = c(60, 20, 40, 50, 50),
    itd_seq = c("AAAA", "CCCC", "GGGG", "TTTT", "AAAA"),
    stringsAsFactors = FALSE
  )
  got <- .deduplicate_support_rows(rows, target_bp = 100)
  expect_equal(got$read_name, c("read1", "read2"))
  expect_equal(got$breakpoint, c(101, 100))
})

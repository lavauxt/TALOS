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


test_that("coverage read deduplication keeps both mates but collapses duplicate alignments", {
  dedup <- .deduplicate_coverage_read_rows
  rows <- data.frame(
    read_name = c("frag1", "frag1", "frag1", "frag2"),
    read_seq = c("AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC",
                 "AAAAAAAAAAAAAAAAAAAA", "GGGGGGGGGGGGGGGGGGGG"),
    read_start = c(100L, 180L, 100L, 200L),
    read_end = c(119L, 199L, 119L, 219L),
    cigar = c("20M", "20M", "20M", "20M"),
    is_reverse = c(FALSE, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
  got <- dedup(rows)
  expect_equal(nrow(got), 3L)
  expect_equal(sum(got$read_name == "frag1"), 2L)
  expect_setequal(got$read_seq, c("AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC",
                                  "GGGGGGGGGGGGGGGGGGGG"))
})

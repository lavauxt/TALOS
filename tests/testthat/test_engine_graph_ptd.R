test_that("PTD graph identifies transcript-order back edges", {
  edges <- data.frame(
    from = c("Exon_1", "Exon_2", "Exon_3", "Exon_3"),
    to = c("Exon_2", "Exon_3", "Exon_1", "Exon_2"),
    weight = c(0.01, 0.01, 0.2, 0.3),
    support_count = c(100L, 100L, 4L, 2L),
    qnames = c("", "", "read3", "read2"),
    evidence_type = c(
      "canonical",
      "canonical",
      "split_read",
      "discordant_pair"
    ),
    stringsAsFactors = FALSE
  )

  cycles <- .detect_ptd_cycles(edges, c("Exon_1", "Exon_2", "Exon_3"))

  expect_length(cycles, 2L)
  expect_setequal(
    vapply(cycles, `[[`, character(1L), "support_qnames"),
    c("read3", "read2")
  )
})


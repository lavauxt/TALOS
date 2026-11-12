test_that("variant hover uses current report fields and tolerates missing optional fields", {
  row <- data.frame(
    Sample = "HORIZON", Gene = "FLT3", GenomicPosition = 28608047,
    Length = 300, AlleleFrequency = 0.12,
    ITDReadCoverage_Direct = 81, ITDReadCoverage_Combined = 95,
    RefMatch_Observed = 100, stringsAsFactors = FALSE
  )
  hover <- TALOS:::.variant_hover(row, 1L, "12.0%")
  expect_type(hover, "character")
  expect_match(hover, "RefMatch: 100.0%")
  expect_match(hover, "81.0% \\(combined 95.0%\\)")
  expect_match(hover, "Independent QNAMEs: N/A")
})

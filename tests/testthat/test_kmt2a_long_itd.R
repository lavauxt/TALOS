test_that("unique anchor matching detects overlapping repeats", {
  expect_equal(.unique_fixed_match("AAA", "AAAAA"), c(1L, 2L))
  expect_equal(.unique_fixed_match("ACG", "TTACGTT"), 3L)
  expect_length(.unique_fixed_match("CCC", "TTACGTT"), 0L)
})

test_that("anchor-refined lengths follow the configured long-ITD policy", {
  expect_equal(
    .apply_itd_length_limit(30000L, "anchor_ext", 30000L, TRUE),
    list(length = 30000L, type = "anchor_ext", converted = FALSE)
  )
  expect_equal(
    .apply_itd_length_limit(30001L, "anchor_ext", 30000L, TRUE),
    list(length = 0L, type = "ptd_clip", converted = TRUE)
  )
  expect_null(.apply_itd_length_limit(30001L, "anchor_ext", 30000L, FALSE))
})

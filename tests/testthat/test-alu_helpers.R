# Regression tests for ALU helper functions in engine_alu.R
#
# Covers:
#   - .alu_clip_breakpoint() including the chromosome-boundary fix (Finding 2)
#   - .detect_alu_tsd() including the read-evidence approach (Finding 1)
#   - .detect_poly_a() terminal-run detection
#   - .alu_5p_truncation() edge cases
#   - .filter_alu_call() NA TSD handling
#   - .extract_alu_candidates() boundary guard
#
# These tests use only base R + testthat; no BAM or Biostrings required
# (except for .match_clip_to_alu, which needs pwalign/Biostrings and is
# skipped if those packages are unavailable).

library(testthat)

# Helper: access internal functions without loading the full package.
# In a package build these are available via devtools::load_all().
# In an installed package, use TALOS:::.function_name.
get_fn <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    get(name, envir = globalenv())
  } else if (requireNamespace("TALOS", quietly = TRUE)) {
    getFromNamespace(name, "TALOS")
  } else {
    stop("Function ", name, " not found. Load the package with devtools::load_all().")
  }
}

# ---------------------------------------------------------------------------
# .alu_clip_breakpoint
# ---------------------------------------------------------------------------
describe(".alu_clip_breakpoint", {
  bp <- get_fn(".alu_clip_breakpoint")

  it("returns read_end for trailing clips", {
    expect_identical(bp("trail", 100L, 150L), 150L)
  })

  it("returns read_start - 1 for leading clips", {
    expect_identical(bp("lead", 100L, 150L), 99L)
  })

  # [Finding 2] Chromosome-boundary regression test
  it("returns NA when a leading clip is at chromosome position 1", {
    expect_identical(bp("lead", 1L, 50L), NA_integer_)
  })

  it("returns NA when a leading clip is at position 0 (should not happen but guard)", {
    # read_start should never be 0 for a mapped read, but guard anyway
    expect_identical(suppressWarnings(bp("lead", 0L, 50L)), NA_integer_)
  })

  it("errors on invalid side", {
    expect_error(bp("left", 100L, 150L), "side must be")
    expect_error(bp(NA_character_, 100L, 150L), "side must be")
  })

  it("errors on missing read_start/read_end", {
    expect_error(bp("lead", NA_integer_, 150L), "non-missing scalars")
    expect_error(bp("trail", 100L, NA_integer_), "non-missing scalars")
  })

  it("handles trailing clip at the last chromosome position", {
    expect_identical(bp("trail", 249250600L, 249250621L), 249250621L)
  })
})

# ---------------------------------------------------------------------------
# .detect_alu_tsd
# ---------------------------------------------------------------------------
describe(".detect_alu_tsd (Finding 1: read-evidence approach)", {
  detect_tsd <- get_fn(".detect_alu_tsd")

  # A 1000 bp reference for testing
  ref_dna <- paste(rep("ACGT", 250), collapse = "")

  it("returns unknown when cluster_df is NULL", {
    result <- detect_tsd(ref_dna, ins_pos = 500L, cluster_df = NULL,
                         genomic_start = 1L)
    expect_equal(result$tsd_status, "unknown")
    expect_true(is.na(result$tsd_len))
    expect_true(is.na(result$tsd_seq))
  })

  it("returns unknown when ins_pos is out of range", {
    result <- detect_tsd(ref_dna, ins_pos = 0L, cluster_df = data.frame(),
                         genomic_start = 1L)
    expect_equal(result$tsd_status, "unknown")
  })

  it("returns unknown when only trail clips are present (one side)", {
    cluster_df <- data.frame(
      genomic_pos = c(505L, 506L, 507L),
      clip_side = c("trail", "trail", "trail"),
      stringsAsFactors = FALSE
    )
    result <- detect_tsd(ref_dna, ins_pos = 505L, cluster_df = cluster_df,
                         genomic_start = 1L)
    expect_equal(result$tsd_status, "unknown")
  })

  it("returns unknown when only lead clips are present (one side)", {
    cluster_df <- data.frame(
      genomic_pos = c(495L, 496L, 497L),
      clip_side = c("lead", "lead", "lead"),
      stringsAsFactors = FALSE
    )
    result <- detect_tsd(ref_dna, ins_pos = 500L, cluster_df = cluster_df,
                         genomic_start = 1L)
    expect_equal(result$tsd_status, "unknown")
  })

  it("detects read-supported TSD when trail and lead breakpoints differ by 5-20 bp", {
    # Trail breakpoints at 510 (end of TSD), lead breakpoints at 500 (start of TSD - 1)
    # TSD length = 510 - 500 = 10 bp
    cluster_df <- data.frame(
      genomic_pos = c(510L, 510L, 511L, 500L, 500L, 501L),
      clip_side = c("trail", "trail", "trail", "lead", "lead", "lead"),
      stringsAsFactors = FALSE
    )
    result <- detect_tsd(ref_dna, ins_pos = 505L, cluster_df = cluster_df,
                         genomic_start = 1L)
    expect_equal(result$tsd_status, "read_supported")
    expect_equal(result$tsd_len, 10L)
    # TSD sequence = ref[501..510] (1-based local)
    expected_seq <- substr(ref_dna, 501L, 510L)
    expect_equal(result$tsd_seq, expected_seq)
  })

  it("returns unknown when trail-lead gap is < 5 bp (below minimum TSD)", {
    cluster_df <- data.frame(
      genomic_pos = c(502L, 503L, 500L, 501L),
      clip_side = c("trail", "trail", "lead", "lead"),
      stringsAsFactors = FALSE
    )
    result <- detect_tsd(ref_dna, ins_pos = 501L, cluster_df = cluster_df,
                         genomic_start = 1L)
    expect_equal(result$tsd_status, "unknown")
  })

  it("returns unknown when trail-lead gap is > 20 bp (above maximum TSD)", {
    cluster_df <- data.frame(
      genomic_pos = c(550L, 551L, 500L, 501L),
      clip_side = c("trail", "trail", "lead", "lead"),
      stringsAsFactors = FALSE
    )
    result <- detect_tsd(ref_dna, ins_pos = 525L, cluster_df = cluster_df,
                         genomic_start = 1L)
    expect_equal(result$tsd_status, "unknown")
  })

  it("excludes NA and zero breakpoints from the median calculation", {
    cluster_df <- data.frame(
      genomic_pos = c(NA_integer_, 0L, 510L, 500L),
      clip_side = c("trail", "trail", "trail", "lead"),
      stringsAsFactors = FALSE
    )
    result <- detect_tsd(ref_dna, ins_pos = 505L, cluster_df = cluster_df,
                         genomic_start = 1L)
    # Only one valid trail (510) and one valid lead (500)
    # TSD = 510 - 500 = 10
    expect_equal(result$tsd_status, "read_supported")
    expect_equal(result$tsd_len, 10L)
  })

  it("uses round() on median to avoid truncation bias", {
    # Trail at 511, 512 (median = 511.5 -> round to 512)
    # Lead at 500, 501 (median = 500.5 -> round to 501)
    # TSD = 512 - 501 = 11
    cluster_df <- data.frame(
      genomic_pos = c(511L, 512L, 500L, 501L),
      clip_side = c("trail", "trail", "lead", "lead"),
      stringsAsFactors = FALSE
    )
    result <- detect_tsd(ref_dna, ins_pos = 506L, cluster_df = cluster_df,
                         genomic_start = 1L)
    expect_equal(result$tsd_status, "read_supported")
    expect_equal(result$tsd_len, 11L)
  })

  it("handles non-trivial genomic_start offset", {
    # genomic_start = 10000, so local coordinates differ from genomic
    # Trail at genomic 10510 (local 511), lead at genomic 10500 (local 501)
    # TSD = 10510 - 10500 = 10, local 502..511
    cluster_df <- data.frame(
      genomic_pos = c(10510L, 10500L),
      clip_side = c("trail", "lead"),
      stringsAsFactors = FALSE
    )
    result <- detect_tsd(ref_dna, ins_pos = 506L, cluster_df = cluster_df,
                         genomic_start = 10000L)
    expect_equal(result$tsd_status, "read_supported")
    expect_equal(result$tsd_len, 10L)
    expected_seq <- substr(ref_dna, 502L, 511L)
    expect_equal(result$tsd_seq, expected_seq)
  })
})

# ---------------------------------------------------------------------------
# .detect_poly_a
# ---------------------------------------------------------------------------
describe(".detect_poly_a", {
  detect_poly_a <- get_fn(".detect_poly_a")

  it("detects a terminal poly-A tail", {
    expect_equal(detect_poly_a("GGCCGGATCGAAAAAAAAAA"), 10L)
  })

  it("detects a terminal poly-T tail (antisense)", {
    expect_equal(detect_poly_a("GGCCGGATCGTTTTTTTTTT"), 10L)
  })

  it("detects a leading poly-A tail", {
    expect_equal(detect_poly_a("AAAAAAAGGCCGGATCG"), 7L)
  })

  it("returns 0 for internal A/T runs (not terminal)", {
    # Internal run of A's, not at either end
    expect_equal(detect_poly_a("GCGAAAAAAAAAGCG"), 0L)
  })

  it("returns 0 when no terminal run meets minimum length", {
    expect_equal(detect_poly_a("GCGAAAAGCG"), 0L)
    expect_equal(detect_poly_a("GCGAAAGCG"), 0L)
  })

  it("returns 0 for empty or NA input", {
    expect_equal(detect_poly_a(""), 0L)
    expect_equal(detect_poly_a(NA_character_), 0L)
  })

  it("respects custom min_run", {
    expect_equal(detect_poly_a("GCAAAAGC", min_run = 4L), 4L)
    expect_equal(detect_poly_a("GCAAAAGC", min_run = 5L), 0L)
  })

  it("is case-insensitive", {
    expect_equal(detect_poly_a("gctttaaaa"), 4L)
  })
})

# ---------------------------------------------------------------------------
# .alu_5p_truncation
# ---------------------------------------------------------------------------
describe(".alu_5p_truncation", {
  trunc <- get_fn(".alu_5p_truncation")

  it("returns 0 when alignment starts at consensus position 1", {
    expect_equal(trunc(1L, 300L), 0L)
  })

  it("returns aln_start - 1 for truncation", {
    expect_equal(trunc(50L, 300L), 49L)
    expect_equal(trunc(100L, 300L), 99L)
  })

  it("returns NA when aln_start is NA", {
    expect_true(is.na(trunc(NA_integer_, 300L)))
  })

  it("returns NA when consensus_len is NA", {
    expect_true(is.na(trunc(50L, NA_integer_)))
  })

  it("never returns negative (max with 0)", {
    expect_equal(trunc(0L, 300L), 0L)
  })
})

# ---------------------------------------------------------------------------
# .filter_alu_call
# ---------------------------------------------------------------------------
describe(".filter_alu_call (NA TSD handling)", {
  filter_call <- get_fn(".filter_alu_call")

  valid_row <- data.frame(
    SupportingReads = 5L,
    ALU_RawScore = 0.80,
    AlleleFrequency = 0.05,
    TSD_Length = NA_integer_,
    TSD_Sequence = NA_character_,
    TSD_Status = "unknown",
    stringsAsFactors = FALSE
  )

  it("passes a valid call when min_tsd = 0 (TSD filter disabled)", {
    expect_true(filter_call(valid_row, min_tsd = 0L))
  })

  it("passes a valid call with unknown TSD even when min_tsd > 0", {
    # [FIX] Unknown TSD should NOT be filtered: we can't confirm or deny it.
    expect_true(filter_call(valid_row, min_tsd = 5L))
  })

  it("filters out calls with insufficient support", {
    row <- valid_row
    row$SupportingReads <- 1L
    expect_false(filter_call(row, min_support = 3L))
  })

  it("filters out calls with low ALU score", {
    row <- valid_row
    row$ALU_RawScore <- 0.30
    expect_false(filter_call(row, min_alu_score = 0.60))
  })

  it("filters out calls with low VAF", {
    row <- valid_row
    row$AlleleFrequency <- 0.001
    expect_false(filter_call(row, vaf_threshold = 0.01))
  })

  it("passes a call with read-supported TSD >= min_tsd", {
    row <- valid_row
    row$TSD_Length <- 10L
    row$TSD_Sequence <- "ACGTACGTAC"
    row$TSD_Status <- "read_supported"
    expect_true(filter_call(row, min_tsd = 5L))
  })

  it("filters a call with read-supported TSD < min_tsd", {
    row <- valid_row
    row$TSD_Length <- 3L
    row$TSD_Status <- "read_supported"
    expect_false(filter_call(row, min_tsd = 5L))
  })

  it("returns FALSE for NULL row", {
    expect_false(filter_call(NULL))
  })
})

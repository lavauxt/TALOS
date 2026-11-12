# Regression tests for ALU helper functions in engine_alu.R
#
# Covers:
#   - .alu_clip_breakpoint() including the chromosome-boundary fix (Finding 2)
#   - .detect_alu_tsd() including the read-evidence approach (Finding 1)
#   - .detect_poly_a() terminal-run detection
#   - .alu_5p_truncation() edge cases
#   - .filter_alu_call() NA TSD handling
#   - .extract_alu_candidates() boundary guard
#   - .match_clip_to_alu() pattern alignment positions (ITEM 2)
#   - .extract_alu_boundary_bases() non-ALU clip boundary recovery (ITEM 2)
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

  it("rejects a call with unknown TSD when an explicit minimum is requested", {
    # An unknown TSD cannot satisfy an explicit minimum-length requirement.
    expect_false(filter_call(valid_row, min_tsd = 5L))
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

# ---------------------------------------------------------------------------
# .extract_alu_boundary_bases (ITEM 2: pattern alignment positions)
# ---------------------------------------------------------------------------
describe(".extract_alu_boundary_bases (ITEM 2)", {
  bbf <- get_fn(".extract_alu_boundary_bases")

  # Trail clip = upstream_TSD + ALU_start.
  # The TSD (non-ALU bases) is the PREFIX of the clip (before pattern_start).
  it("recovers the prefix bases before pattern_start for a trail clip", {
    # Clip = upstream_TSD + ALU_start
    # pattern_start = 9 -> boundary = clip[1..8] = "ACGTACGT" (upstream TSD)
    clip <- "ACGTACGTGGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTG"
    expect_equal(bbf(clip, 9L, 50L, "trail"), "ACGTACGT")
  })

  # Lead clip = ALU_end + downstream_TSD.
  # The TSD (non-ALU bases) is the SUFFIX of the clip (after pattern_end).
  it("recovers the suffix bases after pattern_end for a lead clip", {
    # Clip = ALU_end + downstream_TSD
    # pattern_end = 40 -> boundary = clip[41..50] = "TTTTTTTTTT" (downstream TSD)
    clip <- paste0("GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTG", "TTTTTTTTTT")
    expect_equal(bbf(clip, 1L, 40L, "lead"), "TTTTTTTTTT")
  })

  it("returns empty string when pattern_start == 1 (trail clip has no TSD prefix)", {
    # Trail clip is pure ALU from position 1 -> no upstream TSD
    expect_equal(bbf("ACGTACGT", 1L, 8L, "trail"), "")
  })

  it("returns empty string when pattern_end == nchar(clip) (lead clip has no TSD suffix)", {
    # Lead clip is pure ALU through to the end -> no downstream TSD
    clip <- "ACGTACGT"
    expect_equal(bbf(clip, 1L, nchar(clip), "lead"), "")
  })

  it("returns NA for invalid clip_side", {
    expect_true(is.na(bbf("ACGT", 1L, 4L, "middle")))
  })

  it("returns NA for NA inputs", {
    expect_true(is.na(bbf(NA_character_, 1L, 4L, "lead")))
    expect_true(is.na(bbf("ACGT", NA_integer_, 4L, "lead")))
    expect_true(is.na(bbf("ACGT", 1L, NA_integer_, "trail")))
  })

  it("returns NA when pattern positions are inconsistent", {
    # pattern_start > pattern_end
    expect_true(is.na(bbf("ACGTACGT", 5L, 3L, "lead")))
    # pattern_end beyond clip length
    expect_true(is.na(bbf("ACGT", 1L, 99L, "trail")))
    # pattern_start < 1
    expect_true(is.na(bbf("ACGT", 0L, 4L, "lead")))
  })

  it("returns NA for empty clip", {
    expect_equal(bbf("", 1L, 0L, "lead"), NA_character_)
  })

  it("handles a single boundary base (trail clip with 1bp TSD prefix)", {
    # clip = "A" (1bp TSD) + ALU starting at position 2
    expect_equal(bbf("AGGCCGGG", 2L, 8L, "trail"), "A")
  })

  it("handles a single boundary base (lead clip with 1bp TSD suffix)", {
    # clip = ALU ending at position 7 + "T" (1bp downstream TSD)
    expect_equal(bbf("GGCCGGGT", 1L, 7L, "lead"), "T")
  })

  it("returns the full prefix when pattern_start is near the end (trail clip)", {
    # Pathological case: pattern alignment covers only the last 3 bases
    clip <- "ACGTACGTACGGG"
    # pattern_start = 11, pattern_end = 13
    expect_equal(bbf(clip, 11L, 13L, "trail"), "ACGTACGTAC")
  })

  it("returns the full suffix when pattern_end is near the start (lead clip)", {
    # Pathological case: pattern alignment covers only the first 3 bases
    clip <- "GGGACGTACGTAC"
    # pattern_start = 1, pattern_end = 3
    expect_equal(bbf(clip, 1L, 3L, "lead"), "ACGTACGTAC")
  })
})

# ---------------------------------------------------------------------------
# .match_clip_to_alu pattern positions (ITEM 2)
# ---------------------------------------------------------------------------
# This test exercises the full .match_clip_to_alu() with Biostrings + pwalign.
# It is skipped when those packages are not installed, but when they are, it
# confirms that pattern_start and pattern_end are populated and consistent.
describe(".match_clip_to_alu (ITEM 2: pattern positions stored)", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("pwalign")

  match_fn <- get_fn(".match_clip_to_alu")

  # A small synthetic ALU consensus (~60bp fragment) to avoid depending on
  # the bundled FASTA file. ALU elements are A-rich; this fragment is
  # inspired by the AluSx 5' monomer.
  alu_consensus <- Biostrings::DNAStringSet(list(
    AluSx_test = "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATC"
  ))

  it("returns NULL for a clip shorter than min_clip_len", {
    expect_null(match_fn("ACGTACGT", alu_consensus, min_clip_len = 25L))
  })

  it("returns NULL for NA clip", {
    expect_null(match_fn(NA_character_, alu_consensus))
  })

  it("returns NULL when no ALU matches above min_score", {
    # A random non-ALU sequence, 50bp
    non_alu <- "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    expect_null(match_fn(non_alu, alu_consensus, min_score = 0.95))
  })

  it("stores pattern_start and pattern_end alongside aln_start/aln_end", {
    # Exact ALU consensus subsequence -> should match perfectly
    clip <- "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGC"
    result <- match_fn(clip, alu_consensus, min_score = 0.5)
    expect_false(is.null(result))
    expect_true(is.integer(result$pattern_start))
    expect_true(is.integer(result$pattern_end))
    expect_false(is.na(result$pattern_start))
    expect_false(is.na(result$pattern_end))
    # Pattern positions must be within the clip
    expect_gte(result$pattern_start, 1L)
    expect_lte(result$pattern_end, nchar(clip))
    expect_lte(result$pattern_start, result$pattern_end)
  })

  it("normalises by aligned ALU width and maps antisense positions to the forward clip", {
    alu <- as.character(alu_consensus[[1L]])
    rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(alu)))
    clip <- paste0("TTTTTTTT", substr(rc, 1L, 50L), "AAAAAAAA")
    result <- match_fn(clip, alu_consensus, min_score = 0.5)
    expect_false(is.null(result))
    expect_equal(result$strand, "-")
    expect_gte(result$score, 0.9)
    expect_equal(result$pattern_start, 9L)
    expect_equal(result$pattern_end, 58L)
    expect_equal(result$aligned_bases, 50L)
    expect_lt(result$clip_fraction, 1)
  })

  it("pattern_start > 1 when the clip has a non-ALU prefix", {
    # 8 bases of TSD prefix + 50 bases of ALU. This is the structure of a
    # TRAIL clip from an upstream-junction-spanning read; the prefix is the
    # upstream TSD carried on the read.
    tsd_prefix <- "ACGTACGT"
    alu_part   <- "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGC"
    clip <- paste0(tsd_prefix, alu_part)
    result <- match_fn(clip, alu_consensus, min_score = 0.5)
    expect_false(is.null(result))
    # The alignment should start somewhere after position 1, since the
    # first 8 bases are not ALU-derived.
    expect_gt(result$pattern_start, 1L)
  })

  it("pattern_end < nchar(clip) when the clip has a non-ALU suffix", {
    # 50 bases of ALU + 10 bases of TSD suffix. This is the structure of a
    # LEAD clip from a downstream-junction-spanning read; the suffix is the
    # downstream TSD carried on the read.
    alu_part   <- "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGC"
    tsd_suffix <- "TTTTTTTTTT"
    clip <- paste0(alu_part, tsd_suffix)
    result <- match_fn(clip, alu_consensus, min_score = 0.5)
    expect_false(is.null(result))
    # The alignment should end before the clip's last base, since the
    # trailing 10 bases are not ALU-derived.
    expect_lt(result$pattern_end, nchar(clip))
  })
})

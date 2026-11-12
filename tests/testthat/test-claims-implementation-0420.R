# Regression tests for the implemented v0.4.20 claims (TKD scope decision excluded)
#
# Covers:
#   - Adapter screening: .screen_adapter_seqs / .screen_adapter_reads
#   - Base-quality filtering: .mean_base_quality / .filter_reads_by_base_quality /
#     .bam_has_base_quality and per-read BQ carried into candidate rows
#   - Length-scaled breakpoint clustering: .cluster_breakpoints with
#     lengths / length_scale / max_scaled_tolerance (anti-chaining kept)
#   - Exact QNAME union support merging: .merge_ptd_intervals union semantics
#     and the explicit sum fallback

library(testthat)

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
# Adapter screening
# ---------------------------------------------------------------------------
describe("adapter screening", {
  screen_seqs <- get_fn(".screen_adapter_seqs")
  adapter <- "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" # Illumina TruSeq Read 1

  it("flags reads carrying an adapter prefix of at least min_overlap bases", {
    seqs <- c(
      "AAAAGATCGGAAGAGC",   # contains 12 bases of the adapter prefix
      "TTTTTTTTTTTTTTTTTT", # clean genomic sequence
      NA_character_
    )
    hit <- screen_seqs(seqs, adapter, min_overlap = 10L)
    expect_true(is.logical(hit))
    expect_equal(length(hit), 3L)
    expect_true(hit[[1L]])
    expect_false(hit[[2L]])
    expect_false(hit[[3L]])
  })

  it("respects the min_overlap threshold", {
    seqs <- c("AAAAGATCGGAAGAGC") # longest contiguous adapter prefix: 12
    expect_true(screen_seqs(seqs, adapter, min_overlap = 10L))
    expect_true(screen_seqs(seqs, adapter, min_overlap = 12L))
    expect_false(screen_seqs(seqs, adapter, min_overlap = 13L))
  })

  it("matches adapter prefixes case-insensitively and ignores empty adapters", {
    expect_true(screen_seqs("aaaagatcggaagagc", adapter, min_overlap = 10L))
    expect_equal(sum(screen_seqs(c("ACGTACGTAC", "TTTTTTTTTT"), c("", NA), 10L)), 0L)
  })

  it("reports region-wide counts and unique affected QNAMEs", {
    skip_if_not_installed("GenomicAlignments")
    skip_if_not_installed("S4Vectors")
    skip_if_not_installed("Biostrings")
    screen_reads <- get_fn(".screen_adapter_reads")
    reads <- GenomicAlignments::GAlignments(
      seqnames = S4Vectors::Rle("chr11", 3L),
      pos = as.integer(c(100, 200, 300)),
      cigar = c("16M", "18M", "16M"),
      strand = S4Vectors::Rle(strand("+"), 3L),
      names = c("r1", "r2", "r3")
    )
    S4Vectors::mcols(reads)$seq <- Biostrings::DNAStringSet(c(
      "AAAAGATCGGAAGAGC",
      "TTTTTTTTTTTTTTTTTT",
      "GGGGGGGGGGGGGGGG"
    ))
    out <- screen_reads(reads, c(adapter, "CTGTCTCTTATACACATCT"), 10L)
    expect_equal(out$screened, 3L)
    expect_equal(out$positive, 1L)
    expect_equal(out$qnames, "r1")
    # No adapters configured -> screening disabled, zero hits
    out_off <- screen_reads(reads, character(0), 10L)
    expect_equal(out_off$positive, 0L)
    expect_equal(length(out_off$qnames), 0L)
  })
})

# ---------------------------------------------------------------------------
# Base-quality filtering
# ---------------------------------------------------------------------------
describe("base-quality filtering", {
  mean_bq <- get_fn(".mean_base_quality")
  has_bq <- get_fn(".bam_has_base_quality")
  bq_filter <- get_fn(".filter_reads_by_base_quality")

  it("computes mean Phred+33 quality per read", {
    # 'I' = 73 - 33 = Phred 40; '!' = 0; '#' = 2
    got <- mean_bq(c("IIII", "!!!!", "#!!!", "", NA_character_))
    expect_equal(got[[1L]], 40)
    expect_equal(got[[2L]], 0)
    expect_equal(round(got[[3L]], 3), round((2 + 0 + 0 + 0) / 4, 3))
    expect_true(is.na(got[[4L]]))
    expect_true(is.na(got[[5L]]))
  })

  it("detects whether loaded reads carry a QUAL column", {
    skip_if_not_installed("GenomicAlignments")
    skip_if_not_installed("S4Vectors")
    skip_if_not_installed("Biostrings")
    mk <- function(quals) {
      r <- GenomicAlignments::GAlignments(
        seqnames = S4Vectors::Rle("chr11", 2L),
        pos = c(100L, 200L), cigar = c("10M", "10M"),
        strand = S4Vectors::Rle(strand("+"), 2L), names = c("a", "b")
      )
      if (!is.null(quals)) S4Vectors::mcols(r)$qual <- quals
      r
    }
    expect_false(has_bq(mk(NULL)))
    expect_false(has_bq(mk(Biostrings::PhredQuality(c("", "")))))
    expect_true(has_bq(mk(Biostrings::PhredQuality(c("IIII", "!!!!")))))
  })

  it("drops reads below the mean-quality threshold and keeps unevaluated ones", {
    skip_if_not_installed("GenomicAlignments")
    skip_if_not_installed("S4Vectors")
    skip_if_not_installed("Biostrings")
    reads <- GenomicAlignments::GAlignments(
      seqnames = S4Vectors::Rle("chr11", 3L),
      pos = c(100L, 200L, 300L), cigar = c("10M", "10M", "10M"),
      strand = S4Vectors::Rle(strand("+"), 3L), names = c("good", "bad", "noqual")
    )
    S4Vectors::mcols(reads)$qual <- Biostrings::PhredQuality(c("IIII", "!!!!", ""))
    kept <- bq_filter(reads, 20)
    expect_equal(length(kept), 2L)
    expect_equal(names(kept), c("good", "noqual"))
  })
})

# ---------------------------------------------------------------------------
# Per-read base quality flows into candidate rows
# ---------------------------------------------------------------------------
describe("candidate rows carry read_mean_bq", {
  new_candidate <- get_fn(".new_candidate")
  to_df <- get_fn(".candidates_to_df")

  it("round-trips the base quality through the candidate record", {
    rec <- new_candidate("r1", 10L, 30L, "anchor_ext", 60L, FALSE, "10S90M",
                         "ACGTACGTACGT", itd_seq = "ACGTACGTACGT",
                         read_start = 100L, read_end = 200L, read_mean_bq = 35.5)
    df <- to_df(list(rec, new_candidate("r2", 12L, 30L, "anchor_ext", 60L, FALSE,
                                        "10S90M", "ACGTACGTACGT")))
    expect_true("read_mean_bq" %in% names(df))
    expect_equal(df$read_mean_bq[[1L]], 35.5)
    expect_true(is.na(df$read_mean_bq[[2L]]))
  })
})

# ---------------------------------------------------------------------------
# Length-scaled breakpoint clustering
# ---------------------------------------------------------------------------
describe("length-scaled breakpoint clustering", {
  cluster <- get_fn(".cluster_breakpoints")

  it("keeps the legacy scalar behaviour by default and with length_scale = 0", {
    bps <- c(100, 108, 116, NA_real_)
    expect_equal(cluster(bps, 10), cluster(bps, 10, lengths = c(300, 300, 300, 300), length_scale = 0))
    got <- cluster(c(100, 108, 116), 10, lengths = c(300, 300, 300), length_scale = 0)
    expect_length(got, 2L)
    expect_equal(got[[1L]], c(100, 108))
    expect_equal(got[[2L]], 116)
  })

  it("widens the tolerance for long duplications", {
    # 0.1 * 300 = 30 > cluster_tolerance 10: a 15 bp separation now clusters.
    got <- cluster(c(100, 115), 10, lengths = c(300, 300), length_scale = 0.1,
                   max_scaled_tolerance = Inf)
    expect_length(got, 1L)
    expect_equal(got[[1L]], c(100, 115))
    # Short duplications keep the base tolerance: 0.1 * 30 = 3 < 10.
    got_short <- cluster(c(100, 115), 10, lengths = c(30, 30), length_scale = 0.1,
                         max_scaled_tolerance = Inf)
    expect_length(got_short, 2L)
  })

  it("caps the scaled tolerance at max_scaled_tolerance", {
    # 0.1 * 1000 = 100 -> capped at 50: 55 bp apart stays separate.
    got <- cluster(c(100, 155, 156), 10, lengths = c(1000, 1000, 1000),
                   length_scale = 0.1, max_scaled_tolerance = 50)
    expect_length(got, 2L)
    expect_equal(got[[1L]], 100)
    expect_equal(got[[2L]], c(155, 156))
  })

  it("still bounds the total cluster width while scaling", {
    # Anti-chaining: distances are 8 and 8, but the width from the seed to
    # the last member (16) must not chain beyond the effective tolerance (30).
    got <- cluster(c(100, 108, 116, 124), 10, lengths = rep(300, 4),
                   length_scale = 0.1, max_scaled_tolerance = 30)
    expect_length(got, 1L) # width 24 <= 30
    got2 <- cluster(c(100, 108, 116, 124, 132), 10, lengths = rep(300, 5),
                    length_scale = 0.1, max_scaled_tolerance = 30)
    expect_length(got2, 2L) # width 32 > 30 -> the tail splits off
  })

  it("handles NA lengths and mismatched length vectors gracefully", {
    got <- cluster(c(100, 115), 10, lengths = c(NA_real_, NA_real_),
                   length_scale = 0.1, max_scaled_tolerance = 50)
    expect_length(got, 2L) # no usable lengths -> base tolerance 10
    got2 <- cluster(c(100, 115), 10, lengths = c(300), length_scale = 0.1)
    expect_length(got2, 2L) # length vector ignored when its size mismatches
  })

  it("rejects invalid scaling parameters", {
    expect_error(cluster(c(100, 115), 10, length_scale = 1.5), "between 0 and 1")
    expect_error(cluster(c(100, 115), 10, length_scale = -0.1), "between 0 and 1")
    expect_error(cluster(c(100, 115), 10, max_scaled_tolerance = 5), "cluster_tolerance")
  })
})

# ---------------------------------------------------------------------------
# Exact QNAME-union support merging
# ---------------------------------------------------------------------------
describe("exact QNAME union support merging", {
  merge_fn <- get_fn(".merge_ptd_intervals")
  qname_str <- get_fn(".qname_set_to_string")
  quiet_log <- function(...) invisible(NULL)

  it("serializes QNAME sets deterministically and truncates oversized sets", {
    expect_equal(qname_str(c("r3", "r1", "r2", "r1")), "r1;r2;r3")
    expect_true(is.na(qname_str(character(0))))
    expect_true(is.na(qname_str(c(NA, NA))))
    big <- paste0("read", seq_len(5000L))
    s <- qname_str(big, max_chars = 1000L)
    expect_true(endsWith(s, "...TRUNCATED"))
  })

  mk_df <- function(with_sets = TRUE) {
    df <- data.frame(
      GenomicPosition = c(100, 121),
      Length = c(20, 20),
      SupportingReads = c(5, 5),
      WildtypeReads = c(10, 10),
      DepthAtBreakpoint = c(15, 15),
      AlleleFrequency = c(5 / 15, 5 / 15),
      IndependentSupportQNAMEs = c(5, 5),
      AdjustedSupportingReads = c(5, 5),
      SizeCorrectionFactor = c(1, 1),
      CorrectedDepthAtBreakpoint = c(15, 15),
      CorrectedAlleleFrequency = c(5 / 15, 5 / 15),
      SpanMinCoverage = c(NA_real_, NA_real_),
      SpanDepthFoldChange = c(NA_real_, NA_real_),
      stringsAsFactors = FALSE
    )
    if (with_sets) {
      df$SupportQNAMEs <- c("r1;r2;r3", "r2;r4;r5") # r2 shared between both fragments
      df$WildtypeQNAMEs <- c("w1;w2;w3", "w1;w2;w3")
    }
    df
  }

  it("merges support as the exact QNAME union, not the component-row sum", {
    out <- merge_fn(
      mk_df(TRUE), min_ptd_length = 10, merge_gap = 10, max_plausible_fragments = 4,
      thresholds = list(), cov = NULL, chrom = "chr11", coverage_window = 200,
      log_msg = quiet_log
    )
    expect_equal(nrow(out), 1L)
    # Union of r1,r2,r3 and r2,r4,r5 -> 5 (a sum would give 10)
    expect_equal(out$SupportingReads[[1L]], 5L)
    # Union of identical WT sets -> 3 (a sum would give 20)
    expect_equal(out$WildtypeReads[[1L]], 3L)
    expect_equal(out$DepthAtBreakpoint[[1L]], 8L)
    expect_equal(out$AlleleFrequency[[1L]], 5 / 8)
    expect_equal(out$RawAlleleFrequency[[1L]], 5 / 8)
    expect_equal(out$IndependentSupportQNAMEs[[1L]], 5L)
    expect_equal(out$SupportQNAMEs[[1L]], "r1;r2;r3;r4;r5")
    expect_equal(out$WildtypeQNAMEs[[1L]], "w1;w2;w3")
    expect_equal(out$MergedSupportSemantics[[1L]], "exact_qname_union_across_component_rows")
    # Corrected diagnostic recomputed consistently from the union count
    expect_equal(out$AdjustedSupportingReads[[1L]], 5L)
    expect_equal(out$CorrectedDepthAtBreakpoint[[1L]], 8L)
    expect_equal(out$MergedFromN[[1L]], 2L)
  })

  it("falls back to the documented sum when QNAME sets are absent", {
    out <- merge_fn(
      mk_df(FALSE), min_ptd_length = 10, merge_gap = 10, max_plausible_fragments = 4,
      thresholds = list(), cov = NULL, chrom = "chr11", coverage_window = 200,
      log_msg = quiet_log
    )
    expect_equal(out$SupportingReads[[1L]], 10L)
    expect_equal(out$WildtypeReads[[1L]], 20L)
    expect_equal(out$MergedSupportSemantics[[1L]],
                 "sum_of_component_rows; shared_QNAMEs_may_be_double_counted")
  })

  it("falls back to the sum when a QNAME set was truncated", {
    df <- mk_df(TRUE)
    df$SupportQNAMEs <- c("r1;r2;r3;...TRUNCATED", "r2;r4;r5")
    out <- merge_fn(
      df, min_ptd_length = 10, merge_gap = 10, max_plausible_fragments = 4,
      thresholds = list(), cov = NULL, chrom = "chr11", coverage_window = 200,
      log_msg = quiet_log
    )
    expect_equal(out$SupportingReads[[1L]], 10L)
    expect_equal(out$MergedSupportSemantics[[1L]],
                 "sum_of_component_rows; shared_QNAMEs_may_be_double_counted")
  })
})

# ---------------------------------------------------------------------------
# detect_itd() parameter validation for the new claims
# ---------------------------------------------------------------------------
describe("detect_itd validates the new claim parameters", {
  it("rejects invalid min_base_quality / adapter / scaling parameters", {
    expect_error(
      detect_itd("x.bam", list(gene = "FLT3"), min_base_quality = -1),
      "min_base_quality"
    )
    expect_error(
      detect_itd("x.bam", list(gene = "FLT3"), adapter_min_overlap = 0),
      "adapter_min_overlap"
    )
    expect_error(
      detect_itd("x.bam", list(gene = "FLT3"), cluster_length_scale = 2),
      "cluster_length_scale"
    )
    expect_error(
      detect_itd("x.bam", list(gene = "FLT3"), max_cluster_tolerance = 5),
      "max_cluster_tolerance"
    )
  })
})

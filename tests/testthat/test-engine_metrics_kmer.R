# Regression tests for the k-mer ITD coverage feature and the display fold
# (engine_metrics.R)
#
# Covers .kmer_itd_coverage() behind the ITDCoverageRLE_Kmer /
# ITDReadCoverage_Kmer / ITDKmerReads output columns, and
# .fold_display_coverage() behind the ITDCoverageRLE display vector:
#   - per-read union semantics (one read adds at most +1 per base, the
#     v7 fix for the k-mer-match double counting)
#   - depth stacking across reads, bounded by the read count
#   - partial matches, no-match reads, NA / short / empty inputs
#   - case normalisation between assembled ITD sequence and BAM reads
#   - reverse-complement reads via the RC k-mer index
#   - repeated k-mers being excluded as position-ambiguous
#   - junction reads marking both uniquely mappable arms without assigning
#     repeated motifs to every copy
#   - the specificity-ordered fold (direct > PE > k-mer, binary direct
#     declared by the caller as least specific - v8)
#
# The DNA literals below are pre-verified to contain no repeated 11-mers
# and no 11-mer shared with their own reverse complement (and the tandem
# components share no 11-mer with each other), so every expectation below
# is structural, not statistical: no RNG, no flaky seeds.

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

describe(".kmer_itd_coverage", {
  kic <- get_fn(".kmer_itd_coverage")
  k <- 11L

  # 120bp, no repeated 11-mers, no RC-shared 11-mers (verified upstream).
  itd <- paste0(
    "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTG",
    "TTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCG"
  )
  # 120bp sharing no 11-mer with `itd` on either strand (verified upstream).
  nomatch <- paste0(
    "GCAGATCAAGCAGGAGGCGGAATGTAAACAGAAGGTATGCTTAGGTGGATAGGGAGTGAGC",
    "AACAAACGGATCGTTTCTCCCATGCCAAGTTGGCACAGGGAACTACCTGCGGCGGTTTG"
  )
  # Tandem-duplication components: unit repeated twice + unrelated tail,
  # mutually free of shared 11-mers (verified upstream).
  unit <- "ATCCGTTCCTAATAAGGAATGGTGATTCCC"
  tail <- paste0(
    "TGTCATACCAATCTACCCCCTGTTATGCGCGTTTGTCGTTAGACCAATGTCAGCGCAGCG"
  )

  it("returns an empty result for a sequence shorter than k", {
    res <- kic(substr(itd, 1L, k - 1L), c(itd, itd))
    expect_equal(length(res$cov_vec), 0L)
    expect_true(is.na(res$rle_str))
    expect_equal(res$percent, 0)
    expect_equal(res$n_reads, 0L)
  })

  it("returns an empty result for NULL / empty read sequences", {
    res <- kic(itd, character(0))
    expect_equal(length(res$cov_vec), 0L)
    expect_equal(res$n_reads, 0L)
  })

  it("rejects non-integer or non-numeric k-mer settings", {
    expect_error(kic(itd, itd, k = "eleven"), "positive integers")
    expect_error(kic(itd, itd, min_offset_hits = 1.5), "positive integers")
  })

  it("counts a full-length matching read as +1 per base, not +1 per k-mer", {
    # THE v7 regression: one read whose ~110 k-mers all match used to pile
    # ~11 counts onto every base (values > 100 for a single read).
    res <- kic(itd, itd)
    expect_equal(res$n_reads, 1L)
    expect_true(all(res$cov_vec == 1L))
    expect_equal(res$percent, 100)
    expect_equal(res$rle_str, "1:120")
  })

  it("stacks depth across identical reads", {
    res <- kic(itd, c(itd, itd, itd))
    expect_equal(res$n_reads, 3L)
    expect_equal(res$rle_str, "3:120")
  })

  it("never exceeds the number of reads at any base", {
    reads <- rep(substr(itd, 1L, 60L), 12L)
    res <- kic(itd, reads)
    expect_equal(res$n_reads, 12L)
    expect_true(max(res$cov_vec) <= 12L)
  })

  it("covers exactly the read's extent for an internal exact match", {
    # read = itd[51..110]: its 50 k-mers start at ITD positions 51..100, so
    # the covered bases are exactly 51..110.
    read <- substr(itd, 51L, 110L)
    res <- kic(itd, read)
    expect_equal(res$cov_vec, c(integer(50L), rep(1L, 60L), integer(10L)))
    expect_equal(res$percent, 50)
    expect_equal(res$rle_str, "0:50,1:60,0:10")
  })

  it("ignores a read with no k-mer match", {
    res <- kic(itd, nomatch)
    expect_equal(res$n_reads, 0L)
    expect_equal(res$rle_str, "0:120")
    expect_equal(res$percent, 0)
  })

  it("returns one match flag per read end for fragment-level aggregation", {
    res <- kic(itd, c(itd, nomatch, substr(itd, 51L, 110L)))
    expect_identical(res$matched_reads, c(TRUE, FALSE, TRUE))
    expect_equal(sum(res$matched_reads), res$n_reads)
  })

  it("skips NA and too-short reads without failing", {
    res <- kic(itd, c(NA_character_, substr(itd, 1L, k - 1L), itd))
    expect_equal(res$n_reads, 1L)
    expect_equal(res$rle_str, "1:120")
  })

  it("matches reads against a lower-case ITD sequence", {
    # BAM SEQ is uppercase; assembled ITD sequences may not be.
    res <- kic(tolower(itd), itd)
    expect_equal(res$n_reads, 1L)
    expect_equal(res$rle_str, "1:120")
  })

  it("matches reverse-complement reads through the RC index", {
    skip_if_not(requireNamespace("Biostrings", quietly = TRUE))
    rc_read <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(
        substr(itd, 21L, 100L)
      ))
    )
    res <- kic(itd, rc_read)
    expect_equal(res$n_reads, 1L)
    # The RC read marks the same 80bp extent as the forward read would.
    expect_equal(res$cov_vec, c(integer(20L), rep(1L, 80L), integer(20L)))
  })

  it("returns an explicit ambiguous-coverage vector", {
    res <- kic(paste0(unit, unit, tail), unit, k = 7L, min_offset_hits = 2L)
    expect_length(res$ambiguous_cov_vec, nchar(paste0(unit, unit, tail)))
    expect_true(all(res$ambiguous_cov_vec >= 0L))
  })

  it("does not assign a repeated k-mer to every copy of a tandem repeat", {
    # The 30bp unit occurs twice in itd2. Its 11-mers cannot localise a read
    # to either copy, so they must not inflate positional coverage of both.
    itd2 <- paste0(unit, unit, tail)
    res <- kic(itd2, unit)
    expect_equal(res$n_reads, 0L)
    expect_equal(res$percent, 0)
    expect_true(all(res$cov_vec == 0L))
  })

  it("a junction read covers the ITD end and the ITD start in one pass", {
    # Junction read = last 40bp of the ITD + first 40bp again (what a read
    # spanning the tandem junction actually looks like). The 11-mers
    # straddling the junction are novel and must not match anything.
    jread <- paste0(substr(itd, 81L, 120L), substr(itd, 1L, 40L))
    res <- kic(itd, jread)
    expect_equal(res$n_reads, 1L)
    cov <- res$cov_vec
    expect_true(all(cov[1:40] == 1L))
    expect_true(all(cov[41:80] == 0L))
    expect_true(all(cov[81:120] == 1L))
    expect_equal(max(cov), 1L)
  })

  it("emits RLE strings in the existing ITDCoverageRLE format", {
    res <- kic(itd, c(itd, substr(itd, 51L, 110L)))
    expect_equal(res$rle_str, "1:50,2:60,1:10")
    parts <- strsplit(res$rle_str, ",", fixed = TRUE)[[1]]
    vals <- as.integer(sub(":.*$", "", parts))
    lens <- as.integer(sub("^.*:", "", parts))
    expect_equal(sum(lens), 120L)
    expect_equal(sum(lens[vals > 0]), 120L)  # full-length read covers all
  })
})

describe(".fold_display_coverage", {
  fdc <- get_fn(".fold_display_coverage")

  it("returns integer(0) when every source is empty", {
    expect_equal(length(fdc(integer(0), integer(0), integer(0))), 0L)
    expect_equal(length(fdc(NULL, NULL, NULL)), 0L)
  })

  it("returns the single source unchanged (identity fold)", {
    d <- c(1L, 2L, 0L, 0L)
    expect_equal(fdc(d, integer(0), integer(0)), d)
    p <- c(0L, 5L, 6L, 0L)
    expect_equal(fdc(integer(0), p, integer(0)), p)
    km <- c(0L, 0L, 0L, 9L)
    expect_equal(fdc(integer(0), integer(0), km), km)
  })

  it("fills direct-only zeros with PE counts and keeps direct counts", {
    d <- c(2L, 0L, 0L, 0L)
    p <- c(0L, 7L, 8L, 0L)
    expect_equal(fdc(d, p, integer(0)), c(2L, 7L, 8L, 0L))
  })

  it("fills remaining zeros with k-mer counts (specificity order)", {
    d <- c(2L, 0L, 0L, 0L)
    p <- c(0L, 7L, 0L, 0L)
    km <- c(0L, 0L, 0L, 9L)
    expect_equal(fdc(d, p, km), c(2L, 7L, 0L, 9L))
  })

  it("never lets a less specific source overwrite direct evidence", {
    d <- c(1L, 1L, 0L, 0L)
    p <- c(99L, 99L, 99L, 0L)
    expect_equal(fdc(d, p, integer(0)), c(1L, 1L, 99L, 0L))
  })

  it("treats a binary (0/1) direct vector as the least specific layer", {
    # The no-pwalign fallback emits covered-flags, not counts; real PE
    # counts must not be flattened to 1. The fallback is DECLARED by the
    # caller via direct_is_binary = TRUE: the values alone cannot
    # distinguish flags from a low-depth count vector (every base covered
    # by exactly one read), which is why v9 replaced value inference.
    d <- c(1L, 1L, 0L, 0L)
    p <- c(0L, 7L, 8L, 0L)
    expect_equal(fdc(d, p, integer(0), direct_is_binary = TRUE), c(1L, 7L, 8L, 0L))
  })

  it("a binary direct vector still fills bases no other source covers", {
    d <- c(1L, 0L, 1L, 0L)
    p <- c(0L, 0L, 0L, 0L)
    expect_equal(fdc(d, p, integer(0), direct_is_binary = TRUE), c(1L, 0L, 1L, 0L))
    km <- c(0L, 0L, 0L, 4L)
    expect_equal(fdc(d, p, km, direct_is_binary = TRUE), c(1L, 0L, 1L, 4L))
  })

  it("ignores sources whose length disagrees with the longest source", {
    d <- c(2L, 0L, 0L, 0L)
    p_short <- c(5L, 5L)
    km <- c(0L, 0L, 0L, 9L)
    # PE shorter -> dropped; direct + k-mer fold.
    expect_equal(fdc(d, p_short, km), c(2L, 0L, 0L, 9L))
  })

  it("produces the filled track for the validation-sample pattern", {
    # Direct covers 1-50 and 200-300; PE covers 51-199 (depth 40); k-mer
    # covers nothing new -> display has no internal zero gap.
    d <- c(rep(3L, 50L), integer(149L), rep(4L, 101L))
    p <- c(integer(50L), rep(40L, 149L), integer(101L))
    disp <- fdc(d, p, integer(0))
    expect_equal(length(disp), 300L)
    expect_equal(sum(disp == 0L), 0L)  # internal gap filled
    expect_equal(disp[1:50], rep(3L, 50L))      # direct preserved
    expect_equal(disp[51:199], rep(40L, 149L))  # PE fill
    expect_equal(disp[200:300], rep(4L, 101L))  # direct preserved
  })
})


describe("ITD coverage and partial-sequence diagnostics", {
  get_uncovered <- get_fn(".coverage_uncovered_summary")
  local_match <- get_fn(".local_sequence_match")

  it("returns 1-based inclusive intervals for uncovered bases", {
    got <- get_uncovered(c(1L, 0L, 0L, 2L, 0L, 0L, 0L, 1L))
    expect_equal(got$bases, 5L)
    expect_equal(got$intervals, "2-3|5-7")
  })

  it("returns an empty interval string for fully covered sequences", {
    got <- get_uncovered(c(1L, 2L, 1L))
    expect_equal(got$bases, 0L)
    expect_identical(got$intervals, "")
  })

  it("aligns a partial observed fragment at an internal reference offset", {
    skip_if_not(requireNamespace("pwalign", quietly = TRUE))
    skip_if_not(requireNamespace("Biostrings", quietly = TRUE))
    subject <- paste0(unit, unit, tail)
    query <- substr(subject, 66L, 105L)
    got <- local_match(query, subject)
    expect_equal(got$pid, 100)
    expect_equal(got$query_coverage, 1)
    expect_equal(got$aligned_width, nchar(query))
  })
})

describe(".multi_kmer_itd_coverage", {
  mkic <- get_fn(".multi_kmer_itd_coverage")
  itd <- paste0(
    "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTG",
    "TTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCG"
  )

  it("reports every requested k at every nucleotide offset", {
    res <- mkic(itd, c(itd, paste0(substr(itd, 11, 100), "ACGTACGT")))
    expect_identical(names(res$per_k), c("7", "9", "11", "13", "15"))
    expect_length(res$cov_vec, nchar(itd))
    expect_true(all(vapply(res$per_k, function(z) length(z$cov_vec) == nchar(itd), logical(1))))
    expect_true(max(res$cov_vec) <= 2L)
    expect_equal(res$percent, 100 * sum(res$cov_vec > 0L) / nchar(itd))
  })

  it("recovers a read whose start is shifted relative to a codon/frame origin", {
    shifted_read <- substr(itd, 3L, 92L)
    res <- mkic(itd, shifted_read, k_values = c(7L, 9L, 11L),
                min_offset_hits = c(3L, 2L, 2L))
    expect_gt(res$percent, 0)
    expect_true(any(res$cov_vec > 0L))
  })

  it("does not assign a repeated tandem-unit k-mer to either copy", {
    unit <- "ATCCGTTCCTAATAAGGAATGGTGATTCCC"
    res <- mkic(paste0(unit, unit), unit, k_values = c(7L, 9L, 11L),
                min_offset_hits = c(3L, 2L, 2L))
    expect_equal(res$percent, 0)
    expect_equal(res$n_reads, 0L)
    expect_true(all(vapply(res$per_k, function(z) z$percent == 0, logical(1))))
  })

  it("requires unique k values so thresholds cannot be silently misassigned", {
    expect_error(mkic(itd, itd, k_values = c(7L, 7L)), "must not contain duplicates")
  })
})

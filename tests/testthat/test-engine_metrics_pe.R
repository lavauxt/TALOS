# Regression tests for the PE-spanning ITD coverage feature (engine_metrics.R)
#
# Covers .pe_spanning_coverage(), the pure helper behind the new
# ITDCoverageRLE_PE / ITDReadCoverage_PE / ITDPE_SpanningPairs (diagnostic only).
# The helper's percent_combined remains a legacy geometric-union diagnostic;
# the pipeline no longer uses it for allele-specific combined coverage.
# output columns:
#   - mate-inside-segment contribution over the mate's own extent
#   - exclusion of mates not fully inside the duplicated segment
#   - FF / RR orientation skip (FR and RF kept)
#   - support-qname and validity gating
#   - per-fragment union when both mates are inside (overlap counted once)
#   - clamping to seq_len when the ITD sequence is shorter than the span
#   - combined (union) percentage with the direct-coverage vector
#   - depth stacking across pairs
#   - degenerate inputs and the direct_cov length contract
#
# Plain-vector interface only: no BAM / GAlignmentPairs machinery needed.

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

describe(".pe_spanning_coverage", {
  psc <- get_fn(".pe_spanning_coverage")

  # Duplicated reference segment under test: [1000, 1299]
  # (breakpoint = 1000, span_len = 300, seq_len = 300 unless stated).
  # Defaults: mate1 [900, 1049] crosses the breakpoint (not fully inside),
  # mate2 [1100, 1249] lies entirely inside -> covers ITD positions 101-250.
  call <- function(qname = "r1", s1 = 900L, e1 = 1049L,
                   s2 = 1100L, e2 = 1249L,
                   orientation = "FR", valid = TRUE,
                   support_qnames = "r1",
                   breakpoint = 1000L, span_len = 300L, seq_len = 300L,
                   direct_cov = integer(0)) {
    psc(qname = qname, s1 = s1, e1 = e1, s2 = s2, e2 = e2,
        orientation = orientation, valid = valid,
        support_qnames = support_qnames,
        breakpoint = breakpoint, span_len = span_len, seq_len = seq_len,
        direct_cov = direct_cov)
  }

  it("returns an empty result for NULL / empty qname", {
    res <- call(qname = character(0))
    expect_true(is.na(res$percent))
    expect_true(is.na(res$rle_str))
    expect_true(is.na(res$percent_combined))
    expect_equal(res$n_spanning_pairs, 0L)
  })

  it("returns an empty result for NA breakpoint or non-positive span/seq_len", {
    expect_true(is.na(call(breakpoint = NA_integer_)$percent))
    expect_true(is.na(call(span_len = 0L)$percent))
    expect_true(is.na(call(span_len = NA_integer_)$percent))
    expect_true(is.na(call(seq_len = 0L)$percent))
  })

  it("errors when direct_cov does not have length seq_len", {
    expect_error(call(direct_cov = integer(5)), "direct_cov")
  })

  it("counts a mate that lies entirely inside the segment over its own extent", {
    res <- call()  # mate2 [1100, 1249] inside [1000, 1299] -> positions 101-250
    expect_equal(res$rle_str, "0:100,1:150,0:50")
    expect_equal(res$percent, 50)
    expect_equal(res$percent_combined, 50)
    expect_equal(res$n_spanning_pairs, 1L)
    expect_equal(length(res$cov_vec), 300L)
  })

  it("ignores a mate that straddles the segment boundary", {
    res <- call(s2 = 1280L, e2 = 1330L)  # crosses span_hi; mate1 also not inside
    expect_equal(res$rle_str, "0:300")
    expect_equal(res$percent, 0)
    expect_equal(res$n_spanning_pairs, 0L)
  })

  it("ignores a mate that lies entirely outside the segment", {
    res <- call(s1 = 500L, e1 = 649L, s2 = 700L, e2 = 849L)
    expect_equal(res$percent, 0)
    expect_equal(res$n_spanning_pairs, 0L)
  })

  it("skips FF-oriented pairs", {
    res <- call(orientation = "FF")
    expect_equal(res$percent, 0)
    expect_equal(res$n_spanning_pairs, 0L)
  })

  it("skips RR-oriented pairs", {
    res <- call(orientation = "RR")
    expect_equal(res$percent, 0)
    expect_equal(res$n_spanning_pairs, 0L)
  })

  it("keeps RF (everted) pairs", {
    res <- call(orientation = "RF")
    expect_equal(res$n_spanning_pairs, 1L)
    expect_equal(res$percent, 50)
  })

  it("skips pairs whose QNAME is not among the supporting reads", {
    res <- call(support_qnames = "other_read")
    expect_equal(res$percent, 0)
    expect_equal(res$n_spanning_pairs, 0L)
  })

  it("skips invalid (non-primary / low-MAPQ) pairs", {
    res <- call(valid = FALSE)
    expect_equal(res$percent, 0)
    expect_equal(res$n_spanning_pairs, 0L)
  })

  it("unions both mates when both fall inside the segment", {
    # mate1 [1000, 1149] and mate2 [1150, 1299] tile [1000, 1299].
    # Both contribute, giving complete coverage from one fragment.
    res <- call(s1 = 1000L, e1 = 1149L, s2 = 1150L, e2 = 1299L)
    expect_equal(res$rle_str, "1:300")
    expect_equal(res$percent, 100)
    expect_equal(res$n_spanning_pairs, 1L)
  })

  it("counts overlapping mates once per fragment at shared bases", {
    # mate1 covers 1-150 and mate2 covers 101-250. Shared 101-150
    # remains depth 1, rather than being double-counted as depth 2.
    res <- call(s1 = 1000L, e1 = 1149L, s2 = 1100L, e2 = 1249L)
    expect_equal(res$rle_str, "1:250,0:50")
    expect_equal(res$percent, 250 / 300 * 100)
    expect_equal(res$n_spanning_pairs, 1L)
  })

  it("clamps contributions to seq_len when the sequence is shorter than the span", {
    # mate [1150, 1250] -> local [151, 251], clamped to [151, 200] of 200
    res <- call(s2 = 1150L, e2 = 1250L, seq_len = 200L)
    expect_equal(res$rle_str, "0:150,1:50")
    expect_equal(res$percent, 25)
    expect_equal(res$n_spanning_pairs, 1L)
  })

  it("computes the combined percentage as the union with direct coverage", {
    direct <- c(rep(1L, 50), rep(0L, 250))   # direct covers 1-50
    res <- call(direct_cov = direct)          # PE covers 101-250
    expect_equal(res$percent, 50)
    expect_equal(res$percent_combined, 200 / 300 * 100)
  })

  it("combined coverage is 100% when direct already covers everything", {
    direct <- rep(1L, 300)
    res <- call(direct_cov = direct)
    expect_equal(res$percent_combined, 100)
  })

  it("combined coverage falls back to PE-only when direct_cov is not supplied", {
    res <- call()
    expect_equal(res$percent_combined, res$percent)
  })

  it("stacks depth across pairs covering the same region", {
    res <- psc(
      qname = c("r1", "r2"),
      s1 = c(900L, 900L), e1 = c(1049L, 1049L),
      s2 = c(1100L, 1150L), e2 = c(1249L, 1299L),
      orientation = c("FR", "FR"),
      valid = c(TRUE, TRUE),
      support_qnames = c("r1", "r2"),
      breakpoint = 1000L, span_len = 300L, seq_len = 300L
    )
    # r1 covers 101-250 (depth 1), r2 covers 151-300 (depth 1):
    # 101-150 at 1, 151-250 at 2, 251-300 at 1.
    expect_equal(res$rle_str, "0:100,1:50,2:100,1:50")
    expect_equal(res$n_spanning_pairs, 2L)
    # Union of the two mates' extents: 200 of 300 bases.
    expect_equal(res$percent, 200 / 300 * 100)
  })

  it("handles a mate covering the entire segment", {
    res <- call(s2 = 1000L, e2 = 1299L)
    expect_equal(res$rle_str, "1:300")
    expect_equal(res$percent, 100)
    expect_equal(res$n_spanning_pairs, 1L)
  })

  it("emits RLE strings in the existing ITDCoverageRLE format", {
    res <- call()
    # Existing format: value:length runs joined by commas, as produced by
    # the direct-coverage block of .compute_variant_metrics().
    parts <- strsplit(res$rle_str, ",", fixed = TRUE)[[1]]
    vals <- as.integer(sub(":.*$", "", parts))
    lens <- as.integer(sub("^.*:", "", parts))
    expect_equal(sum(lens), 300L)            # runs tile the whole vector
    expect_equal(sum(lens[vals > 0]), 150L)  # covered base count
  })
})

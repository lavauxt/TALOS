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

describe("candidate clustering and extension support", {
  cluster <- get_fn(".cluster_breakpoints")
  collapse <- get_fn(".collapse_extension_candidates_by_read")
  mode_extension <- get_fn(".mode_extension_length")

  it("bounds total cluster width instead of allowing single-link chaining", {
    got <- cluster(c(100, 108, 116, NA_real_), cluster_tolerance = 10)
    expect_length(got, 2L)
    expect_equal(got[[1L]], c(100, 108))
    expect_equal(got[[2L]], 116)
  })

  it("keeps equal breakpoints together when tolerance is zero", {
    got <- cluster(c(100, 100, 101), cluster_tolerance = 0)
    expect_length(got, 2L)
    expect_equal(got[[1L]], c(100, 100))
    expect_equal(got[[2L]], 101)
  })

  it("rejects invalid cluster tolerances", {
    expect_error(cluster(c(100, 101), -1), "non-negative")
    expect_error(cluster(c(100, 101), c(1, 2)), "single")
  })

  it("prefers anchor-backed length when cluster-level modes are tied", {
    got <- mode_extension(
      c(210L, 210L, 300L, 300L),
      c("kmer_fallback", "kmer_fallback", "anchor", "anchor")
    )
    expect_equal(got, 300L)
  })

  it("uses the smaller observed mode when tied sizes have no anchor support", {
    expect_equal(mode_extension(c(210L, 210L, 300L, 300L)), 210L)
  })

  it("prefers anchor-backed size when one QNAME has tied size hypotheses", {
    x <- data.frame(
      read_name = c("r1", "r1"),
      breakpoint_refined = c(100, 100),
      itdsize = c(210L, 300L),
      itdsize_source = c("kmer_fallback", "anchor"),
      stringsAsFactors = FALSE
    )
    got <- collapse(x)
    expect_equal(got$itdsize, 300L)
    expect_equal(got$itdsize_source, "anchor")
  })

  it("handles missing extension-source annotations when breaking ties", {
    x <- data.frame(
      read_name = c("r1", "r1"),
      breakpoint_refined = c(100, 100),
      itdsize = c(210L, 300L),
      itdsize_source = c(NA_character_, "kmer_fallback"),
      stringsAsFactors = FALSE
    )
    got <- collapse(x)
    expect_true(!is.na(got$itdsize))
    expect_equal(got$itdsize, 210L)
  })

  it("reduces multiple candidate rows to one modal length per QNAME", {
    x <- data.frame(
      read_name = c("r1", "r1", "r1", "r2", "r3", "r3"),
      breakpoint_refined = c(100, 100, 101, 100, 101, 101),
      itdsize = c(300L, 300L, 210L, 300L, 300L, 300L),
      itdsize_source = c("anchor", "anchor", "kmer_fallback", "kmer_fallback", "anchor", "anchor"),
      stringsAsFactors = FALSE
    )
    got <- collapse(x)
    expect_equal(nrow(got), 3L)
    expect_equal(got$read_name, c("r1", "r2", "r3"))
    expect_equal(got$itdsize, c(300L, 300L, 300L))
    expect_equal(got$itdsize_source, c("anchor", "kmer_fallback", "anchor"))
  })
})

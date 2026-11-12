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

describe("repeat-aware ITD coverage placement", {
  align <- get_fn(".align_reads_to_itd_template")
  ambiguity <- get_fn(".align_reads_to_itd_template_ambiguity")

  it("reports conservative and possible coverage for a repeated read placement", {
    unit <- "ATCCGTTCCTAATAAGGAATGGTGATTCCC"
    template <- paste0(unit, unit)
    reads <- unit
    best <- align(reads, template, min_pid = 95, min_aligned = 20L)
    res <- ambiguity(reads, template, min_pid = 95, min_aligned = 20L,
                     max_alternatives = 4L, score_delta = 0,
                     best_result = best, seed_k = 9L)
    expect_equal(res$multi_placement_read_ends, 1L)
    expect_equal(res$certain_percent, 0)
    expect_equal(res$possible_percent, 100)
    expect_equal(res$ambiguous_percent, 100)
  })

  it("reports when alternative-placement enumeration reaches its cap", {
    unit <- "ATCCGTTCCTAATAAGGAATGGTGATTCCC"
    template <- paste0(unit, unit, unit)
    best <- align(unit, template, min_pid = 95, min_aligned = 20L)
    res <- ambiguity(unit, template, min_pid = 95, min_aligned = 20L,
                     max_alternatives = 2L, score_delta = 0,
                     best_result = best, seed_k = 9L)
    expect_equal(res$alternative_search_limit_read_ends, 1L)
    expect_equal(res$multi_placement_read_ends, 1L)
  })

  it("uses the best alignment directly when its template seeds are unique", {
    template <- paste0(
      "AGACTTTCAAAGATATGCTGGGTAGAGGTCGAGGTTATTATTTGTTACCAATTCTCATTGTG",
      "TTTCGGAACTTGCGTTTTAGGTATGTCTTAGTGACTCTAAATACCAAGGCAGTCCTCG"
    )
    read <- substr(template, 1L, 80L)
    best <- align(read, template, min_pid = 95, min_aligned = 20L)
    res <- ambiguity(read, template, min_pid = 95, min_aligned = 20L,
                     best_result = best, seed_k = 9L)
    expect_equal(res$certain_cov, res$possible_cov)
    expect_equal(res$multi_placement_read_ends, 0L)
    expect_gt(res$certain_percent, 0)
  })
})

describe("compute_pe_itd_length support-restricted estimate", {
  estimate <- get_fn("compute_pe_itd_length")

  make_pairs <- function(qnames, first_pos, second_pos) {
    n <- length(qnames)
    first <- GenomicAlignments::GAlignments(
      seqnames = S4Vectors::Rle(factor(rep("chr1", n))),
      pos = as.integer(first_pos), cigar = rep("50M", n),
      seqlengths = c(chr1 = 50000L)
    )
    last <- GenomicAlignments::GAlignments(
      seqnames = S4Vectors::Rle(factor(rep("chr1", n))),
      pos = as.integer(second_pos), cigar = rep("50M", n),
      seqlengths = c(chr1 = 50000L)
    )
    S4Vectors::mcols(first)$mapq <- rep(60L, n)
    S4Vectors::mcols(last)$mapq <- rep(60L, n)
    S4Vectors::mcols(first)$flag <- rep(0L, n)
    S4Vectors::mcols(last)$flag <- rep(16L, n)
    GenomicAlignments::GAlignmentPairs(first, last, names = qnames)
  }

  it("reports all-pair and support-QNAME PE estimates separately", {
    support_names <- paste0("support", seq_len(6L))
    other_span_names <- paste0("other_span", seq_len(6L))
    background_names <- paste0("background", seq_len(10L))
    qnames <- c(support_names, other_span_names, background_names)
    first_pos <- c(rep(9800L, 6L), rep(9750L, 6L), rep(10100L, 10L))
    second_pos <- c(rep(10150L, 6L), rep(10199L, 6L), rep(10250L, 10L))
    pairs <- make_pairs(qnames, first_pos, second_pos)

    all_est <- estimate(pairs, breakpoint = 10000L, nominal_read_len = 150L)
    support_est <- estimate(pairs, breakpoint = 10000L, nominal_read_len = 150L,
                             support_qnames = support_names)

    expect_equal(all_est$n_spanning, 12L)
    expect_equal(support_est$n_spanning, 6L)
    expect_equal(support_est$n_spanning_all, 12L)
    expect_equal(support_est$n_support_qname_matched, 6L)
    expect_equal(all_est$length_pe, 250)
    expect_equal(support_est$length_pe, 200)
    expect_gt(support_est$n_support_insert_pairs, 0L)
    expect_true(support_est$support_insert_status %in% c("computed", "insufficient_support_insert_pairs"))
  })
})

describe("tandem-expanded ambiguity search", {
  expanded <- get_fn(".align_reads_to_tandem_expanded_template")
  it("projects placements from a two-copy template back to one ITD unit", {
    skip_if_not(requireNamespace("pwalign", quietly = TRUE))
    skip_if_not(requireNamespace("Biostrings", quietly = TRUE))
    unit <- "ATCCGTTCCTAATAAGGAATGGTGATTCCC"
    z <- expanded(unit, unit, min_pid = 95, min_aligned = 20L,
                   max_alternatives = 4L, score_delta = 0)
    expect_length(z$certain_cov, nchar(unit))
    expect_length(z$possible_cov, nchar(unit))
    expect_identical(z$template_mode, "single_copy_plus_two_sided_junction_screen")
  })
})

describe("error-tolerant tandem junction candidate screen", {
  expanded <- get_fn(".align_reads_to_tandem_expanded_template")

  it("retains a junction candidate with one substitution in its crossing seed", {
    skip_if_not(requireNamespace("pwalign", quietly = TRUE))
    skip_if_not(requireNamespace("Biostrings", quietly = TRUE))
    unit <- "ACGTTGCAAGTCCTAGGATCCGATGCTAAC"
    # 9-mer spanning the unit|unit junction, with one deliberate mismatch.
    junction_read <- paste0(substr(unit, nchar(unit) - 5L, nchar(unit)),
                            substr(unit, 1L, 8L))
    junction_read <- paste0(substr(junction_read, 1L, 4L),
                            ifelse(substr(junction_read, 5L, 5L) == "A", "C", "A"),
                            substr(junction_read, 6L, nchar(junction_read)))
    got <- expanded(junction_read, unit, min_pid = 70, min_aligned = 12L,
                    max_alternatives = 3L, score_delta = 1, seed_k = 9L)
    expect_gte(got$junction_candidate_read_ends, 1L)
    expect_identical(got$template_mode,
                     "single_copy_plus_selective_two_sided_junction_window")
  })
})

test_that("combined coverage depth is the per-base maximum, not gap-only fill", {
  fold <- get_fn(".fold_display_coverage")
  direct <- c(3L, 0L, 1L, 0L)
  pe <- c(8L, 5L, 0L, 2L)
  kmer <- c(0L, 7L, 4L, 9L)
  got <- fold(direct, pe, kmer)
  expect_equal(got, c(8L, 7L, 4L, 9L))
  expect_gte(sum(got), sum(direct))
  expect_gte(sum(got), sum(pe))
  expect_gte(sum(got), sum(kmer))
  expect_equal(fold(c(0L, 1L), c(4L, 0L), integer(0)), c(4L, 1L))
  expect_equal(fold(c(0L, 1L), c(4L, 0L), integer(0), direct_is_binary = TRUE),
               c(4L, 1L))
})

describe("selective tandem-junction screening", {
  find_candidates <- get_fn(".find_tandem_junction_candidates")

  it("requires a two-sided junction window rather than an isolated short seed", {
    template <- "ACGTTGCAACCTGATCGTAGGCTAACGATCGTACCTGAGCTTACGATCGGA"
    flank <- 12L
    junction <- paste0(substr(template, nchar(template) - flank + 1L,
                              nchar(template)),
                       substr(template, 1L, flank))
    read_spanning <- paste0("TTTT", junction, "AAAA")
    read_left_only <- paste0(substr(junction, 1L, flank), "TTTTGGGGCCCC")
    expect_equal(find_candidates(c(read_spanning, read_left_only), template,
                                 anchor_width = flank, max_mismatch = 0L), 1L)
    reverse_read <- as.character(Biostrings::reverseComplement(
      Biostrings::DNAString(read_spanning)))
    reverse_candidates <- find_candidates(reverse_read, template,
                                           anchor_width = flank, max_mismatch = 0L)
    expect_equal(reverse_candidates, 1L)
    expect_identical(attr(reverse_candidates, "orientation"), "reverse")
    forward_candidates <- find_candidates(read_spanning, template,
                                           anchor_width = flank, max_mismatch = 0L)
    expect_identical(attr(forward_candidates, "orientation"), "forward")
  })

  it("tolerates a small number of junction substitutions", {
    template <- "ACGTTGCAACCTGATCGTAGGCTAACGATCGTACCTGAGCTTACGATCGGA"
    flank <- 12L
    junction <- paste0(substr(template, nchar(template) - flank + 1L,
                              nchar(template)),
                       substr(template, 1L, flank))
    chars <- strsplit(junction, "", fixed = TRUE)[[1L]]
    chars[c(4L, 19L)] <- c("A", "C")
    read_spanning <- paste0("TTTT", paste(chars, collapse = ""), "AAAA")
    expect_equal(find_candidates(read_spanning, template,
                                 anchor_width = flank, max_mismatch = 2L), 1L)
  })
})

describe("per-read-end junction coverage merging", {
  merge_coverage <- get_fn(".merge_junction_read_coverage")

  it("adds distinct read ends but never counts one read end twice at a position", {
    base <- list(
      per_read_certain_positions = list(c(10L, 90L), 10L),
      per_read_possible_positions = list(c(10L, 90L), 10L)
    )
    junction <- list(
      per_read_certain_positions = list(seq_len(24L)),
      per_read_possible_positions = list(seq_len(24L))
    )
    out <- merge_coverage(base, junction, junction_idx = 1L,
                          template_length = 100L, anchor_width = 12L,
                          n_read_ends = 2L)
    expect_equal(out$certain_cov[10L], 2L)
    expect_equal(out$possible_cov[10L], 2L)
    expect_equal(out$certain_cov[90L], 1L)
    expect_equal(out$certain_cov[89L], 1L)
    expect_equal(out$certain_cov[1L], 1L)
    expect_equal(out$certain_cov[12L], 1L)
    expect_equal(out$certain_cov[100L], 1L)
  })
})

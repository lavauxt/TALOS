test_that("exon annotation accepts ALU insertion-site coordinates", {
  exons <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(100L, 200L), end = c(150L, 250L))
  )
  exons$exon_idx <- c(1L, 2L)
  alu_results <- data.frame(InsertionSite = c(120L, 180L))

  annotated <- .annotate_exonic_region(
    alu_results,
    exons,
    pos_col = "InsertionSite"
  )

  expect_equal(annotated$Region, c("exonic", "intronic"))
  expect_equal(annotated$ExonNumber, c(1L, NA_integer_))
})

test_that("exon annotation keeps ITD coordinate defaults", {
  exons <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = 100L, end = 150L)
  )
  results <- data.frame(GenomicPosition = 120L)

  expect_equal(
    .annotate_exonic_region(results, exons)$Region,
    "exonic"
  )
})

test_that("ALU clip matching scores both orientations consistently", {
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("pwalign")
  consensus <- Biostrings::DNAStringSet(
    setNames(paste(rep("A", 40L), collapse = ""), "AluTest")
  )

  hit <- .match_clip_to_alu(strrep("T", 30L), consensus, min_score = 0.6)

  expect_equal(hit$subtype, "AluTest")
  expect_equal(hit$strand, "-")
  expect_gte(hit$score, 0.6)
})

test_that("wild-type depth counts unique non-supporting reads at a breakpoint", {
  reads <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(100L, 105L, 100L, 100L),
                     end = c(110L, 110L, 110L, 110L))
  )
  wt_info <- list(
    gr = reads,
    qnames = c("support1", "support2", "wildtype1", "wildtype1"),
    is_primary = rep(TRUE, 4L)
  )

  expect_equal(
    .count_wildtype_at_position(wt_info, 105L,
                                support_qnames = c("support1", "support2")),
    1L
  )
})

test_that("ALU clusters deduplicate reads and compute a finite VAF", {
  candidates <- data.frame(
    read_name = c("support1", "support1", "support2", "support3"),
    genomic_pos = c(105L, 105L, 106L, 105L),
    alu_score = c(0.9, 0.7, 0.8, 0.85),
    is_reverse = c(FALSE, FALSE, TRUE, FALSE),
    poly_a_len = c(6L, 6L, 0L, 6L),
    alu_subtype = rep("AluTest", 4L),
    alu_strand = rep("+", 4L),
    aln_start = rep(1L, 4L),
    aln_end = rep(30L, 4L),
    clip_seq = rep(strrep("A", 30L), 4L),
    mapq = rep(60L, 4L)
  )
  reads <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = c(100L, 100L, 100L, 100L),
                             end = c(110L, 110L, 110L, 110L))
  )
  wt_info <- list(
    gr = reads,
    qnames = c("support1", "support2", "support3", "wildtype1"),
    is_primary = rep(TRUE, 4L)
  )
  config <- list(genomic_start = 100L)
  alu_seqs <- Biostrings::DNAStringSet(
    setNames(strrep("A", 40L), "AluTest")
  )

  result <- .summarise_alu_cluster(
    candidates, strrep("A", 30L), config, wt_info, alu_seqs, min_support = 3L
  )

  expect_equal(result$SupportingReads, 3L)
  expect_equal(result$WildtypeReads, 1L)
  expect_equal(result$AlleleFrequency, 0.75)
})

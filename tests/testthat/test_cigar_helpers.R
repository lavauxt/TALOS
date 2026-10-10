test_that("terminal soft clips are found beside hard clips", {
  expect_equal(
    .get_softclips("5H3S10M", "AAACCCCCCCCCC")$lead,
    "AAA"
  )
  expect_equal(
    .get_softclips("10M3S5H", "CCCCCCCCCCAAA")$trail,
    "AAA"
  )
  expect_equal(.get_softclips("3S10M", "AAACCCCCCCCCC")$lead, "AAA")
  expect_equal(.get_softclip_lengths("3S"), list(lead = 3L, trail = 0L))
  expect_equal(
    .get_softclips("10M", "CCCCCCCCCC"),
    list(lead = NA_character_, trail = NA_character_)
  )
  expect_equal(
    .has_softclip(c("5H3S10M", "10M3S5H", "10M")),
    c(TRUE, TRUE, FALSE)
  )
})

test_that("soft-clip filtering handles adjacent hard clips", {
  reads <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle("chr1", 2L),
    pos = c(100L, 200L),
    cigar = c("5H3S10M", "10M"),
    strand = S4Vectors::Rle(c("+", "+"))
  )
  S4Vectors::mcols(reads)$mapq <- c(60L, 60L)
  S4Vectors::mcols(reads)$flag <- c(0L, 0L)

  filtered <- .filter_reads_by_cigar(reads, min_mapq = 20L, min_ins_filter = 3L)

  expect_equal(length(filtered), 1L)
  expect_equal(GenomicAlignments::cigar(filtered), "5H3S10M")
})

test_that("CIGAR insertion parser returns each reference junction", {
  insertions <- .parse_cigar_insertions(
    "5M2I3M3I5M",
    read_start = 100L,
    genomic_start = 90L,
    ref_len = 100L,
    min_size = 2L
  )

  expect_equal(insertions$local_breakpoint, c(15L, 18L))
  expect_equal(insertions$length, c(2L, 3L))
})

test_that("CIGAR candidate extraction keeps separate insertions", {
  reads <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle("chr1", 1L),
    pos = 100L,
    cigar = "5M2I3M3I5M",
    strand = S4Vectors::Rle("+", 1L)
  )
  S4Vectors::mcols(reads)$qname <- "read1"
  S4Vectors::mcols(reads)$mapq <- 60L
  S4Vectors::mcols(reads)$flag <- 0L
  S4Vectors::mcols(reads)$seq <- Biostrings::DNAStringSet(strrep("A", 18L))

  candidates <- .extract_candidates_cigar(
    reads,
    genomic_start = 90L,
    ref_len = 100L,
    min_size = 1L
  )

  expect_equal(length(candidates), 2L)
  expect_equal(
    vapply(candidates, `[[`, integer(1L), "local_breakpoint"),
    c(15L, 18L)
  )
  expect_equal(vapply(candidates, `[[`, integer(1L), "length"), c(2L, 3L))
})

test_that("CIGAR extraction accepts terminal clips beside hard clips", {
  reads <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle("chr1", 1L),
    pos = 100L,
    cigar = "5H3S10M",
    strand = S4Vectors::Rle("+", 1L)
  )
  S4Vectors::mcols(reads)$qname <- "read_hardclip"
  S4Vectors::mcols(reads)$mapq <- 60L
  S4Vectors::mcols(reads)$flag <- 0L
  S4Vectors::mcols(reads)$seq <- Biostrings::DNAStringSet("AAACCCCCCCCCC")

  candidates <- .extract_candidates_cigar(
    reads,
    genomic_start = 90L,
    ref_len = 100L
  )

  expect_length(candidates, 1L)
  expect_equal(candidates[[1L]]$local_breakpoint, 11L)
  expect_equal(
    .get_softclips(candidates[[1L]]$cigar, candidates[[1L]]$read_seq)$lead,
    "AAA"
  )
})

test_that("CIGAR parsing preserves short insertion candidates", {
  insertions <- .parse_cigar_insertions(
    "5M1I5M",
    read_start = 100L,
    genomic_start = 90L,
    ref_len = 100L
  )

  expect_equal(insertions$length, 1L)
})

test_that("CIGAR insertion candidates use their own minimum size", {
  reads <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle("chr1", 1L),
    pos = 100L,
    cigar = "5M1I5M",
    strand = S4Vectors::Rle("+", 1L)
  )
  S4Vectors::mcols(reads)$qname <- "read_short_ins"
  S4Vectors::mcols(reads)$mapq <- 60L
  S4Vectors::mcols(reads)$flag <- 0L
  S4Vectors::mcols(reads)$seq <- Biostrings::DNAStringSet(strrep("A", 11L))
  names(reads) <- "read_short_ins"

  candidates <- .extract_candidates_cigar(
    reads,
    genomic_start = 90L,
    ref_len = 100L,
    min_size = 10L
  )

  expect_length(candidates, 1L)
  expect_equal(candidates[[1L]]$length, 1L)
})

test_that("CIGAR insertions account for deletions and skips", {
  insertions <- .parse_cigar_insertions(
    "5M2D3M2I4N5M",
    read_start = 100L,
    genomic_start = 90L,
    ref_len = 100L
  )

  expect_equal(insertions$local_breakpoint, 20L)
})

test_that("coverage helper returns the same span summary and depth ratio", {
  coverage <- c(rep(10L, 20L), rep(20L, 10L), rep(10L, 20L))
  summary <- .compute_span_coverage_summary(
    coverage,
    chrom = "chr1",
    start = 21L,
    span_len = 10L,
    flank = 10L
  )

  expect_equal(summary$min_cov, 20)
  expect_equal(summary$mean_cov, 20)
  expect_equal(summary$depth_fold_change, 2)
  expect_equal(
    compute_span_depth_fold_change(
      coverage,
      "chr1",
      21L,
      10L,
      flank = 10L
    ),
    summary$depth_fold_change
  )
})

test_that("exonic breakpoint lookup uses each available exon field", {
  exons <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(100L, 200L), width = 20L)
  )

  expect_identical(.is_breakpoint_exonic(105L, "chr1", exons), TRUE)
  expect_identical(.is_breakpoint_exonic(150L, "chr1", exons), FALSE)
  expect_identical(.get_gene_exons(list(target_exons = exons)), exons)
  expect_identical(.get_gene_exons(list(all_exons = exons)), exons)
})

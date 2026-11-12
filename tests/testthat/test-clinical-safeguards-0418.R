test_that("raw ITD AF remains distinct from exploratory size-corrected AF", {
  wt_info <- list(
    is_primary = rep(TRUE, 10L),
    gr = GenomicRanges::GRanges(
      seqnames = rep("chr1", 10L),
      ranges = IRanges::IRanges(start = rep(90L, 10L), end = rep(110L, 10L))
    ),
    qnames = paste0("wt", seq_len(10L))
  )
  got <- TALOS:::.calculate_vaf_and_depth(
    raw_support = 10L,
    wildtype_info = wt_info,
    len_specific_bp = 100L,
    support_qnames = paste0("itd", seq_len(10L)),
    best_len = 300L,
    nominal_read_len = 150L,
    max_correction = 2
  )
  expect_equal(got$raw_support, 10L)
  expect_equal(got$corrected_support, 15L)
  expect_equal(got$size_correction_factor, 1.5)
  expect_equal(got$raw_depth, 20L)
  expect_equal(got$depth_at_bp, 20L)
  expect_equal(got$raw_af, 0.5)
  expect_equal(got$corrected_af, 15 / 25)
})

test_that("size correction defaults to no inflation and rejects invalid caps", {
  wt_info <- list(is_primary = logical(), gr = GenomicRanges::GRanges(), qnames = character())
  got <- TALOS:::.calculate_vaf_and_depth(
    10L, wt_info, 100L, character(), 300L, 150L, 1
  )
  expect_equal(got$size_correction_factor, 1)
  expect_equal(got$corrected_support, 10L)
  expect_error(
    TALOS:::.calculate_vaf_and_depth(10L, wt_info, 100L, character(), 300L, 150L, NA_real_),
    "max_correction"
  )
})

test_that("CIGAR insertion extraction honours the caller's minimum size", {
  one <- TALOS:::.parse_cigar_insertions("50M5I50M", 100L, 1L, 500L, min_size = 1L)
  ten <- TALOS:::.parse_cigar_insertions("50M5I50M", 100L, 1L, 500L, min_size = 10L)
  expect_equal(nrow(one), 1L)
  expect_equal(one$length, 5L)
  expect_equal(nrow(ten), 0L)
})

test_that("microhomology is measured against the flank adjacent to each read anchor", {
  ref <- paste0("AAAAAATTAC", "ACGT", paste(rep("C", 30L), collapse = ""))
  rows <- data.frame(
    cigar = c("6S20M", "20M4S"),
    read_seq = c("GGTTAC", paste0(strrep("A", 20L), "ACGT")),
    read_start = c(110L, 90L),
    read_end = c(129L, 109L),
    stringsAsFactors = FALSE
  )
  expect_equal(TALOS:::compute_microhomology(rows, ref, 110L, 100L), 4L)
})

test_that("VCF header omits unknown contig length instead of writing invalid '?'", {
  hdr <- TALOS::build_vcf_header("S1", "chr13", genome_build = "hg38", contig_length = NULL)
  expect_true(any(grepl("^##contig=<ID=chr13>$", hdr)))
  expect_false(any(grepl("length=\\?", hdr)))
})

test_that("the public ITD detector defaults to no size-bias inflation", {
  expect_equal(eval(formals(TALOS::detect_itd)$max_correction), 1.0)
})

test_that("report coverage loader returns its actual coordinate origin", {
  precomputed <- S4Vectors::Rle(rep(2L, 100L))
  got <- TALOS:::.plot_coverage(NULL, list(), from = 10L, to = 20L,
                                precomputed_cov = precomputed)
  expect_equal(got$start, 1L)
  expect_identical(got$cov, precomputed)
  empty <- TALOS:::.plot_coverage(NULL, list(), from = 10L, to = 20L)
  expect_equal(empty$start, 10L)
  expect_null(empty$cov)
})

test_that("variant detail summaries label raw support separately from corrected diagnostics", {
  row <- data.frame(
    Sample = "S1", Gene = "FLT3", AlleleFrequency = 0.2,
    CorrectedAlleleFrequency = 0.4, SupportingReads = 10L,
    AdjustedSupportingReads = 20L, WildtypeReads = 40L,
    DepthAtBreakpoint = 50L, CorrectedDepthAtBreakpoint = 60L,
    SizeCorrectionFactor = 2, InFrame = TRUE,
    stringsAsFactors = FALSE
  )
  fields <- TALOS:::.variant_fields(row, genome_build = "hg38", mode = "itd")
  expect_identical(unname(fields[["Supporting Reads - raw"]]), "10")
  expect_identical(unname(fields[["Supporting Reads - size-corrected diagnostic"]]), "20")
  expect_identical(unname(fields[["VAF (%) - raw"]]), "20.0")
  expect_identical(unname(fields[["VAF (%) - size-corrected diagnostic"]]), "40.0")
  expect_identical(unname(fields[["In frame (diagnostic)"]]), "TRUE")
})

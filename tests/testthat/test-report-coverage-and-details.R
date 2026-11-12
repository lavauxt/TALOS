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

test_that("report coverage falls back to the source-specific RLE vectors", {
  build_display <- get_fn(".build_itd_cov_display")
  row <- data.frame(
    GenomicPosition = 100L,
    ITDCoverageRLE = NA_character_,
    ITDCoverageRLE_Direct = "2:2,0:2",
    ITDCoverageRLE_PE = "0:2,3:2",
    ITDCoverageRLE_Kmer = "0:3,4:1",
    stringsAsFactors = FALSE
  )
  axis <- list(transform_pos = function(x) as.numeric(x))
  out <- build_display(row, 103L, axis)
  # PE reference-coordinate projection is diagnostic only and cannot fill
  # positions in the combined allele-specific display track.
  expect_equal(out$count, c(2, 2, 0, 4))
  expect_equal(out$raw_idx, 1:4)
  expect_equal(out$disp_pos, 100:103)
})

test_that("merged intervals retain their per-base RLE coverage", {
  merge_rle <- get_fn(".merge_coverage_rle_column")
  rows <- data.frame(
    GenomicPosition = c(100L, 102L),
    ITDCoverageRLE = c("1:3", "2:3"),
    stringsAsFactors = FALSE
  )
  expect_equal(
    merge_rle(rows, "ITDCoverageRLE", merged_start = 100L, merged_len = 5L),
    "1:2,2:3"
  )
})


test_that("genome-view plot context includes flanking exon 16", {
  build_context <- get_fn(".build_plot_exon_context")
  # Negative-strand transcript order: E16 is exon_idx 1 and is immediately
  # upstream of the configured E13-E15 target set in transcript order.
  all_exons <- data.frame(
    start = c(28602315L, 28608024L, 28608219L, 28608438L),
    end = c(28602425L, 28608128L, 28608351L, 28608544L),
    exon_idx = 1:4
  )
  target_exons <- all_exons[2:4, , drop = FALSE]
  ctx <- build_context(
    list(all_exons = all_exons, target_exons = target_exons),
    flank_exons = 1L,
    flank_bp = 50L,
    span = c(28608047L, 28608346L)
  )

  expect_true(1L %in% ctx$exons_plot$exon_idx)
  expect_lte(ctx$plot_start, 28602315L)
  expect_gte(ctx$plot_end, 28608544L)
})


test_that("variant details disclose a coverage-only reference template", {
  variant_fields <- get_fn(".variant_fields")
  row <- data.frame(
    Sample = "HORIZON",
    Gene = "FLT3",
    CoverageTemplateSource = "reference_segment_for_coverage_only",
    SequenceSource = "partial_read_unresolved",
    stringsAsFactors = FALSE
  )
  fields <- variant_fields(row, genome_build = "hg19", mode = "itd")
  expect_identical(
    unname(fields[["Coverage Template Source"]]),
    "reference_segment_for_coverage_only"
  )
  expect_identical(unname(fields[["Sequence Source"]]), "partial_read_unresolved")
})


test_that("BAM report coverage is indexed relative to the full plot window", {
  coverage_window <- get_fn(".coverage_window_from_ranges")
  sample_coverage <- get_fn(".sample_compressed_coverage")
  from <- 28602310L
  to <- 28608052L
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr13", "chr13", "chr12"),
    ranges = IRanges::IRanges(
      start = c(28602315L, 28608047L, 28602315L),
      width = c(3L, 4L, 10L)
    )
  )

  cov <- coverage_window(gr, "chr13", from, to)
  expect_length(cov, to - from + 1L)
  expect_equal(as.numeric(cov[6L:8L]), c(1, 1, 1))
  idx_e15 <- 28608047L - from + 1L
  expect_equal(as.numeric(cov[idx_e15:(idx_e15 + 3L)]), rep(1, 4L))
  expect_equal(as.numeric(cov[1L]), 0)

  axis <- list(
    segments = data.frame(
      raw_start = c(28602315L, 28608047L),
      raw_end = c(28602317L, 28608050L),
      type = c("exon", "exon")
    ),
    transform_pos = function(x) as.numeric(x)
  )
  plotted <- sample_coverage(cov, axis, coverage_start = from)
  expect_equal(plotted$depth, rep(1, 7L))
})


test_that("an empty report-window query returns a zero Rle over the requested window", {
  coverage_window <- get_fn(".coverage_window_from_ranges")
  empty_gr <- GenomicRanges::GRanges()
  cov <- coverage_window(empty_gr, "chr13", 100L, 199L)
  expect_length(cov, 100L)
  expect_equal(sum(as.numeric(cov)), 0)
})


test_that("an explicit zero ITD coverage RLE remains drawable", {
  build_display <- get_fn(".build_itd_cov_display")
  row <- data.frame(
    GenomicPosition = 500L,
    ITDCoverageRLE = "0:300",
    stringsAsFactors = FALSE
  )
  axis <- list(transform_pos = function(x) as.numeric(x))
  out <- build_display(row, 799L, axis)
  expect_equal(nrow(out), 300L)
  expect_equal(out$count, rep(0, 300L))
  expect_equal(range(out$disp_pos), c(500, 799))
})


test_that("HTML variant detail payload retains every raw result variable", {
  raw_values <- get_fn(".raw_variant_row_values")
  row <- data.frame(
    Sample = "HORIZON",
    ITDReadCoverage = 66,
    ITDCoverageReadsProcessed = 298L,
    CoverageTemplateSource = "reference_segment_for_coverage_only",
    ITD_Sequence = NA_character_,
    InternalDiagnostic = "retained",
    stringsAsFactors = FALSE
  )
  got <- raw_values(row)
  expect_identical(got$keys, names(row))
  expect_identical(got$values[match("InternalDiagnostic", got$keys)], "retained")
  expect_identical(got$values[match("ITD_Sequence", got$keys)], NA_character_)
})


test_that("Variant Details exposes the new read-end and fragment-level k-mer metrics", {
  variant_fields <- get_fn(".variant_fields")
  row <- data.frame(
    Sample = "HORIZON", Gene = "FLT3",
    ITDKmerReads = 12L, ITDKmerTestedReads = 20L,
    ITDKmerTestedFragments = 15L, ITDKmerMatchedFragments = 10L,
    ITDKmerMatchedFragmentPct = 66.7,
    stringsAsFactors = FALSE
  )
  fields <- variant_fields(row, genome_build = "hg19", mode = "itd")
  expect_true(all(c("K-mer Matched Read Ends", "K-mer Tested Read Ends",
                    "K-mer Tested Fragments", "K-mer Matched Fragments",
                    "K-mer Matched Fragments (%)") %in% names(fields)))
})

test_that("each result table ships its own column chooser, quick filters and density toggle", {
  js <- getFromNamespace(".TALOS_REPORT_JS", "TALOS")
  # Reusable toolbar widgets instantiated by all three result tables.
  expect_match(js, "function makeColumnChooser(", fixed = TRUE)
  expect_match(js, "function makeQuickFilters(", fixed = TRUE)
  expect_match(js, "function makeDensitySeg(", fixed = TRUE)
  # Independent persisted selections per table.
  expect_match(js, "talos-cols-details-", fixed = TRUE)
  expect_match(js, "talos-cols-all-", fixed = TRUE)
  # The legacy per-variant cards still follow the Results selection.
  expect_match(js, "detailColumnSelection", fixed = TRUE)
  expect_match(js, "refreshDetailTables", fixed = TRUE)
  expect_match(js, "Find any variable", fixed = TRUE)
  expect_match(js, "Select all", fixed = TRUE)
  expect_match(js, "T.cols.filter", fixed = TRUE)
  # The details and all-values tables gain the Results-table features:
  # variable selection, compact mode, quick filters, legend, CSV export.
  expect_match(js, "groupLegendEl(chooser.order())", fixed = TRUE)
  expect_match(js, "view.setDensity(currentDensity())", fixed = TRUE)
  expect_match(js, "TALOS_variant_details.csv", fixed = TRUE)
  expect_match(js, "TALOS_all_result_values.csv", fixed = TRUE)
})

test_that("the PDF compact variant-details set is a curated subset of the full fields", {
  compact_fn <- get_fn(".compact_variant_fields")
  full_fn <- get_fn(".variant_fields")
  row <- data.frame(
    Sample = "HORIZON", Gene = "FLT3", GenomicPosition = 28608047L, Length = 300L,
    AlleleFrequency = 0.041, SupportingReads = 84L, WildtypeReads = 1966L,
    DepthAtBreakpoint = 2050L, HGVS_cDNA = "c.1710_1919dup",
    HGVS_Protein = "p.Phe570_Gln640dup", Hotspot = TRUE, HotspotName = "FLT3_JM_ITD_region",
    LengthPE_AllPairs_Eligibility = paste("Very long eligibility note that used to wrap",
      "across several drawn table lines and pushed the old PDF summary pages out of",
      "the page frame; it must not appear in the compact set."),
    stringsAsFactors = FALSE
  )
  full <- full_fn(row, genome_build = "hg19", mode = "itd")
  compact <- compact_fn(row, genome_build = "hg19", mode = "itd")
  expect_true(all(names(compact) %in% names(full)))
  expect_identical(unname(compact[["HGVS cDNA"]]), "c.1710_1919dup")
  expect_identical(unname(compact[["VAF (%) - raw"]]), "4.1")
  # The compact set keeps the clinical anchors but drops the long diagnostics.
  expect_true(all(c("Sample", "Gene", "Genomic Position", "Length (bp)", "VAF (%) - raw",
                    "Supporting Reads - raw", "Wildtype Reads", "HGVS cDNA", "HGVS Protein",
                    "Hotspot", "Evidence") %in% names(compact)))
  expect_false("PE all-pair eligibility" %in% names(compact))
  expect_false("PE support-insert eligibility" %in% names(compact))
  expect_lte(length(compact), 36L)
  # ALU mode keeps its own curated labels.
  alu_row <- data.frame(Sample = "S1", Gene = "ALU", stringsAsFactors = FALSE)
  alu <- compact_fn(alu_row, genome_build = "hg19", mode = "alu")
  expect_true("VAF (%)" %in% names(alu))
  expect_false("VAF (%) - raw" %in% names(alu))
})

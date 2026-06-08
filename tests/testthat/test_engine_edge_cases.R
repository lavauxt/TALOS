# tests/testthat/test_engine_edge_cases.R

library(testthat)
library(GenomicRanges)
library(IRanges)

# ----------------------------------------------------------------------
# 1. .filter_reads_by_cigar() – missing CIGAR operation columns
# ----------------------------------------------------------------------
test_that(".filter_reads_by_cigar handles reads without soft clips", {
  # Create minimal GAlignments object with no 'S' in CIGAR
  ga <- GenomicAlignments::GAlignments(
    seqnames = Rle("chr1"),
    pos = c(100L, 200L),
    cigar = c("50M", "30M"),
    strand = Rle("+"),
    names = c("read1", "read2")
  )
  S4Vectors::mcols(ga)$mapq <- c(30L, 30L)
  S4Vectors::mcols(ga)$flag <- c(0L, 0L)
  
  filtered <- .filter_reads_by_cigar(ga, min_mapq = 20, min_ins_filter = 3)
  expect_equal(length(filtered), 0L)  # No softclip, no insertion -> filtered out
  
  # Add a read with insertion
  ga2 <- ga
  GenomicAlignments::cigar(ga2)[1] <- "50M2I30M"
  S4Vectors::mcols(ga2)$flag <- 0L
  filtered2 <- .filter_reads_by_cigar(ga2, min_mapq = 20, min_ins_filter = 2)
  expect_equal(length(filtered2), 1L)  # insertion passes filter
})

# ----------------------------------------------------------------------
# 2. .extract_candidates_standard() – NA sequences and short reads
# ----------------------------------------------------------------------
test_that(".extract_candidates_standard skips NA reads and short reads", {
  ref_dna <- paste(rep("A", 500), collapse = "")
  ref_kmers <- .prepare_kmers(ref_dna, k = 11)
  
  # Mock GAlignments with one NA seq, one too short
  ga <- GenomicAlignments::GAlignments(
    seqnames = Rle("chr1"),
    pos = c(100L, 200L),
    cigar = c("150M", "5M"),
    strand = Rle("+")
  )
  S4Vectors::mcols(ga)$seq <- DNAStringSet(c(NA, "AAAAA"))
  S4Vectors::mcols(ga)$mapq <- c(30L, 30L)
  S4Vectors::mcols(ga)$flag <- c(0L, 0L)
  names(ga) <- c("read1", "read2")
  
  cand <- .extract_candidates_standard(ga, ref_kmers, ptd_mode = FALSE,
                                        min_size = 10, max_missing_kmers = 0.5,
                                        refine_bp = FALSE, use_cigar_bp = FALSE,
                                        genomic_start = 1L, ref_dna = ref_dna)
  expect_length(cand, 0L)  # both rejected
})

# ----------------------------------------------------------------------
# 3. .cluster_breakpoints() – edge cases
# ----------------------------------------------------------------------
test_that(".cluster_breakpoints works with empty, single, and exact tolerance", {
  expect_equal(.cluster_breakpoints(numeric(0), 10), list())
  expect_equal(.cluster_breakpoints(c(500), 10), list(c(500)))
  bp <- c(100, 105, 200, 210, 215)
  clusters <- .cluster_breakpoints(bp, 10)
  expect_equal(length(clusters), 2L)
  expect_equal(clusters[[1]], c(100, 105))
  expect_equal(clusters[[2]], c(200, 210, 215))
})

# ----------------------------------------------------------------------
# 4. compute_hgvs_annotations() – with full_cdna missing or invalid
# ----------------------------------------------------------------------
test_that("compute_hgvs_annotations returns NA when transcript missing", {
  cfg <- list(transcript = NA_character_, full_cdna = "ATGC", all_exons = GRanges(),
              cds_offset = 0L)
  res <- compute_hgvs_annotations(cfg, genomic_pos = 1000, dup_len = 30)
  expect_equal(res$c_notation, "no transcript reference")
  expect_equal(res$p_notation, NA_character_)
})

test_that("compute_hgvs_annotations rebuilds full_cdna from exons if needed", {
  # Create a simple two-exon gene on + strand
  exons <- GRanges("chr1", IRanges(c(100, 200), c(150, 250)), strand = "+")
  exons$exon_idx <- 1:2
  bs <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19  # or mock
  cfg <- list(
    transcript = "NM_001",
    all_exons = exons,
    bsgenome_obj = bs,
    strand = "+",
    cds_offset = 0L,
    full_cdna = NULL   # force rebuild
  )
  res <- compute_hgvs_annotations(cfg, genomic_pos = 120, dup_len = 10)
  expect_true(!is.na(res$c_notation) && nchar(res$c_notation) > 0)
})

# ----------------------------------------------------------------------
# 5. detect_orientation() – fallback when pwalign missing
# ----------------------------------------------------------------------
test_that("detect_orientation uses exact match when alignment packages missing", {
  # Simulate missing pwalign (temporarily unload)
  with_mocked_bindings(
    requireNamespace = function(pkg, ...) pkg != "pwalign",
    {
      res <- detect_orientation("ATCG", "ATCG")
      expect_equal(res, "tandem")
      res2 <- detect_orientation("ATCG", as.character(reverseComplement(DNAString("ATCG"))))
      expect_equal(res2, "inverted")
      res3 <- detect_orientation("ATCG", "TTTT")
      expect_true(is.na(res3))
    }
  )
})

# ----------------------------------------------------------------------
# 6. PTD plotting – zero-length duplication handling (mock Gviz)
# ----------------------------------------------------------------------
test_that("plot_talos_report skips zero-length ITDs or adjusts", {
  itd_row <- data.frame(
    GenomicPosition = 5000, Length = 0, AlleleFrequency = 0.05,
    Orientation = "tandem", stringsAsFactors = FALSE
  )
  gene_config <- list(
    chrom = "chr1", genomic_start = 4000, genomic_end = 6000,
    target_exons = GRanges(), build = "hg19", gene = "TEST"
  )
  # We can't test graphics output, but ensure function does not error
  # Replace Gviz functions with mocks to avoid actual PDF creation
  with_mocked_bindings(
    pdf = function(...) NULL,
    dev.off = function(...) NULL,
    plotTracks = function(...) NULL,
    grid.text = function(...) NULL,
    grid.newpage = function(...) NULL,
    {
      expect_error(
        plot_talos_report(itd_row, gene_config, bam_path = "test_data/mini_flt3.bam",
                          sample_name = "sample", output_pdf = tempfile()),
        NA   # no error
      )
    }
  )
})

# We need a mini bam with FLT3
test_that("detect_itd runs without crashing on minimal BAM", {
  skip_if_not(file.exists("test_data/mini_flt3.bam"))
  cfg <- get_gene_config("FLT3", build = "hg19", use_db = FALSE)
  res <- detect_itd("test_data/mini_flt3.bam", cfg,
                    min_support = 2, verbose = FALSE,
                    write_vcf = FALSE, plot = FALSE)
  expect_s3_class(res, "data.frame")
  if (nrow(res) > 0) {
    expect_true(all(c("GenomicPosition", "Length", "AlleleFrequency") %in% names(res)))
  }
})
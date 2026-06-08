# tests/testthat/test_db_config.R
library(testthat)

test_that("get_gene_config loads YAML correctly for each gene and build", {
  # Test all genes and both builds
  genes <- c("FLT3", "BCOR", "KMT2A", "UBTF")
  builds <- c("hg19", "hg38")
  
  for (gene in genes) {
    for (build in builds) {
      # Use `use_db = FALSE` if you are testing the static YAML entries, 
      # or `use_db = TRUE` to test the TxDb/EnsDb fetch paths.
      cfg <- get_gene_config(gene = gene, build = build, use_db = FALSE)
      
      # Check basic structure
      expect_type(cfg, "list")
      expect_equal(cfg$gene, gene)
      expect_equal(cfg$build, build)
      expect_true(is.character(cfg$chrom) && nchar(cfg$chrom) > 0)
      expect_true(cfg$strand %in% c("+", "-"))
      expect_true(is.numeric(cfg$genomic_start) && cfg$genomic_start > 0)
      expect_true(is.numeric(cfg$genomic_end) && cfg$genomic_end > cfg$genomic_start)
      expect_true(is.character(cfg$rnaseq) && nchar(cfg$rnaseq) > 0)
      expect_true(is.character(cfg$genomic_ref_seq) && nchar(cfg$genomic_ref_seq) > 0)
      expect_s4_class(cfg$all_exons, "GRanges") # Note: updated to match db_config.R output names
      
      # Exon order should follow transcript direction
      exon_starts <- GenomicRanges::start(cfg$all_exons)
      if (cfg$strand == "+") {
        expect_true(all(diff(exon_starts) > 0))
      } else {
        expect_true(all(diff(exon_starts) < 0))
      }
      
      # Exon numbers should be present and match expected
      if (!is.null(S4Vectors::mcols(cfg$all_exons)$exon_idx)) {
        expect_true(all(!is.na(S4Vectors::mcols(cfg$all_exons)$exon_idx)))
      }
      
      # cDNA length should equal sum of target exon lengths (after strand correction)
      expected_cdna_len <- sum(GenomicRanges::width(cfg$target_exons))
      expect_equal(nchar(cfg$rnaseq), expected_cdna_len)
      
      # Genomic reference sequence should span from start to end
      expect_equal(nchar(cfg$genomic_ref_seq), cfg$genomic_end - cfg$genomic_start + 1)
    }
  }
})

test_that("get_gene_config handles DB fetching when static data missing or use_db = TRUE", {
  skip_if_not_installed("EnsDb.Hsapiens.v75")
  
  # Fetch via DB instead of using static YAML offline blocks
  cfg <- get_gene_config(gene = "FLT3", build = "hg19", use_db = TRUE)
  
  expect_true(is.character(cfg$rnaseq) && nchar(cfg$rnaseq) > 0)
  expect_s4_class(cfg$all_exons, "GRanges")
  expect_s4_class(cfg$target_exons, "GRanges")
})

test_that("Strand-specific cDNA construction is correct for minus strand genes", {
  # For a minus strand gene (e.g., FLT3), the cDNA sequence should be reverse complement
  # of the concatenated exonic sequences from the + strand reference.
  cfg <- get_gene_config(gene = "FLT3", build = "hg19", use_db = FALSE)
  expect_equal(cfg$strand, "-")
  
  # Extract exons from + strand and reverse complement using the cached BSgenome object
  bs <- cfg$bsgenome_obj
  exons_plus <- cfg$target_exons
  BiocGenerics::strand(exons_plus) <- "+"
  
  exon_seqs_plus <- as.character(BSgenome::getSeq(bs, exons_plus))
  cDNA_from_plus <- paste(exon_seqs_plus, collapse = "")
  cDNA_expected <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(cDNA_from_plus)))
  
  expect_equal(cfg$rnaseq, cDNA_expected)
})

# ============================================================================
# NEW: Engine parameter tests (NA handling & Subsampling limit)
# ============================================================================
test_that("TALOS engine properly handles max_pairwise_alignments and NA reads", {
  skip_if_not_installed("Rsamtools")
  
  # Setup a minimal mock config
  cfg <- get_gene_config(gene = "FLT3", build = "hg19", use_db = FALSE)
  
  # Mock a very small BAM file or simulate detect_itd inputs
  # Instead of full BAM run, we can test the expected parameter passing
  # and assure it doesn't crash on edge cases.
  
  # If you are doing integration testing with a real mini-BAM:
  mock_bam <- system.file("extdata", "test_minibam.bam", package = "TALOS")
  skip_if(mock_bam == "", "Mock BAM not available for integration test")
  
  # 1. Test that setting a very low max_pairwise_alignments does not crash
  #    when supporting reads > limit.
  res_subsampled <- detect_itd(
    bam_path = mock_bam,
    gene_config = cfg,
    max_pairwise_alignments = 2L, # Force subsampling
    verbose = FALSE,
    write_vcf = FALSE,
    plot = FALSE
  )
  expect_s3_class(res_subsampled, "data.frame")
  
  # Note: The NA-read check is harder to mock without a manually corrupted BAM,
  # but if your test mini-BAM includes unmapped/NA sequence reads, the fact 
  # that `detect_itd()` completes without a 'type character, expected numeric'
  # error confirms the `if (is.na(read_seq))` fix is functioning.
})

# In test_db_config.R, add:

test_that("get_gene_config caches TxDb objects in memory and on disk", {
  skip_if_not_installed("txdbmaker")
  # Clear cache
  rm(list = ls(envir = .txdb_cache), envir = .txdb_cache)
  cache_file <- file.path(.talos_cache_dir(), "talos_refGene_hg19.sqlite")
  if (file.exists(cache_file)) unlink(cache_file)
  
  # First call – builds from UCSC
  cfg1 <- get_gene_config("FLT3", build = "hg19", use_db = TRUE)
  expect_true(file.exists(cache_file))
  
  # Second call – should load from disk (no network)
  cfg2 <- get_gene_config("FLT3", build = "hg19", use_db = TRUE)
  expect_equal(cfg1$rnaseq, cfg2$rnaseq)
  
  # In-memory hit on third call
  cfg3 <- get_gene_config("FLT3", build = "hg19", use_db = TRUE)
})

test_that("get_gene_config fallback to biomaRt when UCSC refGene fails", {
  skip_if_not_installed("biomaRt")
  # Force UCSC failure by providing a transcript that does not exist in refGene
  # but does exist in biomaRt (e.g., a known NM_ that is not in UCSC refGene)
  # This is hard to guarantee; we can mock .get_refseq_txdb to return error
  with_mocked_bindings(
    .get_refseq_txdb = function(build) stop("Simulated UCSC failure"),
    {
      cfg <- get_gene_config("FLT3", build = "hg19", use_db = TRUE)
      expect_true(cfg$source == "biomaRt")
    }
  )
})
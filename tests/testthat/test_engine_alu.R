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

  forward_hit <- .match_clip_to_alu(strrep("A", 30L), consensus,
                                    min_score = 0.6)
  expect_equal(forward_hit$strand, "+")
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

test_that("bundled UBTF hotspots match NM_014233.4 exon 13 by build", {
  hotspots <- utils::read.csv(
    testthat::test_path("../../inst/extdata/hotspots.csv"),
    stringsAsFactors = FALSE
  )
  ubtf <- hotspots[hotspots$Gene == "UBTF", ]
  expect_equal(ubtf$Build, c("hg19", "hg38"))
  expect_equal(ubtf$Start, c(42288160L, 44210792L))
  expect_equal(ubtf$End, c(42288315L, 44210947L))
})

test_that("bundled BCOR hg38 hotspot matches NM_001123385.2 exon 15", {
  hotspots <- utils::read.csv(
    testthat::test_path("../../inst/extdata/hotspots.csv"),
    stringsAsFactors = FALSE
  )
  bcor_hg38 <- hotspots[hotspots$Gene == "BCOR" & hotspots$Build == "hg38", ]

  expect_equal(bcor_hg38$Start, 40051246L)
  expect_equal(bcor_hg38$End, 40052400L)
})

test_that("ALU YAML aliases retain the canonical hotspot gene symbols", {
  config <- yaml::read_yaml(
    testthat::test_path("../../inst/extdata/gene_config.yaml")
  )
  expect_equal(config$FLT3_ALU$gene_symbol, "FLT3")
  expect_equal(config$KMT2A_ALU$gene_symbol, "KMT2A")
})

test_that("UBTF hotspot annotations use build-specific exon coordinates", {
  hotspots_path <- testthat::test_path("../../inst/extdata/hotspots.csv")
  cases <- data.frame(
    Gene = c("UBTF", "UBTF"),
    GenomicPosition = c(42288160L, 44210792L),
    Length = c(1L, 1L)
  )

  hg19 <- annotate_hotspots(cases[1L, ], db_path = hotspots_path,
                            genome_build = "hg19")
  hg38 <- annotate_hotspots(cases[2L, ], db_path = hotspots_path,
                            genome_build = "hg38")

  expect_true(hg19$Hotspot)
  expect_true(hg38$Hotspot)
  expect_equal(hg19$HotspotName, "UBTF_exon13")
  expect_equal(hg38$HotspotName, "UBTF_exon13")
})

test_that("poly-A run detection follows the configured minimum", {
  expect_equal(.detect_poly_a("CCAAAAAGG", min_run = 5L), 5L)
  expect_equal(.detect_poly_a("CCAAAAAGG", min_run = 6L), 0L)
})

test_that("LRU cache promotes hits and evicts the least-recently-used key", {
  cache <- .lru_cache(max_size = 2L)
  cache$set("a", 1L)
  cache$set("b", 2L)
  expect_equal(cache$get("a"), 1L)
  cache$set("c", 3L)

  expect_equal(cache$get("a"), 1L)
  expect_null(cache$get("b"))
  expect_equal(cache$get("c"), 3L)
  expect_equal(cache$size(), 2L)
})

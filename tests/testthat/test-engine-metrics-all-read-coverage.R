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

test_that("direct ITD coverage processes every support read across chunks", {
  skip_if_not_installed("pwalign")
  skip_if_not_installed("Biostrings")
  align_all <- get_fn(".align_reads_to_itd_template")
  set.seed(31)
  template <- paste(sample(c("A", "C", "G", "T"), 300L, replace = TRUE), collapse = "")
  starts <- c(1L, 61L, 121L, 181L, 241L)
  reads <- vapply(starts, function(s) substr(template, s, s + 59L), character(1L))
  reads_rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(reads)))
  # The short incidental fragment must not count as a qualifying read.
  got <- align_all(c(reads, reads_rc, "ACGTACGTAC"), template, chunk_size = 2L)

  expect_equal(got$reads_processed, 11L)
  expect_equal(got$reads_matched, 10L)
  expect_equal(got$percent, 100)
  expect_length(got$cov_vec, 300L)
  expect_equal(sum(got$cov_vec > 0L), 300L)
  rle_parts <- strsplit(got$rle_str, ",", fixed = TRUE)[[1L]]
  rle_widths <- vapply(rle_parts, function(x) as.integer(sub("^.*:", "", x)), integer(1L))
  expect_equal(sum(rle_widths), 300L)
})

test_that("direct ITD coverage leaves deleted reference bases uncovered", {
  skip_if_not_installed("pwalign")
  skip_if_not_installed("Biostrings")
  align_all <- get_fn(".align_reads_to_itd_template")
  set.seed(41)
  template <- paste(sample(c("A", "C", "G", "T"), 180L, replace = TRUE), collapse = "")
  # Remove 5 internal reference bases; the local alignment must not bridge
  # the deletion as though those reference positions had been sequenced.
  read <- paste0(substr(template, 1L, 70L), substr(template, 76L, 180L))
  got <- align_all(read, template, chunk_size = 1L, min_pid = 90, min_aligned = 20L)
  expect_equal(got$reads_matched, 1L)
  expect_equal(sum(got$cov_vec > 0L), 175L)
  expect_true(all(got$cov_vec[71L:75L] == 0L))
})


test_that("mate read sequences are extracted only for supporting QNAMEs", {
  skip_if_not_installed("GenomicAlignments")
  skip_if_not_installed("Biostrings")
  extract_mates <- get_fn(".extract_support_mate_coverage_rows")
  reads <- GenomicAlignments::GAlignments(
    seqnames = S4Vectors::Rle(factor(c("chr1", "chr1", "chr1"))),
    pos = c(10L, 100L, 200L),
    cigar = c("8M", "8M", "8M"),
    strand = S4Vectors::Rle(BiocGenerics::strand(c("+", "-", "+"))),
    names = c("fragA", "fragA", "fragB")
  )
  S4Vectors::mcols(reads)$seq <- Biostrings::DNAStringSet(c("ACGTACGT", "TTTTCCCC", "GGGGAAAA"))
  S4Vectors::mcols(reads)$flag <- c(99L, 147L, 99L)
  S4Vectors::mcols(reads)$cigar <- c("8M", "8M", "8M")

  got <- extract_mates(reads, "fragA")
  expect_equal(nrow(got), 2L)
  expect_setequal(got$read_seq, c("ACGTACGT", "TTTTCCCC"))
  expect_true(all(got$read_name == "fragA"))
})

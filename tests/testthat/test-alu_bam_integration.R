# Full BAM-integration test for ALU detection (ITEM 1)
#
# This test exercises the entire detect_alu() pipeline against a synthetic
# BAM containing a known ALU insertion. Unlike the unit tests in
# test-alu_helpers.R, which exercise individual helpers in isolation, this
# test validates:
#
#   - BAM I/O (.load_bam_data_streaming reads the synthetic BAM correctly)
#   - Soft-clip extraction (.get_softclips recovers the ALU-bearing clips)
#   - ALU matching (.match_clip_to_alu finds the ALU consensus)
#   - Candidate extraction (.extract_alu_candidates populates the data.frame)
#   - Breakpoint clustering (.cluster_breakpoints groups reads at one site)
#   - Wildtype counting (.count_wildtype_at_position sees the WT reads)
#   - VAF calculation (SupportingReads / (SupportingReads + WildtypeReads))
#   - TSD detection (.detect_alu_tsd reports "read_supported" when both
#     lead and trail clips are present, "unknown" otherwise)
#   - Pattern position storage (ITEM 2: pattern_start/pattern_end columns
#     are populated in the candidate data.frame)
#
# The test is skipped when Rsamtools / Biostrings / pwalign /
# GenomicAlignments are not installed, since the BAM build + ALU alignment
# require those packages. When they are installed, the test runs end-to-end
# and asserts numerical properties of the result (insertion site within a
# tolerance of the truth, VAF close to the design, supporting/WT counts
# close to the design).

library(testthat)

# ---------------------------------------------------------------------------
# Helper: build a synthetic BAM with a known ALU insertion
# ---------------------------------------------------------------------------
#
# Layout (1kb synthetic reference, "chrT"):
#
#   ref:  [1 ........................................................................... 1000]
#                                  [--- TSD: ref[a..b] ---]
#                                          |
#                              insertion site (centre of the TSD)
#
# STANDARD TSD MODEL (the donor duplicates the reference target site on
# both sides of the ALU element):
#
#   reference:  ... ref[a-1] | ref[a..b] (target site, single copy) | ref[b+1] ...
#   donor:      ... ref[a-1] | ref[a..b] (upstream TSD copy)        | ALU |
#                       ref[a..b] (downstream TSD copy) | ref[b+1] ...
#
# Equivalently:
#   donor = ref[1..b] + ALU + ref[a..N]
#
# where a = ins_pos - tsd_length/2 + 1, b = a + tsd_length - 1, and
# tsd_seq = ref[a..b].
#
# This is the model .detect_alu_tsd() expects: trail-clip reads align
# through the upstream TSD copy (their last aligned reference base is at
# position b, the END of the TSD in the reference), and lead-clip reads
# align from the downstream TSD copy (their first aligned reference base is
# at position a, the START of the TSD in the reference). The trail-lead
# breakpoint gap is therefore b - (a-1) = tsd_length, which is exactly the
# TSD length .detect_alu_tsd() reports as "read_supported".
#
# Read classes generated:
#   1. n_trail supporting reads spanning the UPSTREAM junction (TSD end -> ALU).
#      Aligned to ref[upstream..b], trail clip = ALU[1..alu_in_read].
#   2. n_lead supporting reads spanning the DOWNSTREAM junction (ALU -> TSD start).
#      Aligned to ref[a..downstream], lead clip = ALU[end-alu_in_read+1..end].
#   3. n_wt wildtype reads covering the insertion site with no clip.
#
# All reads are forward-strand, MAPQ=60, primary, no duplicates. The test
# controls the VAF precisely by choosing n_trail, n_lead, n_wt.
#
# Returns a list with:
#   bam_path    - path to the indexed BAM
#   ref_seq     - the 1000bp reference (character)
#   ins_pos     - insertion site (integer, = 500)
#   tsd_seq     - the 10bp TSD sequence (character)
#   tsd_a       - reference start of the TSD (integer)
#   tsd_b       - reference end of the TSD (integer)
#   alu_element - the 280bp ALU element (character)
#   gene_config - a minimal gene_config list consumable by detect_alu()
#   consensus_fa - path to a small ALU consensus FASTA

.build_synthetic_alu_bam <- function(
  n_trail = 8L,
  n_lead = 8L,
  n_wt = 24L,
  alu_length = 280L,
  tsd_length = 10L,
  ref_length = 1000L,
  ins_pos = 500L
) {
  set.seed(2026L)

  # --- Synthetic reference (deterministic for reproducibility) --------------
  bases <- c("A", "C", "G", "T")
  ref_seq <- paste(sample(bases, ref_length, replace = TRUE), collapse = "")

  # --- ALU element: 280bp inspired by AluSx --------------------------------
  # Repeats a 60bp core to reach the desired length. This is NOT a real ALU
  # but it's a stable target the .match_clip_to_alu() pairwise alignment can
  # hit at score >= 0.6 when the consensus FASTA contains the same core.
  alu_core <- "GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATC"
  reps <- ceiling(alu_length / nchar(alu_core))
  alu_element <- substr(
    paste(rep(alu_core, reps), collapse = ""),
    1L, alu_length
  )

  # --- TSD: the target site in the reference, centred on ins_pos -----------
  # a = ins_pos - tsd_length/2 + 1, b = a + tsd_length - 1.
  # For ins_pos=500, tsd_length=10: a=496, b=505. TSD = ref[496..505].
  tsd_a <- ins_pos - tsd_length %/% 2L + 1L
  tsd_b <- tsd_a + tsd_length - 1L
  tsd_seq <- substr(ref_seq, tsd_a, tsd_b)
  stopifnot(nchar(tsd_seq) == tsd_length)

  # --- Donor allele (standard TSD model) ------------------------------------
  # donor = ref[1..b] + ALU + ref[a..N]
  # The upstream TSD copy (ref[a..b]) is at donor positions a..b (same as ref).
  # The downstream TSD copy (ref[a..b]) is at donor positions b+1+ALU..b+ALU+tsd_length.
  # The donor's ref[a..N] portion starts right after the ALU, so ref[a..b] is
  # the downstream TSD copy, and ref[b+1..N] continues from there.
  donor_allele <- paste0(
    substr(ref_seq, 1L, tsd_b),          # ref[1..b] (includes upstream TSD copy at ref[a..b])
    alu_element,                          # ALU element
    substr(ref_seq, tsd_a, ref_length)    # ref[a..N] (downstream TSD copy at the start, then ref[b+1..])
  )

  # --- Build reads ----------------------------------------------------------
  # Trail-clip reads (upstream junction, donor TSD end -> ALU):
  #   Aligned to ref[read_start..tsd_b] (width = tsd_b - read_start + 1)
  #   Trail clip = ALU[1..alu_in_read] (length = alu_in_read)
  #   In donor coordinates: read covers donor[read_start .. tsd_b + alu_in_read]
  trail_reads <- list()
  alu_in_read <- 30L  # how many ALU bases appear in the trail clip
  for (i in seq_len(n_trail)) {
    # Vary the read start a bit so reads aren't identical.
    read_start <- tsd_b - (30L + i)  # ~475..466 for tsd_b=505
    aligned_width <- tsd_b - read_start + 1L  # ~31..40
    clip_len <- alu_in_read
    read_seq <- paste0(
      substr(ref_seq, read_start, tsd_b),  # ref portion (aligned, ends at tsd_b)
      substr(alu_element, 1L, alu_in_read)  # ALU start (in clip)
    )
    trail_reads[[i]] <- list(
      qname = sprintf("support_trail_%02d", i),
      flag  = 0L,
      pos   = read_start,
      cigar = paste0(aligned_width, "M", clip_len, "S"),
      seq   = read_seq
    )
  }

  # Lead-clip reads (downstream junction, donor ALU -> TSD start):
  #   Aligned to ref[tsd_a..read_end] (width = read_end - tsd_a + 1)
  #   Lead clip = ALU[end - alu_in_read + 1 .. end]
  #   In donor coordinates: read covers donor[alu_end_pos - alu_in_read + 1 .. alu_end_pos + aligned_width]
  #   where alu_end_pos = tsd_b + alu_length (the last ALU base in donor)
  lead_reads <- list()
  alu_in_read_lead <- 25L
  ref_in_read <- 35L
  for (i in seq_len(n_lead)) {
    # Small variation in aligned width so reads aren't identical.
    aligned_width <- ref_in_read + (i %% 3L) - 1L
    aligned_start <- tsd_a  # always start at the TSD start in reference
    alu_end_portion <- substr(
      alu_element,
      nchar(alu_element) - alu_in_read_lead + 1L,
      nchar(alu_element)
    )
    lead_clip <- alu_end_portion  # pure ALU end, no TSD (TSD is in the aligned portion)
    read_seq <- paste0(
      lead_clip,
      substr(ref_seq, aligned_start, aligned_start + aligned_width - 1L)
    )
    lead_reads[[i]] <- list(
      qname = sprintf("support_lead_%02d", i),
      flag  = 0L,
      pos   = aligned_start,
      cigar = paste0(nchar(lead_clip), "S", aligned_width, "M"),
      seq   = read_seq
    )
  }

  # Wildtype reads: cover the insertion site with no clip
  wt_reads <- list()
  for (i in seq_len(n_wt)) {
    read_start <- ins_pos - 35L + (i %% 5L)  # spread around the insertion site
    if (read_start < 1L) read_start <- 1L
    read_width <- 75L
    if (read_start + read_width - 1L > ref_length) {
      read_width <- ref_length - read_start + 1L
    }
    read_seq <- substr(ref_seq, read_start, read_start + read_width - 1L)
    wt_reads[[i]] <- list(
      qname = sprintf("wt_%02d", i),
      flag  = 0L,
      pos   = read_start,
      cigar = paste0(read_width, "M"),
      seq   = read_seq
    )
  }

  all_reads <- c(trail_reads, lead_reads, wt_reads)

  # --- Build SAM text -------------------------------------------------------
  sam_header <- c(
    "@HD\tVN:1.6\tSO:unsorted",
    sprintf("@SQ\tSN:chrT\tLN:%d", ref_length)
  )
  sam_body <- vapply(all_reads, function(r) {
    qual <- strrep("I", nchar(r$seq))
    sprintf("%s\t%d\tchrT\t%d\t60\t%s\t*\t0\t0\t%s\t%s",
            r$qname, r$flag, r$pos, r$cigar, r$seq, qual)
  }, character(1L))
  sam_text <- c(sam_header, sam_body)

  # --- Write SAM, convert to BAM, index -------------------------------------
  tmp_dir <- tempfile("talos_alu_integration_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  sam_path <- file.path(tmp_dir, "synthetic_alu.sam")
  bam_path <- file.path(tmp_dir, "synthetic_alu.bam")
  writeLines(sam_text, sam_path)

  # asBam requires Rsamtools; it writes the BAM next to the SAM.
  bam_obj <- Rsamtools::asBam(sam_path, filepath = tmp_dir,
                              overwrite = TRUE, indexDestination = TRUE)
  # asBam returns the BAM path; rename for predictability.
  if (!file.exists(bam_obj)) stop("asBam did not produce a BAM file.")
  if (bam_obj != bam_path) {
    file.rename(bam_obj, bam_path)
    # Re-index at the new path (indexDestination = TRUE indexed bam_obj).
    Rsamtools::indexBam(bam_path)
  }

  # --- Build the ALU consensus FASTA (small, just the AluSx core) -----------
  consensus_fa <- file.path(tmp_dir, "alu_consensus.fa")
  writeLines(c(
    ">AluSx",
    alu_core,
    ">AluY",
    alu_core,
    ">AluJb",
    alu_core
  ), consensus_fa)

  # --- Build a minimal gene_config consumable by detect_alu() --------------
  gene_config <- list(
    gene            = "TEST_GENE",
    chrom           = "chrT",
    strand          = "+",
    transcript      = NA_character_,
    genomic_start   = 1L,
    genomic_end     = ref_length,
    rnaseq          = ref_seq,
    genomic_ref_seq = ref_seq,
    full_cdna       = ref_seq,
    build           = "hg19",
    target_exons    = NULL,
    all_exons       = NULL,
    cds_offset      = 0L,
    bsgenome_obj    = NULL,
    gene_settings   = NULL,
    alu_settings    = NULL,
    graph_settings  = NULL
  )

  list(
    bam_path     = bam_path,
    ref_seq      = ref_seq,
    ins_pos      = ins_pos,
    tsd_seq      = tsd_seq,
    tsd_a        = tsd_a,
    tsd_b        = tsd_b,
    alu_element  = alu_element,
    gene_config  = gene_config,
    consensus_fa = consensus_fa,
    tmp_dir      = tmp_dir,
    n_trail      = n_trail,
    n_lead       = n_lead,
    n_wt         = n_wt
  )
}

# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------
describe("ITEM 1: detect_alu() end-to-end on synthetic BAM", {
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("Biostrings")
  skip_if_not_installed("pwalign")
  skip_if_not_installed("GenomicAlignments")
  skip_if_not_installed("GenomicRanges")

  # Helper to look up detect_alu() (exported) and the internal helpers we
  # verify along the way.
  get_fn <- function(name) {
    if (exists(name, envir = globalenv(), inherits = FALSE)) {
      get(name, envir = globalenv())
    } else if (requireNamespace("TALOS", quietly = TRUE)) {
      getFromNamespace(name, "TALOS")
    } else {
      stop("Function ", name, " not found. Load the package with devtools::load_all().")
    }
  }

  it("detects an ALU insertion at the designed position with correct VAF", {
    setup <- .build_synthetic_alu_bam(
      n_trail = 8L, n_lead = 8L, n_wt = 24L, alu_length = 280L
    )
    on.exit(unlink(setup$tmp_dir, recursive = TRUE), add = TRUE)

    result <- detect_alu(
      bam_path        = setup$bam_path,
      gene_config     = setup$gene_config,
      consensus_fa    = setup$consensus_fa,
      min_support     = 3L,
      min_alu_score   = 0.55,    # the synthetic core is short; lower the bar
      min_clip_len    = 25L,
      min_mapq        = 20L,
      cluster_tolerance = 15L,
      vaf_threshold   = 0.05,
      min_tsd         = 0L,
      do_annotate_hotspots = FALSE,  # no hotspot DB in the synthetic test
      output_prefix   = NULL,        # skip TSV / HTML writing
      html_report     = FALSE,
      verbose         = FALSE
    )

    # We expect exactly one insertion call.
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 1L,
                 info = "expected exactly 1 ALU insertion call")

    # The insertion site should be at (or very close to) the designed 500.
    # Tolerance of cluster_tolerance (15bp) accounts for the trail/lead
    # breakpoint scatter built into the synthetic reads.
    expect_true(abs(result$InsertionSite - setup$ins_pos) <= 15L,
                info = sprintf("insertion site %d not within 15bp of %d",
                               result$InsertionSite, setup$ins_pos))

    # Supporting reads: should be the deduplicated count of trail + lead
    # reads. We built n_trail + n_lead = 16 unique qnames, no duplicates.
    expect_gte(result$SupportingReads, 8L)  # at least half (one side) should pass
    expect_lte(result$SupportingReads, 16L)

    # Wildtype reads: there are n_wt = 24 WT reads covering the insertion
    # site. .count_wildtype_at_position() deduplicates by qname and excludes
    # supporting reads, so the count should be in [20, 24].
    expect_gte(result$WildtypeReads, 18L)
    expect_lte(result$WildtypeReads, 24L)

    # VAF = SupportingReads / DepthAtBreakpoint.
    expected_depth <- result$SupportingReads + result$WildtypeReads
    expect_equal(result$DepthAtBreakpoint, as.integer(expected_depth))
    expected_vaf <- result$SupportingReads / expected_depth
    expect_equal(result$AlleleFrequency, round(expected_vaf, 4L))

    # ALU subtype should be one of the consensus entries.
    expect_in(result$ALUSubtype, c("AluSx", "AluY", "AluJb"))

    # Score should meet the min_alu_score threshold.
    expect_gte(result$ALU_RawScore, 0.55)
  })

  it("reports TSD_Status as 'read_supported' when both lead and trail clips are present", {
    # With n_trail > 0 AND n_lead > 0, .detect_alu_tsd() should be able to
    # compute a TSD length from the trail-lead breakpoint gap. The status
    # is "read_supported" if the gap is in [min_tsd, max_tsd] = [5, 20].
    setup <- .build_synthetic_alu_bam(
      n_trail = 10L, n_lead = 10L, n_wt = 20L
    )
    on.exit(unlink(setup$tmp_dir, recursive = TRUE), add = TRUE)

    result <- detect_alu(
      bam_path        = setup$bam_path,
      gene_config     = setup$gene_config,
      consensus_fa    = setup$consensus_fa,
      min_support     = 3L,
      min_alu_score   = 0.55,
      min_clip_len    = 25L,
      min_mapq        = 20L,
      cluster_tolerance = 15L,
      vaf_threshold   = 0.05,
      min_tsd         = 0L,
      do_annotate_hotspots = FALSE,
      output_prefix   = NULL,
      html_report     = FALSE,
      verbose         = FALSE
    )

    expect_equal(nrow(result), 1L)
    # With the standard TSD model, trail_bp = tsd_b and lead_bp = tsd_a - 1,
    # so TSD length = tsd_b - (tsd_a - 1) = tsd_length = 10. The gap is in
    # [5, 20], so .detect_alu_tsd() should report "read_supported".
    expect_equal(result$TSD_Status, "read_supported",
                 info = sprintf("expected read_supported, got %s (TSD_Length=%s)",
                                result$TSD_Status,
                                ifelse(is.na(result$TSD_Length), "NA",
                                       as.character(result$TSD_Length))))
    expect_equal(result$TSD_Length, 10L)
    expect_equal(result$TSD_Sequence, setup$tsd_seq)
  })

  it("reports TSD_Status as 'unknown' when only trail clips are present", {
    # Only trail-side reads -> only one side of the breakpoint is represented.
    # .detect_alu_tsd() returns "unknown" in this case (Finding 1 fix).
    setup <- .build_synthetic_alu_bam(
      n_trail = 10L, n_lead = 0L, n_wt = 20L
    )
    on.exit(unlink(setup$tmp_dir, recursive = TRUE), add = TRUE)

    result <- detect_alu(
      bam_path        = setup$bam_path,
      gene_config     = setup$gene_config,
      consensus_fa    = setup$consensus_fa,
      min_support     = 3L,
      min_alu_score   = 0.55,
      min_clip_len    = 25L,
      min_mapq        = 20L,
      cluster_tolerance = 15L,
      vaf_threshold   = 0.05,
      min_tsd         = 0L,
      do_annotate_hotspots = FALSE,
      output_prefix   = NULL,
      html_report     = FALSE,
      verbose         = FALSE
    )

    # May or may not pass the min_support filter with only trail reads.
    if (nrow(result) >= 1L) {
      expect_equal(result$TSD_Status, "unknown")
    }
  })

  it("returns an empty data.frame when no ALU-bearing reads exist", {
    # Only WT reads -> no ALU-positive candidates -> empty result.
    setup <- .build_synthetic_alu_bam(
      n_trail = 0L, n_lead = 0L, n_wt = 30L
    )
    on.exit(unlink(setup$tmp_dir, recursive = TRUE), add = TRUE)

    result <- detect_alu(
      bam_path        = setup$bam_path,
      gene_config     = setup$gene_config,
      consensus_fa    = setup$consensus_fa,
      min_support     = 3L,
      min_alu_score   = 0.55,
      do_annotate_hotspots = FALSE,
      output_prefix   = NULL,
      html_report     = FALSE,
      verbose         = FALSE
    )

    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 0L)
  })

  it("candidate data.frame includes the ITEM 2 pattern_start/pattern_end columns", {
    # This is a pipeline-level check that the ITEM 2 columns survive the
    # .extract_alu_candidates() -> .summarise_alu_cluster() journey. The
    # columns are internal to the candidate data.frame (not in the final
    # output) so we call .extract_alu_candidates() directly.
    skip_if_not_installed("Biostrings")
    skip_if_not_installed("pwalign")

    setup <- .build_synthetic_alu_bam(
      n_trail = 6L, n_lead = 6L, n_wt = 10L
    )
    on.exit(unlink(setup$tmp_dir, recursive = TRUE), add = TRUE)

    # Look up internal helpers (engine_bam.R + engine_alu.R)
    load_bam_fn  <- get_fn(".load_bam_data_streaming")
    extract_fn   <- get_fn(".extract_alu_candidates")

    # Load BAM reads the same way detect_alu() does.
    bam_data <- load_bam_fn(
      setup$bam_path,
      setup$gene_config,
      compute_pairs = FALSE,
      max_reads = 200000L,
      verbose = FALSE
    )
    expect_true(length(bam_data$reads) > 0L,
                info = "no reads loaded from synthetic BAM")

    alu_seqs <- Biostrings::readDNAStringSet(setup$consensus_fa)

    cands <- extract_fn(
      reads          = bam_data$reads,
      alu_seqs       = alu_seqs,
      genomic_start  = setup$gene_config$genomic_start,
      min_clip_len   = 25L,
      min_alu_score  = 0.55,
      verbose        = FALSE
    )

    # We should have at least one trail and one lead candidate.
    expect_gte(nrow(cands), 2L)
    expect_in(c("trail", "lead"), unique(cands$clip_side))

    # ITEM 2: pattern_start / pattern_end / boundary_bases columns present
    expect_in(c("pattern_start", "pattern_end", "boundary_bases"),
              names(cands))
    # All candidate rows should have non-NA pattern positions.
    expect_false(any(is.na(cands$pattern_start)))
    expect_false(any(is.na(cands$pattern_end)))
    # pattern_start <= pattern_end (within the clip)
    expect_true(all(cands$pattern_start <= cands$pattern_end))
    # pattern_end <= nchar(clip_seq) (the alignment can't exceed the clip)
    expect_true(all(cands$pattern_end <= nchar(cands$clip_seq)))
    expect_true(all(cands$pattern_start >= 1L))

    # For trail clips, the clip = ALU[1..alu_in_read] (pure ALU, no TSD
    # prefix in this standard-model construction), so pattern_start should
    # be 1 and boundary_bases should be "" (empty prefix).
    # For lead clips, the clip = ALU[end-alu_in_read+1..end] (pure ALU, no
    # TSD suffix), so pattern_end should equal nchar(clip_seq) and
    # boundary_bases should be "" (empty suffix).
    trail_rows <- cands[cands$clip_side == "trail", , drop = FALSE]
    lead_rows  <- cands[cands$clip_side == "lead",  , drop = FALSE]
    if (nrow(trail_rows) > 0L) {
      expect_true(all(trail_rows$pattern_start == 1L),
                  info = "trail clip is pure ALU; pattern_start should be 1")
      expect_true(all(trail_rows$boundary_bases == "" |
                      is.na(trail_rows$boundary_bases)),
                  info = "trail clip has no TSD prefix; boundary should be empty")
    }
    if (nrow(lead_rows) > 0L) {
      expect_true(all(lead_rows$pattern_end == nchar(lead_rows$clip_seq)),
                  info = "lead clip is pure ALU; pattern_end should equal clip length")
      expect_true(all(lead_rows$boundary_bases == "" |
                      is.na(lead_rows$boundary_bases)),
                  info = "lead clip has no TSD suffix; boundary should be empty")
    }

    # The boundary_bases column should be populated (not NA) for every row.
    expect_false(any(is.na(cands$boundary_bases)))
  })

  it("respects min_support: filters out insertions with too few reads", {
    # Build a setup with only 2 supporting reads (below min_support=3).
    setup <- .build_synthetic_alu_bam(
      n_trail = 1L, n_lead = 1L, n_wt = 20L
    )
    on.exit(unlink(setup$tmp_dir, recursive = TRUE), add = TRUE)

    result <- detect_alu(
      bam_path        = setup$bam_path,
      gene_config     = setup$gene_config,
      consensus_fa    = setup$consensus_fa,
      min_support     = 3L,
      min_alu_score   = 0.55,
      do_annotate_hotspots = FALSE,
      output_prefix   = NULL,
      html_report     = FALSE,
      verbose         = FALSE
    )

    # Only 2 supporting reads -> below min_support=3 -> no calls.
    expect_equal(nrow(result), 0L)
  })

  it("VAF scales correctly with WT read count", {
    # Run two scenarios with different WT counts and confirm VAF changes.
    setup_a <- .build_synthetic_alu_bam(
      n_trail = 8L, n_lead = 8L, n_wt = 16L  # ~50% VAF
    )
    setup_b <- .build_synthetic_alu_bam(
      n_trail = 8L, n_lead = 8L, n_wt = 64L  # ~20% VAF
    )
    on.exit({
      unlink(setup_a$tmp_dir, recursive = TRUE)
      unlink(setup_b$tmp_dir, recursive = TRUE)
    }, add = TRUE)

    common_args <- list(
      min_alu_score   = 0.55,
      min_clip_len    = 25L,
      min_mapq        = 20L,
      cluster_tolerance = 15L,
      vaf_threshold   = 0.05,
      min_tsd         = 0L,
      do_annotate_hotspots = FALSE,
      output_prefix   = NULL,
      html_report     = FALSE,
      verbose         = FALSE
    )

    res_a <- do.call(detect_alu, c(
      list(bam_path = setup_a$bam_path, gene_config = setup_a$gene_config,
           consensus_fa = setup_a$consensus_fa, min_support = 3L),
      common_args
    ))
    res_b <- do.call(detect_alu, c(
      list(bam_path = setup_b$bam_path, gene_config = setup_b$gene_config,
           consensus_fa = setup_b$consensus_fa, min_support = 3L),
      common_args
    ))

    # Both should produce one call (the insertion is the same design).
    expect_equal(nrow(res_a), 1L)
    expect_equal(nrow(res_b), 1L)

    # The VAF in setup_b (4x more WT) should be lower than in setup_a.
    expect_lt(res_b$AlleleFrequency, res_a$AlleleFrequency,
              info = sprintf("VAF should drop with more WT reads: a=%g, b=%g",
                             res_a$AlleleFrequency, res_b$AlleleFrequency))
  })
})

# ============================================================================
# TALOS – Variant metrics, VAF computation, and filter application
# ============================================================================

# ---------------------------------------------------------------------------
# WT read count spanning the breakpoint (buffer bp on both sides)
# ---------------------------------------------------------------------------
.compute_left_right_coverage <- function(wildtype_info, len_specific_bp,
                                          support_qnames, buffer = 10L) {
  strict_idx <- which(
    wildtype_info$is_primary &
    start(wildtype_info$gr) <= (len_specific_bp - buffer) &
    end(wildtype_info$gr)   >= (len_specific_bp + buffer)
  )
  strict_qnames <- unique(wildtype_info$qnames[strict_idx])
  pure_wt       <- setdiff(strict_qnames, support_qnames)
  wt_count      <- length(pure_wt)
  list(left = wt_count, right = wt_count, avg = wt_count)
}


# ---------------------------------------------------------------------------
# Size-bias corrected VAF and depth
# ---------------------------------------------------------------------------
.calculate_vaf_and_depth <- function(raw_support, wildtype_info, len_specific_bp,
                                      support_qnames, best_len,
                                      nominal_read_len, max_correction,
                                      buffer = 10L) {
  lr        <- .compute_left_right_coverage(wildtype_info, len_specific_bp,
                                            support_qnames, buffer)
  strict_wt <- lr$avg
  len_used  <- if (is.na(best_len)) 0L else best_len
  size_scalar       <- min(1 + (0.25 * len_used / nominal_read_len), max_correction)
  corrected_support <- round(raw_support * size_scalar)
  depth_at_bp       <- corrected_support + strict_wt
  corrected_af      <- if (depth_at_bp == 0L) 0 else corrected_support / depth_at_bp
  list(
    corrected_support = corrected_support, depth_at_bp = depth_at_bp,
    raw_af = corrected_af, raw_depth = raw_support + strict_wt,
    left_wt = strict_wt, right_wt = strict_wt, avg_wt = strict_wt
  )
}


# ---------------------------------------------------------------------------
# Compute all variant-level metrics for one breakpoint cluster / length group
# ---------------------------------------------------------------------------
.compute_variant_metrics <- function(cluster_bps, best_len, support_rows,
                                      genomic_start, ref_dna, gene_config,
                                      all_reads_cov, all_pairs, wildtype_info,
                                      ptd_mode, min_support, min_wt_reads,
                                      nominal_read_len, max_correction,
                                      vaf_threshold,
                                      do_alignment_score, do_support_bases,
                                      do_consistency, do_itd_coverage,
                                      do_coverage_drop, do_microhomology,
                                      do_repeat_entropy, do_discordant_ratio,
                                      do_detect_orientation, do_hgvs,
                                      max_pairwise_alignments,
                                      length_ext = NA_integer_,
                                      pre_assembled_seq = NA_character_,
                                      debug = FALSE, verbose = FALSE) {

  if (verbose) message(
    "[TALOS] Computing metrics for cluster at ",
    paste(range(cluster_bps), collapse = "-"), " length=", best_len
  )

  # ---- Initialise all metric variables to NA ----
  consistency_score    <- NA_real_
  itd_coverage_percent <- NA_real_
  itd_coverage_rle     <- NA_character_
  ref_match_observed   <- NA_real_
  ref_match_total      <- NA_real_
  alignment_score_total<- NA_real_
  total_support_bases  <- NA_integer_
  coverage_drop        <- NA_real_
  local_coverage        <- NA_real_
  median_microhomology <- NA_real_
  repeat_entropy       <- NA_real_
  discordant_ratio     <- NA_real_
  length_pe            <- NA_real_
  length_pe_nspanning  <- NA_integer_
  pe_softclip_support  <- NA_integer_
  pe_softclip_event_pairs <- NA_integer_
  pe_softclip_long_pairs  <- NA_integer_
  pe_orientation_fr    <- NA_integer_
  pe_orientation_rf    <- NA_integer_
  pe_orientation_ff    <- NA_integer_
  pe_orientation_rr    <- NA_integer_
  pe_orientation_other <- NA_integer_
  pe_orientation_dominant <- NA_character_
  orientation          <- NA_character_
  hgvs <- list(c_notation = NA_character_, p_notation = NA_character_)

  len_specific_bp <- round(median(cluster_bps, na.rm = TRUE))
  if (is.na(len_specific_bp) || is.infinite(len_specific_bp)) return(NULL)
  support_qnames  <- unique(support_rows$read_name)
  raw_support     <- length(support_qnames)

  # ---- Soft-clip side counts ----
  left_sc_count  <- 0L; right_sc_count <- 0L; both_sc_count <- 0L
  if (nrow(support_rows) > 0L && !is.null(support_rows$cigar)) {
    for (cig in support_rows$cigar) {
      if (is.na(cig)) next
      has_left  <- grepl("^\\d+S", cig, perl = TRUE)
      has_right <- grepl("\\d+S$",  cig, perl = TRUE)
      if (has_left && has_right) {
        both_sc_count  <- both_sc_count + 1L
        left_sc_count  <- left_sc_count  + 1L
        right_sc_count <- right_sc_count + 1L
      } else if (has_left)  { left_sc_count  <- left_sc_count  + 1L
      } else if (has_right) { right_sc_count <- right_sc_count + 1L }
    }
  }

  # ---- ITD sequence extraction ----
  itd_seq <- NA_character_; imputed <- TRUE; sequence_partial <- FALSE
  observed_seq <- NA_character_; observed_len <- 0L

  valid_seqs <- support_rows$itd_seq[!is.na(support_rows$itd_seq) &
                                      nchar(support_rows$itd_seq) > 0L]
  if (length(valid_seqs) > 0L) {
    itd_seq      <- names(sort(table(valid_seqs), decreasing = TRUE))[1L]
    imputed      <- FALSE
    observed_seq <- itd_seq
    observed_len <- nchar(itd_seq)
  }

  if (is.na(itd_seq)) {
    extracted_seqs <- vector("character", 2L * nrow(support_rows))
    seq_idx <- 1L
    for (i in seq_len(nrow(support_rows))) {
      sc <- .get_softclips(support_rows$cigar[i], support_rows$read_seq[i])
      if (!is.na(sc$lead))  { extracted_seqs[seq_idx] <- sc$lead;  seq_idx <- seq_idx + 1L }
      if (!is.na(sc$trail)) { extracted_seqs[seq_idx] <- sc$trail; seq_idx <- seq_idx + 1L }
    }
    extracted_seqs <- extracted_seqs[seq_len(seq_idx - 1L)]
    if (length(extracted_seqs) > 0L) {
      itd_seq      <- names(sort(table(extracted_seqs), decreasing = TRUE))[1L]
      imputed      <- FALSE
      observed_seq <- itd_seq
      observed_len <- nchar(itd_seq)
    }
  }

  if (is.na(itd_seq) && !is.na(pre_assembled_seq) && nchar(pre_assembled_seq) > 0L) {
    itd_seq      <- pre_assembled_seq
    imputed      <- FALSE
    observed_seq <- itd_seq
    observed_len <- nchar(itd_seq)
  }

  if (!is.na(itd_seq) && !is.na(best_len) && best_len > 0L &&
      nchar(itd_seq) > best_len) {
    itd_seq      <- substr(itd_seq, 1L, best_len)
    if (!is.na(observed_seq) && nchar(observed_seq) > best_len)
      observed_seq <- substr(observed_seq, 1L, best_len)
    observed_len <- min(observed_len, best_len)
  }

  if (!is.na(itd_seq) && !is.na(best_len) && best_len > 0L &&
      nchar(itd_seq) < best_len && !ptd_mode) {
    partial_obs_len <- nchar(itd_seq)
    ref_local_start <- len_specific_bp - genomic_start + partial_obs_len + 1L
    ref_local_end   <- ref_local_start + (best_len - partial_obs_len) - 1L
    if (ref_local_end <= nchar(ref_dna)) {
      itd_seq          <- paste0(itd_seq, substr(ref_dna, ref_local_start, ref_local_end))
      sequence_partial <- TRUE; imputed <- TRUE
      observed_len     <- partial_obs_len
    }
  }

  if (is.na(itd_seq) && !is.na(best_len) && best_len > 0L && !ptd_mode) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end   <- local_start + best_len - 1L
    if (!is.na(local_end) && local_end <= nchar(ref_dna)) {
      itd_seq          <- substr(ref_dna, local_start, local_end)
      imputed          <- TRUE; sequence_partial <- FALSE
      observed_seq     <- NA_character_; observed_len <- 0L
    }
  }


  # ---- Consistency / ITD coverage ----
  if ((do_consistency || do_itd_coverage) &&
      !is.na(itd_seq) && nchar(itd_seq) > 0L && nrow(support_rows) > 0L) {
    valid_mask <- !is.na(support_rows$read_seq) & nchar(support_rows$read_seq) > 0L
    v_seqs     <- support_rows$read_seq[valid_mask]
    if (length(v_seqs) > 0L) {
      if (length(v_seqs) > max_pairwise_alignments) {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
          .GlobalEnv$.Random.seed else NULL
        set.seed(42L)
        v_seqs <- v_seqs[sample(length(v_seqs), max_pairwise_alignments)]
        if (!is.null(old_seed)) .GlobalEnv$.Random.seed <- old_seed
        else rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      }
      if (.has_pwalign && .has_biostrings) {
        aln  <- pwalign::pairwiseAlignment(
          pattern = Biostrings::DNAStringSet(v_seqs),
          subject = Biostrings::DNAString(itd_seq), type = "local"
        )
        pids <- pwalign::pid(aln)
        if (do_consistency)
          consistency_score <- sum(pids >= 90, na.rm = TRUE) / length(v_seqs) * 100
        if (do_itd_coverage) {
          good_idx <- which(pids >= 90)
          if (length(good_idx) > 0L) {
            subj_ranges <- pwalign::subject(aln[good_idx])
            cov_vec <- integer(nchar(itd_seq))
            for (ki in seq_along(good_idx))
              cov_vec[start(subj_ranges)[ki]:end(subj_ranges)[ki]] <-
                cov_vec[start(subj_ranges)[ki]:end(subj_ranges)[ki]] + 1L
            itd_coverage_percent <- sum(cov_vec > 0L) / length(cov_vec) * 100
            rle_res <- rle(cov_vec)
            itd_coverage_rle <- paste(paste0(rle_res$values, ":", rle_res$lengths),
                                       collapse = ",")
          } else {
            itd_coverage_percent <- 0
            itd_coverage_rle     <- paste0("0:", nchar(itd_seq))
          }
        }
      } else {
        if (do_consistency)
          consistency_score <- sum(grepl(itd_seq, v_seqs, fixed = TRUE)) /
                               length(v_seqs) * 100
        if (do_itd_coverage) {
          itd_coverage_percent <- if (any(grepl(itd_seq, v_seqs, fixed = TRUE))) 100 else 0
          itd_coverage_rle     <- paste0(as.integer(itd_coverage_percent > 0),
                                          ":", nchar(itd_seq))
        }
      }
    }
  }

  # ---- Alignment scores ----
  if (do_alignment_score && !is.na(itd_seq) && !is.na(best_len) && best_len > 0L) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end   <- local_start + best_len - 1L
    if (!is.na(local_end) && local_end <= nchar(ref_dna)) {
      ref_segment  <- substr(ref_dna, local_start, local_end)
      min_len_tot  <- min(nchar(itd_seq), nchar(ref_segment))
      if (min_len_tot > 0L) {
        alignment_score_total <-
          if (itd_seq == ref_segment) 1.0
          else sum(charToRaw(substr(itd_seq, 1L, min_len_tot)) ==
                   charToRaw(substr(ref_segment, 1L, min_len_tot))) / min_len_tot
        ref_match_total <- round(alignment_score_total * 100, 1)
        if (!is.na(observed_seq) && observed_len > 0L) {
          obs_len <- min(observed_len, min_len_tot)
          if (obs_len > 0L) {
            obs_sub <- substr(observed_seq, 1L, obs_len)
            ref_sub <- substr(ref_segment, 1L, obs_len)
            alignment_score_obs <- if (obs_sub == ref_sub) 1.0
              else sum(charToRaw(obs_sub) == charToRaw(ref_sub)) / obs_len
            ref_match_observed <- round(alignment_score_obs * 100, 1)
          }
        }
      }
    }
  }

  # ---- Support bases (total soft-clip bases) ----
  if (do_support_bases && nrow(support_rows) > 0L) {
    cigars  <- support_rows$cigar; cigars[is.na(cigars)] <- ""
    lead_S  <- as.integer(regmatches(cigars, regexec("^(\\d+)S", cigars)) |>
                          sapply(function(m) if (length(m) > 1L) m[2L] else "0"))
    trail_S <- as.integer(regmatches(cigars, regexec("(\\d+)S$", cigars)) |>
                          sapply(function(m) if (length(m) > 1L) m[2L] else "0"))
    total_support_bases <- sum(lead_S, trail_S, na.rm = TRUE)
  }

  # ---- Minimum support check ----
  if (raw_support < min_support) return(NULL)

  # ---- VAF / depth ----
  vaf_metrics       <- .calculate_vaf_and_depth(
    raw_support, wildtype_info, len_specific_bp, support_qnames,
    best_len, nominal_read_len, max_correction, buffer = 10L
  )
  corrected_support <- vaf_metrics$corrected_support
  depth_at_bp       <- vaf_metrics$depth_at_bp
  raw_af            <- vaf_metrics$raw_af
  wildtype_reads    <- round(vaf_metrics$avg_wt)

  if (wildtype_reads < min_wt_reads) return(NULL)
  if (raw_af < vaf_threshold)        return(NULL)

  # ---- WT-relative soft-clip fractions ----
  left_sc_pct_wt  <- if (wildtype_reads > 0L) (left_sc_count  / wildtype_reads) * 100 else NA_real_
  right_sc_pct_wt <- if (wildtype_reads > 0L) (right_sc_count / wildtype_reads) * 100 else NA_real_

  # ---- Read-level diagnostic metrics ----
  strand_bias       <- if (nrow(support_rows) > 0L)
    round(mean(support_rows$is_reverse, na.rm = TRUE), 4) else NA_real_
  mean_support_mapq <- if (nrow(support_rows) > 0L)
    round(mean(support_rows$mapq, na.rm = TRUE), 1) else NA_real_
  bp_spread         <- if (nrow(support_rows) > 1L)
    as.integer(max(support_rows$breakpoint) - min(support_rows$breakpoint)) else 0L
  softclip_frac     <- if (nrow(support_rows) > 0L)
    round(mean(grepl("S", support_rows$cigar, fixed = TRUE), na.rm = TRUE), 4) else NA_real_
  unique_bps        <- length(unique(support_rows$breakpoint))
  left_sc_pct_support  <- if (corrected_support > 0L)
    (left_sc_count  / corrected_support) * 100 else NA_real_
  right_sc_pct_support <- if (corrected_support > 0L)
    (right_sc_count / corrected_support) * 100 else NA_real_

  # ---- Coverage drop ----
  if (do_coverage_drop && !is.na(best_len) && best_len > 0L) {
    cov_ok <- FALSE
    if (!is.null(all_reads_cov)) {
      coverage_drop <- compute_coverage_drop(all_reads_cov, gene_config$chrom,
                                             len_specific_bp, flank = 200L,
                                             verbose = debug)
      if (!is.na(coverage_drop)) cov_ok <- TRUE
    }
    if (!cov_ok && !is.null(wildtype_info$gr) && length(wildtype_info$gr) > 0L) {
      if (debug) message("Coverage drop: recomputing from wildtype reads")
      reg_gr <- GenomicRanges::GRanges(
        seqnames = gene_config$chrom,
        ranges   = IRanges::IRanges(start = gene_config$genomic_start,
                                    end   = gene_config$genomic_end)
      )
      ov <- GenomicRanges::findOverlaps(wildtype_info$gr, reg_gr)
      if (length(ov) > 0L) {
        reads_sub    <- wildtype_info$gr[unique(S4Vectors::queryHits(ov))]
        cov_fallback <- GenomicAlignments::coverage(reads_sub)[[gene_config$chrom]]
        if (!is.null(cov_fallback))
          coverage_drop <- compute_coverage_drop(cov_fallback, gene_config$chrom,
                                                 len_specific_bp, flank = 200L,
                                                 verbose = debug)
      }
    }
    if (debug && is.na(coverage_drop)) message("Coverage drop: still NA after fallback")
  }

  # ---- Absolute local coverage (no fallback -- see compute_local_coverage) ----
  if (!is.null(all_reads_cov)) {
    local_coverage <- compute_local_coverage(all_reads_cov, gene_config$chrom,
                                             len_specific_bp, buffer = 10L,
                                             verbose = debug)
  }

  # ---- Microhomology ----
  if (do_microhomology)
    median_microhomology <- compute_microhomology(support_rows, ref_dna,
                                                  len_specific_bp, genomic_start,
                                                  debug = debug)

  # ---- Repeat entropy ----
  if (do_repeat_entropy)
    repeat_entropy <- compute_repeat_entropy(ref_dna, len_specific_bp, genomic_start)

  # ---- Discordant ratio ----
  if (do_discordant_ratio && !is.null(all_pairs))
    discordant_ratio <- compute_discordant_ratio(all_pairs, len_specific_bp)

  # ---- Paired-end ITD length estimate (also for PTD mode) ----
  if (!is.null(all_pairs)) {
    pe_est <- compute_pe_itd_length(
      all_pairs, len_specific_bp,
      nominal_read_len = nominal_read_len
    )
    length_pe           <- pe_est$length_pe
    length_pe_nspanning <- as.integer(pe_est$n_spanning)
    if (verbose && !is.na(length_pe))
      message("[TALOS] PE insert-size ITD length estimate: ", length_pe,
              " bp (", length_pe_nspanning, " spanning pairs)")

    pe_sc <- compute_pe_softclip_metrics(
      all_pairs = all_pairs,
      breakpoint = len_specific_bp,
      support_rows = support_rows,
      best_len = best_len
    )
    pe_softclip_support <- pe_sc$pe_softclip_support
    pe_softclip_event_pairs <- pe_sc$pe_softclip_event_pairs
    pe_softclip_long_pairs <- pe_sc$pe_softclip_long_pairs
    pe_orientation_fr <- pe_sc$pe_orientation_fr
    pe_orientation_rf <- pe_sc$pe_orientation_rf
    pe_orientation_ff <- pe_sc$pe_orientation_ff
    pe_orientation_rr <- pe_sc$pe_orientation_rr
    pe_orientation_other <- pe_sc$pe_orientation_other
    pe_orientation_dominant <- pe_sc$pe_orientation_dominant
  }

  # ---- Orientation ----
  if (do_detect_orientation && !is.na(itd_seq) && !is.na(best_len) && best_len > 0L) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end   <- local_start + best_len - 1L
    if (!is.na(local_end) && local_end <= nchar(ref_dna))
      orientation <- detect_orientation(itd_seq,
                                        substr(ref_dna, local_start, local_end))
  }

  # ---- HGVS ----
  is_exonic <- FALSE
  if (!is.null(gene_config$exons)) {
    bp_gr <- GenomicRanges::GRanges(
      seqnames = gene_config$chrom,
      ranges   = IRanges::IRanges(start = len_specific_bp, end = len_specific_bp)
    )
    is_exonic <- length(GenomicRanges::findOverlaps(bp_gr, gene_config$exons,
                                                    type = "within")) > 0L
  }
  if (do_hgvs) {
    dup_len_for_hgvs <- if (!is.na(itd_seq) && nchar(itd_seq) > 0L) nchar(itd_seq)
                        else if (!is.na(best_len))                   best_len
                        else                                         NA_integer_
    if (!is.null(gene_config$all_exons) || !is.null(gene_config$exons)) {
      hgvs <- compute_hgvs_annotations(gene_config, len_specific_bp,
                                       itd_seq, dup_len_for_hgvs, debug)
      if (!is_exonic && !is.na(hgvs$p_notation)) {
        hgvs$p_notation <- "intronic"
        if (is.na(hgvs$c_notation)) hgvs$c_notation <- "c.?+?"
      }
    }
  }

  # ---- Artifact flag: fully-imputed sequence that is self-fulfilling
  # (imputed from reference then compared against the same reference).
  # Restored from old code: only flag when fully imputed (not partial
  # read+ref) and sequence identity is ~100%.
  artifact_suspect <- isTRUE(imputed) && !isTRUE(sequence_partial) &&
                      !is.na(ref_match_total) && ref_match_total >= 99.0

  if (debug && artifact_suspect)
    message(sprintf("  [ARTIFACT] bp=%d fully imputed + 100%% ref match", len_specific_bp))

  if (verbose) message("[TALOS] Metrics computed: support=", corrected_support,
                       ", AF=", round(raw_af, 4), ", length=", best_len)

  list(
    Sample = wildtype_info$sample_name, Gene = gene_config$gene,
    Genome = if (is.null(gene_config$build)) "custom" else gene_config$build,
    GenomicPosition = len_specific_bp, ITD_Sequence = itd_seq, Length = best_len,
    SupportingReads = corrected_support, WildtypeReads = wildtype_reads,
    DepthAtBreakpoint = depth_at_bp, AlleleFrequency = raw_af,
    HGVS_cDNA = hgvs$c_notation, HGVS_Protein = hgvs$p_notation,
    StrandBias = strand_bias, MeanSupportMAPQ = mean_support_mapq,
    BreakpointSpread = bp_spread, SoftclipFraction = softclip_frac,
    UniqueBreakpoints = unique_bps, CoverageDrop = coverage_drop,
    LocalCoverage = local_coverage,
    MedianMicrohomology = median_microhomology, DiscordantRatio = discordant_ratio,
    RepeatEntropy = repeat_entropy,
    SequenceImputed = imputed, SequencePartial = sequence_partial,
    SequenceSource = if (!imputed) "observed"
                     else if (sequence_partial) "partial_read+ref"
                     else "ref_imputed",
    RefMatch_Observed = ref_match_observed, RefMatch_Total = ref_match_total,
    ITDReadCoverage = itd_coverage_percent, ITDCoverageRLE = itd_coverage_rle,
    SupportConsistency = consistency_score, AlignmentScore = alignment_score_total,
    TotalSupportBases = total_support_bases, Orientation = orientation,
    LeftSoftclipCount = left_sc_count, RightSoftclipCount = right_sc_count,
    BothSoftclipCount = both_sc_count,
    LeftSoftclipPctSupport  = left_sc_pct_support,
    RightSoftclipPctSupport = right_sc_pct_support,
    LeftSoftclipPctWT  = left_sc_pct_wt,
    RightSoftclipPctWT = right_sc_pct_wt,
    LengthPE           = length_pe,
    LengthPE_NSpanning = length_pe_nspanning,
    LengthExt          = as.integer(length_ext),
    PESoftclipSupport  = pe_softclip_support,
    PESoftclipEventPairs = pe_softclip_event_pairs,
    PESoftclipLongPairs  = pe_softclip_long_pairs,
    PEOrientationFR = pe_orientation_fr,
    PEOrientationRF = pe_orientation_rf,
    PEOrientationFF = pe_orientation_ff,
    PEOrientationRR = pe_orientation_rr,
    PEOrientationOther = pe_orientation_other,
    PEOrientationDominant = pe_orientation_dominant,
    ArtifactSuspect = artifact_suspect
  )
}


# ---------------------------------------------------------------------------
# Apply all hard-threshold filters; returns TRUE if the call should be kept
# Modified to support asymmetric PTD mode and PTD-specific microhomology check
# ---------------------------------------------------------------------------
.apply_filters <- function(metrics, thresholds,
                            min_length = NULL, max_length = NULL,
                            ptd_mode = FALSE, ptd_allow_asymmetric = TRUE) {
  with(thresholds, {

    # ---- Restored from old code: reject self-fulfilling imputed artifacts ----
    # ArtifactSuspect is TRUE only when sequence was fully imputed from the
    # reference AND the resulting sequence matches the reference exactly —
    # i.e. the call is circular evidence of itself.
    if (isTRUE(metrics$ArtifactSuspect)) return(FALSE)

    # ---- Restored from old code: low alignment score (read-derived only) ----
    # AlignmentScore is NA for imputed sequences; the !is.na() guard means
    # this filter is bypassed for imputed calls, exactly as in the old code.
    # Guard additionally for SequenceImputed to be safe under the new metric.
    if (!is.na(metrics$AlignmentScore) &&
        !isTRUE(metrics$SequenceImputed) &&
        metrics$AlignmentScore < min_alignment_score) return(FALSE)

    if (!is.null(min_length) && !is.na(metrics$Length) &&
        metrics$Length > 0L && metrics$Length < min_length)
      return(FALSE)
    if (!is.null(max_length) && !is.na(metrics$Length) &&
        metrics$Length > 0L && metrics$Length > max_length)
      return(FALSE)
    if (!is.na(metrics$CoverageDrop)         && metrics$CoverageDrop         < min_coverage_drop)       return(FALSE)
    if (!is.na(metrics$LocalCoverage)        && metrics$LocalCoverage        < min_local_coverage)      return(FALSE)
    if (!is.na(metrics$MedianMicrohomology)  && metrics$MedianMicrohomology  < min_microhomology)       return(FALSE)
    if (!is.na(metrics$DiscordantRatio)      && metrics$DiscordantRatio      < min_discordant_ratio)    return(FALSE)
    if (!is.na(metrics$RepeatEntropy)        && metrics$RepeatEntropy        < min_entropy)              return(FALSE)
    if (!is.na(metrics$StrandBias)           && (metrics$StrandBias < min_strand_bias ||
                                                  metrics$StrandBias > max_strand_bias))                return(FALSE)
    if (!is.na(metrics$MeanSupportMAPQ)      && metrics$MeanSupportMAPQ      < min_mean_support_mapq)  return(FALSE)
    if (!is.na(metrics$BreakpointSpread)     && metrics$BreakpointSpread     > max_breakpoint_spread)   return(FALSE)
    if (!is.na(metrics$SoftclipFraction)     && metrics$SoftclipFraction     < min_softclip_fraction)  return(FALSE)
    if (!is.na(metrics$UniqueBreakpoints)    && metrics$UniqueBreakpoints    < min_unique_breakpoints)  return(FALSE)
    if (!is.na(metrics$ITDReadCoverage)      && metrics$ITDReadCoverage      < min_itd_read_coverage)  return(FALSE)

    # NOTE: The broad imputed-length filter (SequenceImputed && Length 41-500)
    # that was present in the first refactor has been removed. It was not in
    # the reference (old) code and caused significant sensitivity loss.
    # Self-fulfilling imputed calls are now caught exclusively by ArtifactSuspect.

    # ---- Soft-clip side filters (ITD vs PTD) ----
    left  <- metrics$LeftSoftclipCount
    right <- metrics$RightSoftclipCount
    if (!is.na(left) && !is.na(right)) {
      if (ptd_mode && ptd_allow_asymmetric) {
        # PTD asymmetric mode: require at least one side to have sufficient reads
        max_side <- max(left, right)
        if (max_side < min_side_softclip_reads) return(FALSE)
        # Check dominant side's percentage of support and WT
        if (left >= right) {
          if (!is.na(metrics$LeftSoftclipPctSupport) && metrics$LeftSoftclipPctSupport < min_softclip_pct_side) return(FALSE)
          if (!is.na(metrics$LeftSoftclipPctWT) && metrics$LeftSoftclipPctWT < min_left_softclip_pct_wt) return(FALSE)
        } else {
          if (!is.na(metrics$RightSoftclipPctSupport) && metrics$RightSoftclipPctSupport < min_softclip_pct_side) return(FALSE)
          if (!is.na(metrics$RightSoftclipPctWT) && metrics$RightSoftclipPctWT < min_right_softclip_pct_wt) return(FALSE)
        }
      } else {
        # Original symmetric ITD/PTD filter
        if (left < min_side_softclip_reads || right < min_side_softclip_reads) return(FALSE)
        lw <- metrics$LeftSoftclipPctWT;  rw <- metrics$RightSoftclipPctWT
        if (!is.na(lw) && lw < min_left_softclip_pct_wt)  return(FALSE)
        if (!is.na(rw) && rw < min_right_softclip_pct_wt) return(FALSE)
        if (!is.na(metrics$Length) && metrics$Length > 0L) {
          # ITD: additional checks already covered
        } else {
          # PTD symmetric (original)
          # BUGFIX: when both sides are 0, min(left,right)==0 -> 0/0 == NaN,
          # and `if (NaN > x)` throws "missing value where TRUE/FALSE needed".
          # Treat "both sides zero support" as a straightforward fail rather
          # than crashing the whole detection run.
          if (max(left, right) == 0L) return(FALSE)
          if (min(left, right) == 0L) return(FALSE)
          if (max(left, right) / min(left, right) > max_side_ratio) return(FALSE)
          if (left < min_abs_side_softclip || right < min_abs_side_softclip) return(FALSE)
          lp <- metrics$LeftSoftclipPctSupport;  rp <- metrics$RightSoftclipPctSupport
          if (!is.na(lp) && lp < min_softclip_pct_side) return(FALSE)
          if (!is.na(rp) && rp < min_softclip_pct_side) return(FALSE)
        }
      }
    }
    TRUE
  })
}

# ---------------------------------------------------------------------------
# Helper: mode of numeric vector (returns most frequent value)
# ---------------------------------------------------------------------------
.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# ---------------------------------------------------------------------------
# Debug helper: render one candidate's full metrics as a single log line.
# Called from detect_itd() (main.R) right after .compute_variant_metrics() /
# .apply_filters(), when debug = TRUE, so that EVERY candidate cluster is
# written to the log -- not just the ones that survive filtering. This is
# what lets you see, post-hoc, exactly why a real variant got dropped.
# ---------------------------------------------------------------------------
.format_candidate_debug <- function(metrics, passed) {
  fields <- c(
    "GenomicPosition", "Length", "LengthExt", "SupportingReads", "WildtypeReads",
    "DepthAtBreakpoint", "AlleleFrequency", "CoverageDrop", "LocalCoverage", "MedianMicrohomology",
    "DiscordantRatio", "RepeatEntropy", "StrandBias", "MeanSupportMAPQ",
    "BreakpointSpread", "SoftclipFraction", "UniqueBreakpoints", "ITDReadCoverage",
    "AlignmentScore", "SequenceImputed", "SequencePartial", "ArtifactSuspect",
    "LeftSoftclipCount", "RightSoftclipCount", "LeftSoftclipPctSupport",
    "RightSoftclipPctSupport", "LeftSoftclipPctWT", "RightSoftclipPctWT"
  )
  kv <- vapply(fields, function(f) {
    v <- metrics[[f]]
    if (is.null(v) || (length(v) == 1L && is.na(v))) return(paste0(f, "=NA"))
    if (is.numeric(v)) v <- format(round(v, 4), trim = TRUE)
    paste0(f, "=", v)
  }, character(1))
  sprintf("[CANDIDATE] %s -> %s", paste(kv, collapse = " "),
          if (isTRUE(passed)) "KEPT" else "FILTERED")
}
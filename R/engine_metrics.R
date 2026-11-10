.compute_left_right_coverage <- function(
  wildtype_info,
  len_specific_bp,
  support_qnames,
  buffer = 10L
) {
  strict_idx <- which(
    wildtype_info$is_primary &
      GenomicRanges::start(wildtype_info$gr) <= (len_specific_bp - buffer) &
      GenomicRanges::end(wildtype_info$gr) >= (len_specific_bp + buffer)
  )
  strict_qnames <- unique(wildtype_info$qnames[strict_idx])
  pure_wt <- setdiff(strict_qnames, support_qnames)
  wt_count <- length(pure_wt)
  list(left = wt_count, right = wt_count, avg = wt_count)
}


.calculate_vaf_and_depth <- function(
  raw_support,
  wildtype_info,
  len_specific_bp,
  support_qnames,
  best_len,
  nominal_read_len,
  max_correction,
  buffer = 10L
) {
  lr <- .compute_left_right_coverage(
    wildtype_info,
    len_specific_bp,
    support_qnames,
    buffer
  )
  strict_wt <- lr$avg
  len_used <- if (is.na(best_len)) 0L else best_len
  size_scalar <- min(1 + (0.25 * len_used / nominal_read_len), max_correction)
  corrected_support <- round(raw_support * size_scalar)
  depth_at_bp <- corrected_support + strict_wt
  corrected_af <- if (depth_at_bp == 0L) 0 else corrected_support / depth_at_bp
  list(
    corrected_support = corrected_support,
    depth_at_bp = depth_at_bp,
    raw_af = corrected_af,
    raw_depth = raw_support + strict_wt,
    left_wt = strict_wt,
    right_wt = strict_wt,
    avg_wt = strict_wt
  )
}


# Longest duplication on which sequence-level work (reference imputation,
# pairwise ref / orientation alignment) is attempted.
.MAX_SEQUENCE_OPS_LEN <- 2000L

# Keep one representative alignment row per read before calculating
# row-count-based metrics (soft-clip counts, strand bias, MAPQ, sequence
# consensus, etc.). Support is defined by distinct QNAMEs; duplicate
# evidence rows from one read must not inflate those metrics.
.deduplicate_support_rows <- function(support_rows, target_bp) {
  if (is.null(support_rows) || !is.data.frame(support_rows) ||
      !"read_name" %in% names(support_rows)) {
    return(support_rows)
  }
  qname <- as.character(support_rows$read_name)
  valid <- !is.na(qname) & nzchar(qname)
  support_rows <- support_rows[valid, , drop = FALSE]
  qname <- qname[valid]
  if (nrow(support_rows) < 2L || !anyDuplicated(qname)) return(support_rows)

  bp <- if ("breakpoint" %in% names(support_rows)) {
    suppressWarnings(as.numeric(support_rows$breakpoint))
  } else rep(NA_real_, nrow(support_rows))
  distance <- abs(bp - as.numeric(target_bp)[1L])
  distance[!is.finite(distance)] <- Inf

  mapq <- if ("mapq" %in% names(support_rows)) {
    suppressWarnings(as.numeric(support_rows$mapq))
  } else rep(0, nrow(support_rows))
  mapq[!is.finite(mapq)] <- -Inf

  has_seq <- if ("itd_seq" %in% names(support_rows)) {
    !is.na(support_rows$itd_seq) & nzchar(as.character(support_rows$itd_seq))
  } else rep(FALSE, nrow(support_rows))

  # Prefer the closest breakpoint, then highest MAPQ, then a row carrying
  # an observed sequence. The final index makes ties deterministic.
  ord <- order(distance, -mapq, !has_seq, seq_len(nrow(support_rows)))
  chosen <- ord[!duplicated(qname[ord])]
  support_rows[sort(chosen), , drop = FALSE]
}

.compute_variant_metrics <- function(
  cluster_bps,
  best_len,
  support_rows,
  genomic_start,
  ref_dna,
  gene_config,
  all_reads_cov,
  all_pairs,
  wildtype_info,
  ptd_mode,
  min_support,
  min_wt_reads,
  nominal_read_len,
  max_correction,
  vaf_threshold,
  do_alignment_score,
  do_support_bases,
  do_consistency,
  do_itd_coverage,
  do_coverage_drop,
  do_microhomology,
  do_repeat_entropy,
  do_discordant_ratio,
  do_detect_orientation,
  do_hgvs,
  max_pairwise_alignments,
  span_flank = 300L,
  length_ext = NA_integer_,
  pre_assembled_seq = NA_character_,
  max_itd_length = 1000L,
  debug = FALSE,
  verbose = FALSE
) {
  if (verbose) {
    message(
      "[TALOS] Computing metrics for cluster at ",
      paste(range(cluster_bps), collapse = "-"),
      " length=",
      best_len
    )
  }

  consistency_score <- NA_real_
  itd_coverage_percent <- NA_real_
  itd_coverage_rle <- NA_character_
  ref_match_observed <- NA_real_
  ref_match_total <- NA_real_
  alignment_score_total <- NA_real_
  total_support_bases <- NA_integer_
  coverage_drop <- NA_real_
  local_coverage <- NA_real_
  span_min_coverage <- NA_real_
  span_mean_coverage <- NA_real_
  span_depth_fold_change <- NA_real_
  median_microhomology <- NA_real_
  repeat_entropy <- NA_real_
  discordant_ratio <- NA_real_
  length_pe <- NA_real_
  length_pe_nspanning <- NA_integer_
  pe_softclip_support <- NA_integer_
  pe_softclip_event_pairs <- NA_integer_
  pe_softclip_long_pairs <- NA_integer_
  pe_orientation_fr <- NA_integer_
  pe_orientation_rf <- NA_integer_
  pe_orientation_ff <- NA_integer_
  pe_orientation_rr <- NA_integer_
  pe_orientation_other <- NA_integer_
  pe_orientation_dominant <- NA_character_
  orientation <- NA_character_
  hgvs <- list(c_notation = NA_character_, p_notation = NA_character_)

  len_specific_bp <- round(median(cluster_bps, na.rm = TRUE))
  if (is.na(len_specific_bp) || is.infinite(len_specific_bp)) {
    return(NULL)
  }
  support_rows <- .deduplicate_support_rows(support_rows, len_specific_bp)
  if (is.null(support_rows) || nrow(support_rows) == 0L) return(NULL)
  support_qnames <- unique(support_rows$read_name)
  raw_support <- length(support_qnames)
  # Checked before any alignment work (it used to come after the pairwise
  # alignments, which were then thrown away).
  if (raw_support < min_support) {
    return(NULL)
  }
  # PTD-scale events (graph PTDs, multi-kb): the read-level sequence is a
  # junction consensus, not the duplicated segment. Reference imputation of
  # a multi-kb segment, positional ref matching and segment-vs-segment
  # orientation alignment are meaningless (and quadratic) there.
  # Capped at .MAX_SEQUENCE_OPS_LEN: with max_itd_length = 30000 (KMT2A YAML)
  # a 20 kb call would otherwise get a 20 kb imputed sequence and a
  # 20 kb x 20 kb orientation alignment.
  seq_ops_limit <- min(c(max_itd_length, .MAX_SEQUENCE_OPS_LEN), na.rm = TRUE)
  is_ptd_scale <- isTRUE(ptd_mode) || (!is.na(best_len) && best_len > seq_ops_limit)
  local_seg_start <- len_specific_bp - genomic_start + 1L
  seg_in_ref <- !is.na(best_len) && best_len > 0L && local_seg_start >= 1L &&
    local_seg_start + best_len - 1L <= nchar(ref_dna)

  left_sc_count <- 0L
  right_sc_count <- 0L
  both_sc_count <- 0L
  support_cigars <- if (is.null(support_rows$cigar)) {
    rep(NA_character_, nrow(support_rows))
  } else {
    support_rows$cigar
  }
  support_clip_lengths <- lapply(support_cigars, .get_softclip_lengths)
  has_left_clip <- vapply(
    support_clip_lengths,
    function(x) x$lead > 0L,
    logical(1L)
  )
  has_right_clip <- vapply(
    support_clip_lengths,
    function(x) x$trail > 0L,
    logical(1L)
  )
  if (nrow(support_rows) > 0L) {
    both_sc_count <- sum(has_left_clip & has_right_clip, na.rm = TRUE)
    left_sc_count <- sum(has_left_clip, na.rm = TRUE)
    right_sc_count <- sum(has_right_clip, na.rm = TRUE)
  }

  itd_seq <- NA_character_
  imputed <- TRUE
  sequence_partial <- FALSE
  observed_seq <- NA_character_
  observed_len <- 0L

  valid_seqs <- support_rows$itd_seq[
    !is.na(support_rows$itd_seq) &
      nchar(support_rows$itd_seq) > 0L
  ]
  if (length(valid_seqs) > 0L) {
    itd_seq <- names(sort(table(valid_seqs), decreasing = TRUE))[1L]
    imputed <- FALSE
    observed_seq <- itd_seq
    observed_len <- nchar(itd_seq)
  }

  if (is.na(itd_seq)) {
    extracted_seqs <- vector("character", 2L * nrow(support_rows))
    seq_idx <- 1L
    for (i in seq_len(nrow(support_rows))) {
      sc <- .get_softclips(support_rows$cigar[i], support_rows$read_seq[i])
      if (!is.na(sc$lead)) {
        extracted_seqs[seq_idx] <- sc$lead
        seq_idx <- seq_idx + 1L
      }
      if (!is.na(sc$trail)) {
        extracted_seqs[seq_idx] <- sc$trail
        seq_idx <- seq_idx + 1L
      }
    }
    extracted_seqs <- extracted_seqs[seq_len(seq_idx - 1L)]
    if (length(extracted_seqs) > 0L) {
      itd_seq <- names(sort(table(extracted_seqs), decreasing = TRUE))[1L]
      imputed <- FALSE
      observed_seq <- itd_seq
      observed_len <- nchar(itd_seq)
    }
  }

  if (
    is.na(itd_seq) && !is.na(pre_assembled_seq) && nchar(pre_assembled_seq) > 0L
  ) {
    itd_seq <- pre_assembled_seq
    imputed <- FALSE
    observed_seq <- itd_seq
    observed_len <- nchar(itd_seq)
  }

  if (
    !is.na(itd_seq) &&
      !is.na(best_len) &&
      best_len > 0L &&
      nchar(itd_seq) > best_len
  ) {
    itd_seq <- substr(itd_seq, 1L, best_len)
    if (!is.na(observed_seq) && nchar(observed_seq) > best_len) {
      observed_seq <- substr(observed_seq, 1L, best_len)
    }
    observed_len <- min(observed_len, best_len)
  }

  if (
    !is.na(itd_seq) &&
      !is.na(best_len) &&
      best_len > 0L &&
      nchar(itd_seq) < best_len &&
      !is_ptd_scale && seg_in_ref
  ) {
    partial_obs_len <- nchar(itd_seq)
    ref_local_start <- len_specific_bp - genomic_start + partial_obs_len + 1L
    ref_local_end <- ref_local_start + (best_len - partial_obs_len) - 1L
    if (ref_local_end <= nchar(ref_dna)) {
      itd_seq <- paste0(
        itd_seq,
        substr(ref_dna, ref_local_start, ref_local_end)
      )
      sequence_partial <- TRUE
      imputed <- TRUE
      observed_len <- partial_obs_len
    }
  }

  if (is.na(itd_seq) && !is_ptd_scale && seg_in_ref) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end <- local_start + best_len - 1L
    if (!is.na(local_end) && local_end <= nchar(ref_dna)) {
      itd_seq <- substr(ref_dna, local_start, local_end)
      imputed <- TRUE
      sequence_partial <- FALSE
      observed_seq <- NA_character_
      observed_len <- 0L
    }
  }

  if (
    !is.na(itd_seq) &&
      nchar(itd_seq) > 0L &&
      nrow(support_rows) > 0L &&
      (do_consistency || do_itd_coverage)
  ) {
    valid_mask <- !is.na(support_rows$read_seq) &
      nchar(support_rows$read_seq) > 0L
    v_seqs <- support_rows$read_seq[valid_mask]
    if (length(v_seqs) > 0L) {
      if (
        (do_consistency || do_itd_coverage) &&
          .has_pwalign &&
          .has_biostrings &&
          length(v_seqs) > max_pairwise_alignments
      ) {
        old_seed <- if (
          exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        ) {
          .GlobalEnv$.Random.seed
        } else {
          NULL
        }
        set.seed(42L)
        v_seqs <- v_seqs[sample(length(v_seqs), max_pairwise_alignments)]
        if (!is.null(old_seed)) {
          .GlobalEnv$.Random.seed <- old_seed
        } else {
          rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        }
      }
      alignments_needed <- do_consistency || do_itd_coverage
      if (alignments_needed && .has_pwalign && .has_biostrings) {
        aln <- pwalign::pairwiseAlignment(
          pattern = Biostrings::DNAStringSet(v_seqs),
          subject = Biostrings::DNAString(itd_seq),
          type = "local"
        )
        pids <- pwalign::pid(aln)
        if (do_consistency) {
          consistency_score <- sum(pids >= 90, na.rm = TRUE) /
            length(v_seqs) *
            100
        }
        if (do_itd_coverage) {
          good_idx <- which(pids >= 90)
          if (length(good_idx) > 0L) {
            subj_ranges <- pwalign::subject(aln[good_idx])
            cov_vec <- integer(nchar(itd_seq))
            s_st <- BiocGenerics::start(subj_ranges)
            s_en <- BiocGenerics::end(subj_ranges)
            for (ki in seq_along(good_idx)) {
              cov_vec[s_st[ki]:s_en[ki]] <- cov_vec[s_st[ki]:s_en[ki]] + 1L
            }
            itd_coverage_percent <- sum(cov_vec > 0L) / length(cov_vec) * 100
            rle_res <- rle(cov_vec)
            itd_coverage_rle <- paste(
              paste0(rle_res$values, ":", rle_res$lengths),
              collapse = ","
            )
          } else {
            itd_coverage_percent <- 0
            itd_coverage_rle <- paste0("0:", nchar(itd_seq))
          }
        }
      } else {
        if (do_consistency) {
          consistency_score <- sum(grepl(itd_seq, v_seqs, fixed = TRUE)) /
            length(v_seqs) *
            100
        }
        if (do_itd_coverage) {
          itd_coverage_percent <- if (
            any(grepl(itd_seq, v_seqs, fixed = TRUE))
          ) {
            100
          } else {
            0
          }
          itd_coverage_rle <- paste0(
            as.integer(itd_coverage_percent > 0),
            ":",
            nchar(itd_seq)
          )
        }
      }
    }
  }

  if (do_alignment_score && !is.na(itd_seq) && !is_ptd_scale && seg_in_ref) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end <- local_start + best_len - 1L
    if (!is.na(local_end) && local_end <= nchar(ref_dna)) {
      ref_segment <- substr(ref_dna, local_start, local_end)
      min_len_tot <- min(nchar(itd_seq), nchar(ref_segment))
      if (min_len_tot > 0L) {
        alignment_score_total <-
          if (itd_seq == ref_segment) {
            1.0
          } else {
            sum(
              charToRaw(substr(itd_seq, 1L, min_len_tot)) ==
                charToRaw(substr(ref_segment, 1L, min_len_tot))
            ) /
              min_len_tot
          }
        ref_match_total <- round(alignment_score_total * 100, 1)
        if (!is.na(observed_seq) && observed_len > 0L) {
          obs_len <- min(observed_len, min_len_tot)
          if (obs_len > 0L) {
            obs_sub <- substr(observed_seq, 1L, obs_len)
            ref_sub <- substr(ref_segment, 1L, obs_len)
            alignment_score_obs <- if (obs_sub == ref_sub) {
              1.0
            } else {
              sum(charToRaw(obs_sub) == charToRaw(ref_sub)) / obs_len
            }
            ref_match_observed <- round(alignment_score_obs * 100, 1)
          }
        }
      }
    }
  }

  if (do_support_bases && nrow(support_rows) > 0L) {
    lead_S <- vapply(support_clip_lengths, `[[`, integer(1L), "lead")
    trail_S <- vapply(support_clip_lengths, `[[`, integer(1L), "trail")
    total_support_bases <- sum(lead_S, trail_S, na.rm = TRUE)
  }

  vaf_metrics <- .calculate_vaf_and_depth(
    raw_support,
    wildtype_info,
    len_specific_bp,
    support_qnames,
    best_len,
    nominal_read_len,
    max_correction,
    buffer = 10L
  )
  corrected_support <- vaf_metrics$corrected_support
  depth_at_bp <- vaf_metrics$depth_at_bp
  raw_af <- vaf_metrics$raw_af
  wildtype_reads <- round(vaf_metrics$avg_wt)

  if (wildtype_reads < min_wt_reads) {
    return(NULL)
  }
  if (raw_af < vaf_threshold) {
    return(NULL)
  }

  left_sc_pct_wt <- if (wildtype_reads > 0L) {
    (left_sc_count / wildtype_reads) * 100
  } else {
    NA_real_
  }
  right_sc_pct_wt <- if (wildtype_reads > 0L) {
    (right_sc_count / wildtype_reads) * 100
  } else {
    NA_real_
  }

  strand_bias <- if (nrow(support_rows) > 0L) {
    round(mean(support_rows$is_reverse, na.rm = TRUE), 4)
  } else {
    NA_real_
  }
  mean_support_mapq <- if (nrow(support_rows) > 0L) {
    round(mean(support_rows$mapq, na.rm = TRUE), 1)
  } else {
    NA_real_
  }
  bp_spread <- if (nrow(support_rows) > 1L) {
    as.integer(max(support_rows$breakpoint) - min(support_rows$breakpoint))
  } else {
    0L
  }
  softclip_frac <- if (nrow(support_rows) > 0L) {
    round(mean(has_left_clip | has_right_clip, na.rm = TRUE), 4)
  } else {
    NA_real_
  }
  unique_bps <- length(unique(support_rows$breakpoint))
  left_sc_pct_support <- if (corrected_support > 0L) {
    (left_sc_count / corrected_support) * 100
  } else {
    NA_real_
  }
  right_sc_pct_support <- if (corrected_support > 0L) {
    (right_sc_count / corrected_support) * 100
  } else {
    NA_real_
  }

  if (do_coverage_drop && !is.na(best_len) && best_len > 0L) {
    cov_ok <- FALSE
    if (!is.null(all_reads_cov)) {
      coverage_drop <- compute_coverage_drop(
        all_reads_cov,
        gene_config$chrom,
        len_specific_bp,
        flank = 200L,
        verbose = debug
      )
      if (!is.na(coverage_drop)) cov_ok <- TRUE
    }
    if (
      !cov_ok && !is.null(wildtype_info$gr) && length(wildtype_info$gr) > 0L
    ) {
      if (debug) {
        message("Coverage drop: recomputing from wildtype reads")
      }
      reg_gr <- GenomicRanges::GRanges(
        seqnames = gene_config$chrom,
        ranges = IRanges::IRanges(
          start = gene_config$genomic_start,
          end = gene_config$genomic_end
        )
      )
      ov <- GenomicRanges::findOverlaps(wildtype_info$gr, reg_gr)
      if (length(ov) > 0L) {
        reads_sub <- wildtype_info$gr[unique(S4Vectors::queryHits(ov))]
        cov_fallback <- GenomicAlignments::coverage(reads_sub)[[
          gene_config$chrom
        ]]
        if (!is.null(cov_fallback)) {
          coverage_drop <- compute_coverage_drop(
            cov_fallback,
            gene_config$chrom,
            len_specific_bp,
            flank = 200L,
            verbose = debug
          )
        }
      }
    }
    if (debug && is.na(coverage_drop)) {
      message("Coverage drop: still NA after fallback")
    }
  }

  if (!is.null(all_reads_cov)) {
    local_coverage <- compute_local_coverage(
      all_reads_cov,
      gene_config$chrom,
      len_specific_bp,
      buffer = 10L,
      verbose = debug
    )
  }

  if (!is.null(all_reads_cov) && !is.na(best_len) && best_len > 0L) {
    span_cov <- .compute_span_coverage_summary(
      all_reads_cov,
      gene_config$chrom,
      len_specific_bp,
      best_len,
      flank = span_flank
    )
    span_min_coverage <- span_cov$min_cov
    span_mean_coverage <- span_cov$mean_cov
    span_depth_fold_change <- span_cov$depth_fold_change
  }

  if (do_microhomology) {
    median_microhomology <- compute_microhomology(
      support_rows,
      ref_dna,
      len_specific_bp,
      genomic_start,
      debug = debug
    )
  }

  if (do_repeat_entropy) {
    repeat_entropy <- compute_repeat_entropy(
      ref_dna,
      len_specific_bp,
      genomic_start
    )
  }

  if (!is.null(all_pairs)) {
    if (do_discordant_ratio) {
      discordant_ratio <- compute_discordant_ratio(all_pairs, len_specific_bp)
    }

    pe_est <- compute_pe_itd_length(
      all_pairs,
      len_specific_bp,
      nominal_read_len = nominal_read_len
    )
    length_pe <- pe_est$length_pe
    length_pe_nspanning <- as.integer(pe_est$n_spanning)
    if (verbose && !is.na(length_pe)) {
      message(
        "[TALOS] PE insert-size ITD length estimate: ",
        length_pe,
        " bp (",
        length_pe_nspanning,
        " spanning pairs)"
      )
    }

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

  if (do_detect_orientation && !is.na(itd_seq) && !is_ptd_scale && seg_in_ref) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end <- local_start + best_len - 1L
    if (!is.na(local_end) && local_end <= nchar(ref_dna)) {
      orientation <- detect_orientation(
        itd_seq,
        substr(ref_dna, local_start, local_end)
      )
    }
  }

  exon_ranges <- .get_gene_exons(gene_config)
  is_exonic <- .is_breakpoint_exonic(
    len_specific_bp,
    gene_config$chrom,
    exon_ranges
  )
  if (do_hgvs) {
    # The called length wins: for PTD-scale events itd_seq is only the
    # junction consensus, so its length is not the duplication length.
    dup_len_for_hgvs <- if (!is.na(best_len) && best_len > 0L) {
      best_len
    } else if (!is.na(itd_seq) && nchar(itd_seq) > 0L) {
      nchar(itd_seq)
    } else {
      NA_integer_
    }
    if (!is.null(gene_config$all_exons) || !is.null(exon_ranges)) {
      hgvs <- compute_hgvs_annotations(
        gene_config,
        len_specific_bp,
        itd_seq,
        dup_len_for_hgvs,
        debug
      )
      if (!is_exonic && !is.na(hgvs$p_notation)) {
        hgvs$p_notation <- "intronic"
        if (is.na(hgvs$c_notation)) hgvs$c_notation <- "c.?+?"
      }
    }
  }

  dup_exons <- if (!is.null(support_rows$dup_exons)) {
    de <- support_rows$dup_exons[!is.na(support_rows$dup_exons)]
    if (length(de)) .mode(de) else NA_character_
  } else {
    NA_character_
  }
  candidate_types <- if (!is.null(support_rows$type)) {
    paste(sort(unique(stats::na.omit(support_rows$type))), collapse = "+")
  } else {
    NA_character_
  }

  artifact_suspect <- isTRUE(imputed) &&
    !isTRUE(sequence_partial) &&
    !is.na(ref_match_total) &&
    ref_match_total >= 99.0

  if (debug && artifact_suspect) {
    message(sprintf(
      "  [ARTIFACT] bp=%d fully imputed + 100%% ref match",
      len_specific_bp
    ))
  }

  if (verbose) {
    message(
      "[TALOS] Metrics computed: support=",
      corrected_support,
      ", AF=",
      round(raw_af, 4),
      ", length=",
      best_len
    )
  }

  list(
    Sample = wildtype_info$sample_name,
    Gene = gene_config$gene,
    Genome = if (is.null(gene_config$build)) "custom" else gene_config$build,
    GenomicPosition = len_specific_bp,
    ITD_Sequence = itd_seq,
    Length = best_len,
    DupExons = dup_exons,
    CandidateTypes = candidate_types,
    SupportingReads = corrected_support,
    WildtypeReads = wildtype_reads,
    DepthAtBreakpoint = depth_at_bp,
    AlleleFrequency = raw_af,
    HGVS_cDNA = hgvs$c_notation,
    HGVS_Protein = hgvs$p_notation,
    StrandBias = strand_bias,
    MeanSupportMAPQ = mean_support_mapq,
    BreakpointSpread = bp_spread,
    SoftclipFraction = softclip_frac,
    UniqueBreakpoints = unique_bps,
    CoverageDrop = coverage_drop,
    LocalCoverage = local_coverage,
    SpanMinCoverage = span_min_coverage,
    SpanMeanCoverage = span_mean_coverage,
    SpanDepthFoldChange = span_depth_fold_change,
    MedianMicrohomology = median_microhomology,
    DiscordantRatio = discordant_ratio,
    RepeatEntropy = repeat_entropy,
    SequenceImputed = imputed,
    SequencePartial = sequence_partial,
    SequenceSource = if (!imputed) {
      "observed"
    } else if (sequence_partial) {
      "partial_read+ref"
    } else {
      "ref_imputed"
    },
    RefMatch_Observed = ref_match_observed,
    RefMatch_Total = ref_match_total,
    ITDReadCoverage = itd_coverage_percent,
    ITDCoverageRLE = itd_coverage_rle,
    SupportConsistency = consistency_score,
    AlignmentScore = alignment_score_total,
    TotalSupportBases = total_support_bases,
    Orientation = orientation,
    LeftSoftclipCount = left_sc_count,
    RightSoftclipCount = right_sc_count,
    BothSoftclipCount = both_sc_count,
    LeftSoftclipPctSupport = left_sc_pct_support,
    RightSoftclipPctSupport = right_sc_pct_support,
    LeftSoftclipPctWT = left_sc_pct_wt,
    RightSoftclipPctWT = right_sc_pct_wt,
    LengthPE = length_pe,
    LengthPE_NSpanning = length_pe_nspanning,
    LengthExt = as.integer(length_ext),
    PESoftclipSupport = pe_softclip_support,
    PESoftclipEventPairs = pe_softclip_event_pairs,
    PESoftclipLongPairs = pe_softclip_long_pairs,
    PEOrientationFR = pe_orientation_fr,
    PEOrientationRF = pe_orientation_rf,
    PEOrientationFF = pe_orientation_ff,
    PEOrientationRR = pe_orientation_rr,
    PEOrientationOther = pe_orientation_other,
    PEOrientationDominant = pe_orientation_dominant,
    ArtifactSuspect = artifact_suspect
  )
}


# Threshold rules: metric column, threshold name, direction. A rule only
# applies when the metric is not NA (same semantics as before, one table
# instead of seventeen hand-written blocks).
.FILTER_RULES <- list(
  c("CoverageDrop", "min_coverage_drop", "min"),
  c("LocalCoverage", "min_local_coverage", "min"),
  c("SpanMinCoverage", "min_span_coverage", "min"),
  c("SpanDepthFoldChange", "min_span_depth_fold_change", "min"),
  c("SpanDepthFoldChange", "max_span_depth_fold_change", "max"),
  c("MedianMicrohomology", "min_microhomology", "min"),
  c("DiscordantRatio", "min_discordant_ratio", "min"),
  c("RepeatEntropy", "min_entropy", "min"),
  c("StrandBias", "min_strand_bias", "min"),
  c("StrandBias", "max_strand_bias", "max"),
  c("MeanSupportMAPQ", "min_mean_support_mapq", "min"),
  c("BreakpointSpread", "max_breakpoint_spread", "max"),
  c("SoftclipFraction", "min_softclip_fraction", "min"),
  c("UniqueBreakpoints", "min_unique_breakpoints", "min"),
  c("ITDReadCoverage", "min_itd_read_coverage", "min")
)

#' Apply the candidate filters
#' @return TRUE/FALSE; when FALSE the attribute "reason" names the first
#'   failing rule (shown by the debug log).
.apply_filters <- function(
  metrics,
  thresholds,
  min_length = NULL,
  max_length = NULL,
  ptd_mode = FALSE,
  ptd_allow_asymmetric = TRUE
) {
  fail <- function(reason) structure(FALSE, reason = reason)
  num <- function(x) if (is.null(x) || length(x) == 0L) NA_real_ else suppressWarnings(as.numeric(x[[1L]]))
  th <- function(name) thresholds[[name]] %||% (if (startsWith(name, "max")) Inf else 0)

  if (isTRUE(metrics$ArtifactSuspect)) return(fail("ArtifactSuspect"))
  aln <- num(metrics$AlignmentScore)
  if (!is.na(aln) && !isTRUE(metrics$SequenceImputed) && aln < th("min_alignment_score")) {
    return(fail("min_alignment_score"))
  }
  len <- num(metrics$Length)
  if (!is.na(len) && len > 0) {
    if (!is.null(min_length) && len < min_length) return(fail("min_length"))
    if (!is.null(max_length) && len > max_length) return(fail("max_length"))
  }
  for (r in .FILTER_RULES) {
    v <- num(metrics[[r[[1L]]]])
    if (is.na(v)) next
    t <- th(r[[2L]])
    if ((r[[3L]] == "min" && v < t) || (r[[3L]] == "max" && v > t)) return(fail(r[[2L]]))
  }

  left <- num(metrics$LeftSoftclipCount)
  right <- num(metrics$RightSoftclipCount)
  if (is.na(left) || is.na(right)) return(TRUE)
  below <- function(col, thr) { v <- num(metrics[[col]]); !is.na(v) && v < th(thr) }

  if (ptd_mode && ptd_allow_asymmetric) {
    if (max(left, right) < th("min_side_softclip_reads")) return(fail("min_side_softclip_reads"))
    side <- if (left >= right) "Left" else "Right"
    if (below(paste0(side, "SoftclipPctSupport"), "min_softclip_pct_side")) return(fail("min_softclip_pct_side"))
    if (below(paste0(side, "SoftclipPctWT"), paste0("min_", tolower(side), "_softclip_pct_wt"))) {
      return(fail(paste0("min_", tolower(side), "_softclip_pct_wt")))
    }
    return(TRUE)
  }

  if (left < th("min_side_softclip_reads") || right < th("min_side_softclip_reads")) {
    return(fail("min_side_softclip_reads"))
  }
  if (below("LeftSoftclipPctWT", "min_left_softclip_pct_wt")) return(fail("min_left_softclip_pct_wt"))
  if (below("RightSoftclipPctWT", "min_right_softclip_pct_wt")) return(fail("min_right_softclip_pct_wt"))
  # Zero-length (PTD clip) calls must be supported on both sides.
  if (is.na(len) || len <= 0) {
    if (max(left, right) == 0 || min(left, right) == 0) return(fail("one_sided_softclips"))
    if (max(left, right) / min(left, right) > th("max_side_ratio")) return(fail("max_side_ratio"))
    if (left < th("min_abs_side_softclip") || right < th("min_abs_side_softclip")) {
      return(fail("min_abs_side_softclip"))
    }
    if (below("LeftSoftclipPctSupport", "min_softclip_pct_side") ||
        below("RightSoftclipPctSupport", "min_softclip_pct_side")) {
      return(fail("min_softclip_pct_side"))
    }
  }
  TRUE
}

.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

.format_candidate_debug <- function(metrics, passed) {
  fields <- c(
    "GenomicPosition",
    "Length",
    "LengthExt",
    "SupportingReads",
    "WildtypeReads",
    "DepthAtBreakpoint",
    "AlleleFrequency",
    "CoverageDrop",
    "LocalCoverage",
    "SpanMinCoverage",
    "SpanMeanCoverage",
    "SpanDepthFoldChange",
    "MedianMicrohomology",
    "DiscordantRatio",
    "RepeatEntropy",
    "StrandBias",
    "MeanSupportMAPQ",
    "BreakpointSpread",
    "SoftclipFraction",
    "UniqueBreakpoints",
    "ITDReadCoverage",
    "AlignmentScore",
    "SequenceImputed",
    "SequencePartial",
    "ArtifactSuspect",
    "LeftSoftclipCount",
    "RightSoftclipCount",
    "LeftSoftclipPctSupport",
    "RightSoftclipPctSupport",
    "LeftSoftclipPctWT",
    "RightSoftclipPctWT",
    "DupExons",
    "CandidateTypes"
  )
  kv <- vapply(
    fields,
    function(f) {
      v <- metrics[[f]]
      if (is.null(v) || (length(v) == 1L && is.na(v))) {
        return(paste0(f, "=NA"))
      }
      if (is.numeric(v)) {
        v <- format(round(v, 4), trim = TRUE)
      }
      paste0(f, "=", v)
    },
    character(1)
  )
  sprintf(
    "[CANDIDATE] %s -> %s",
    paste(kv, collapse = " "),
    if (isTRUE(passed)) {
      "KEPT"
    } else if (!is.null(attr(passed, "reason"))) {
      paste0("FILTERED (", attr(passed, "reason"), ")")
    } else {
      "FILTERED"
    }
  )
}

# ============================================================================
# TALOS – Detection Engine (Memory + Speed Optimized)
# ============================================================================

.has_fastmatch  <- requireNamespace("fastmatch",  quietly = TRUE)
.has_cigarillo  <- requireNamespace("cigarillo",  quietly = TRUE)
.has_pwalign    <- requireNamespace("pwalign",    quietly = TRUE)
.has_biostrings <- requireNamespace("Biostrings", quietly = TRUE)


.prepare_kmers <- function(ref_dna, k) {
  ref_len <- nchar(ref_dna)
  if (ref_len < k) stop("Reference shorter than k-mer length.")

  if (.has_biostrings) {
    dna_obj <- Biostrings::DNAString(ref_dna)
    v       <- Biostrings::Views(dna_obj, start = seq_len(ref_len - k + 1L), width = k)
    kmers   <- as.character(v)
  } else {
    kmer_starts <- seq_len(ref_len - k + 1L)
    kmers       <- substring(ref_dna, kmer_starts, kmer_starts + k - 1L)
  }

  if (.has_fastmatch) fastmatch::fmatch("", kmers)  # prime hash
  kmers
}

.kmer_match <- function(query, table) {
  if (.has_fastmatch) fastmatch::fmatch(query, table) else match(query, table)
}


# ----------------------------------------------------------------------------
# Streaming BAM loader (with coverage computation) – ENHANCED with logging
# ----------------------------------------------------------------------------
.load_bam_data_streaming <- function(bam_path, gene_config,
                                      compute_pairs = FALSE,
                                      max_reads = NULL, chunk_size = 50000,
                                      verbose = FALSE) {
  if (verbose) message("[TALOS] Loading BAM region: ", gene_config$chrom, ":", gene_config$genomic_start, "-", gene_config$genomic_end)

  target_chrom <- gene_config$chrom
  which_range <- GenomicRanges::GRanges(
    seqnames = target_chrom,
    ranges   = IRanges::IRanges(start = gene_config$genomic_start,
                                end   = gene_config$genomic_end)
  )
  param <- Rsamtools::ScanBamParam(
    which = which_range,
    what  = c("seq", "qname", "cigar", "flag", "mapq", "isize"),
    flag  = Rsamtools::scanBamFlag(isSecondaryAlignment    = FALSE,
                                   isSupplementaryAlignment = FALSE)
  )
  bam_file <- Rsamtools::BamFile(bam_path, yieldSize = chunk_size)

  all_reads_list <- list()
  total_reads    <- 0L
  chunk_count    <- 0L
  open(bam_file)
  repeat {
    chunk <- GenomicAlignments::readGAlignments(bam_file, param = param,
                                                use.names = TRUE)
    if (length(chunk) == 0) break
    chunk_count <- chunk_count + 1L
    total_reads <- total_reads + length(chunk)
    if (verbose) message("[TALOS] Loaded chunk ", chunk_count, " (", length(chunk), " reads, total ", total_reads, ")")
    all_reads_list[[length(all_reads_list) + 1L]] <- chunk
    if (!is.null(max_reads) && total_reads >= max_reads) {
      if (total_reads > max_reads) {
        keep <- max_reads - (total_reads - length(chunk))
        all_reads_list[[length(all_reads_list)]] <- chunk[seq_len(keep)]
        if (verbose) message("[TALOS] Stopping early at max_reads = ", max_reads)
      }
      break
    }
  }
  close(bam_file)
  all_reads <- do.call(c, all_reads_list)
  if (verbose) message("[TALOS] Total reads loaded: ", length(all_reads))
  if (length(all_reads) == 0) {
    if (verbose) message("[TALOS] No reads found in target region.")
    return(list(reads = all_reads, cov = NULL, pairs = NULL))
  }

  # ---- ENHANCED: robust chromosome name normalisation using mapSeqlevels ----
  current_levels <- GenomeInfoDb::seqlevels(all_reads)
  ucsc_name <- target_chrom
  if (!ucsc_name %in% current_levels) {
    converted <- GenomeInfoDb::mapSeqlevels(current_levels, "UCSC")
    if (any(!is.na(converted))) {
      new_levels <- ifelse(is.na(converted), current_levels, converted)
      all_reads <- GenomeInfoDb::renameSeqlevels(all_reads, new_levels)
      if (!ucsc_name %in% GenomeInfoDb::seqlevels(all_reads)) {
        warning("Chromosome name '", ucsc_name,
                "' not found after seqlevel conversion. Coverage will be NA.")
      }
    } else {
      warning("Chromosome name '", ucsc_name, "' not found in BAM seqlevels. Coverage will be NA.")
    }
  }
  # -----------------------------------------------------

  if (verbose) message("[TALOS] Computing coverage...")
  cov_rle_list <- GenomicAlignments::coverage(all_reads)
  cov <- cov_rle_list[[target_chrom]]

  pairs <- NULL
  if (compute_pairs) {
    if (verbose) message("[TALOS] Loading read pairs...")
    bam_file_pairs <- Rsamtools::BamFile(bam_path, yieldSize = chunk_size)
    pairs_list  <- list()
    total_pairs <- 0L
    chunk_count <- 0L
    open(bam_file_pairs)
    repeat {
      p_chunk <- suppressWarnings(GenomicAlignments::readGAlignmentPairs(bam_file_pairs, param = param))
      if (length(p_chunk) == 0) break
      chunk_count <- chunk_count + 1L
      total_pairs <- total_pairs + length(p_chunk)
      if (verbose) message("[TALOS] Loaded pair chunk ", chunk_count, " (", length(p_chunk), " pairs, total ", total_pairs, ")")
      pairs_list[[length(pairs_list) + 1L]] <- p_chunk
      if (!is.null(max_reads) && total_pairs >= max_reads) {
        if (total_pairs > max_reads) {
          keep_n <- max_reads - (total_pairs - length(p_chunk))
          pairs_list[[length(pairs_list)]] <- p_chunk[seq_len(keep_n)]
        }
        break
      }
    }
    close(bam_file_pairs)
    if (length(pairs_list) > 0) pairs <- do.call(c, pairs_list)
    if (verbose) message("[TALOS] Total pairs loaded: ", length(pairs))
  }

  list(reads = all_reads, cov = cov, pairs = pairs, cov_list = cov_rle_list)
}


# ----------------------------------------------------------------------------
# Pre-filter reads by CIGAR/MAPQ (with logging)
# ----------------------------------------------------------------------------
.filter_reads_by_cigar <- function(reads, min_mapq, min_ins_filter, verbose = FALSE) {
  if (verbose) message("[TALOS] Pre-filtering reads: min_mapq=", min_mapq, ", min_ins_filter=", min_ins_filter)
  flags  <- S4Vectors::mcols(reads)$flag
  mapqs  <- S4Vectors::mcols(reads)$mapq
  mapqs[is.na(mapqs)] <- 0L
  cigars  <- GenomicAlignments::cigar(reads)
  unmapped <- bitwAnd(flags, 0x4) != 0

  op_table <- if (.has_cigarillo)
    cigarillo::tabulate_cigar_ops(cigars)
  else
    suppressWarnings(GenomicAlignments::cigarOpTable(cigars))

  has_softclip <- if ("S" %in% colnames(op_table))
    op_table[, "S"] > 0L
  else
    rep(FALSE, nrow(op_table))

  net_ins <- if ("I" %in% colnames(op_table) && "D" %in% colnames(op_table)) {
    op_table[, "I"] - op_table[, "D"]
  } else if ("I" %in% colnames(op_table)) {
    op_table[, "I"]
  } else if ("D" %in% colnames(op_table)) {
    -op_table[, "D"]
  } else {
    rep(0L, nrow(op_table))
  }
  net_ins[is.na(net_ins)] <- 0L

  keep <- (unmapped | has_softclip | (net_ins >= min_ins_filter)) &
          (mapqs >= min_mapq)
  keep[is.na(keep)] <- FALSE
  filtered <- reads[keep]
  if (verbose) message("[TALOS] ", sum(keep), " reads passed pre-filter (", length(reads), " total)")
  filtered
}


.safe_qnames <- function(reads) {
  qn <- as.character(S4Vectors::mcols(reads)$qname)
  if (length(qn) == 0 || all(is.na(qn))) qn <- names(reads)
  na_idx <- which(is.na(qn))
  if (length(na_idx) > 0) qn[na_idx] <- paste0("read_", na_idx)
  qn
}


# ----------------------------------------------------------------------------
# PTD fast path (with logging)
# ----------------------------------------------------------------------------
.extract_candidates_ptd <- function(reads, genomic_start, ref_len, verbose = FALSE) {
  n_reads <- length(reads)
  if (n_reads == 0) return(list())
  if (verbose) message("[TALOS] Extracting PTD candidates from ", n_reads, " reads")

  all_qnames  <- .safe_qnames(reads)
  all_mapqs   <- S4Vectors::mcols(reads)$mapq; all_mapqs[is.na(all_mapqs)] <- 0L
  all_flags   <- S4Vectors::mcols(reads)$flag
  is_reverse  <- bitwAnd(all_flags, 0x10L) != 0L
  read_starts <- BiocGenerics::start(reads)
  read_ends   <- BiocGenerics::end(reads)
  read_cigars <- GenomicAlignments::cigar(reads)
  read_seqs   <- as.character(S4Vectors::mcols(reads)$seq)

  has_lead  <- grepl("^\\d+S", read_cigars, perl = TRUE) & !is.na(read_cigars)
  has_trail <- grepl("\\d+S$",  read_cigars, perl = TRUE) & !is.na(read_cigars)
  idxs <- which(has_lead | has_trail)

  candidates <- vector("list", 2L * length(idxs))
  cand_idx   <- 1L
  for (i in idxs) {
    if (has_lead[i]) {
      local_bp <- read_starts[i] - genomic_start + 1L
      if (local_bp >= 1L && local_bp <= ref_len) {
        candidates[[cand_idx]] <- list(
          read_name = all_qnames[i], local_breakpoint = local_bp,
          length = 0L, type = "ptd_clip", mapq = all_mapqs[i],
          is_reverse = is_reverse[i], cigar = read_cigars[i],
          read_seq = read_seqs[i]
        )
        cand_idx <- cand_idx + 1L
      }
    }
    if (has_trail[i]) {
      local_bp <- read_ends[i] - genomic_start + 1L
      if (local_bp >= 1L && local_bp <= ref_len) {
        candidates[[cand_idx]] <- list(
          read_name = all_qnames[i], local_breakpoint = local_bp,
          length = 0L, type = "ptd_clip", mapq = all_mapqs[i],
          is_reverse = is_reverse[i], cigar = read_cigars[i],
          read_seq = read_seqs[i]
        )
        cand_idx <- cand_idx + 1L
      }
    }
  }
  candidates <- candidates[seq_len(cand_idx - 1L)]
  if (verbose) message("[TALOS] Found ", length(candidates), " PTD candidates")
  candidates
}


# ----------------------------------------------------------------------------
# Standard k-mer tracing (with logging and auto-PTD conversion for long duplications)
# ----------------------------------------------------------------------------
.extract_candidates_standard <- function(reads, ref_kmers, ptd_mode, min_size,
                                          max_missing_kmers, refine_bp,
                                          use_cigar_bp, genomic_start, ref_dna,
                                          max_itd_length = 1000,
                                          convert_long_to_ptd = TRUE,
                                          verbose = FALSE) {
  n_reads <- length(reads)
  k       <- nchar(ref_kmers[1L])
  if (verbose) message("[TALOS] Extracting standard candidates from ", n_reads, " reads (k=", k, ", ptd_mode=", ptd_mode, ", max_itd_length=", max_itd_length, ")")

  all_qnames  <- .safe_qnames(reads)
  all_mapqs   <- S4Vectors::mcols(reads)$mapq; all_mapqs[is.na(all_mapqs)] <- 0L
  all_flags   <- S4Vectors::mcols(reads)$flag
  all_starts  <- BiocGenerics::start(reads)
  all_ends    <- BiocGenerics::end(reads)
  all_cigars  <- GenomicAlignments::cigar(reads)
  all_seqs    <- as.character(S4Vectors::mcols(reads)$seq)

  candidates      <- vector("list", 2L * n_reads)
  cand_idx        <- 1L
  ref_len         <- nchar(ref_dna)
  max_missing_frac <- max_missing_kmers
  reads_processed <- 0L

  for (i in seq_len(n_reads)) {
    read_seq <- all_seqs[i]

    if (is.na(read_seq)) next
    read_len <- nchar(read_seq)
    if (read_len < k || read_len < (k + min_size)) next

    is_reverse <- bitwAnd(all_flags[i], 0x10L) != 0L

    if (.has_biostrings) {
      dna_obj    <- Biostrings::DNAString(read_seq)
      v          <- Biostrings::Views(dna_obj, start = seq_len(read_len - k + 1L), width = k)
      all_kmers  <- as.character(v)
    } else {
      kmer_starts <- seq_len(read_len - k + 1L)
      all_kmers   <- substring(read_seq, kmer_starts, kmer_starts + k - 1L)
    }

    trace           <- .kmer_match(all_kmers, ref_kmers)
    trace[is.na(trace)] <- -1L

    missing <- sum(trace == -1L)
    if (!ptd_mode && (missing / length(trace) > max_missing_frac)) next

    read_qname <- all_qnames[i]
    read_mapq  <- all_mapqs[i]
    read_cigar <- all_cigars[i]
    read_start <- all_starts[i]
    read_end   <- all_ends[i]

    # ----- Strategy 1: backward jump (direct adjacency) -----
    if (length(trace) >= 2L) {
      for (j in seq(2L, length(trace))) {
        if (trace[j] == -1L || trace[j - 1L] == -1L) next
        if (trace[j - 1L] + k - 1L > trace[j]) {
          pos_after  <- trace[j]
          pos_before <- trace[j - 1L]
          dup_len    <- (pos_before + k - 1L) - pos_after + 1L
          if (dup_len >= min_size) {
            # --- Auto-PTD conversion for long duplications ---
            if (!ptd_mode && dup_len > max_itd_length) {
              if (convert_long_to_ptd) {
                if (verbose) message("[TALOS] Long duplication (", dup_len, " bp) converted to PTD at breakpoint ", genomic_start + pos_after - 1L)
                candidates[[cand_idx]] <- list(
                  read_name = read_qname, local_breakpoint = pos_after,
                  itd_seq = NA_character_, read_seq = read_seq,
                  length = 0L, type = "ptd_clip", mapq = read_mapq,
                  is_reverse = is_reverse, cigar = read_cigar
                )
                cand_idx <- cand_idx + 1L
              }
              next
            }
            dup_seq <- NA_character_
            if (!ptd_mode && dup_len <= 500)
              dup_seq <- substr(read_seq, j, j + dup_len - 1L)
            candidates[[cand_idx]] <- list(
              read_name = read_qname, local_breakpoint = pos_after,
              itd_seq = dup_seq, read_seq = read_seq,
              length = if (ptd_mode) 0L else as.integer(dup_len),
              type = "bwd_jump", mapq = read_mapq,
              is_reverse = is_reverse, cigar = read_cigar)
            cand_idx <- cand_idx + 1L
          }
        }
      }
    }

    # ----- Strategy 2: missing k-mer block -----
    in_break    <- FALSE
    break_start <- NULL
    for (j in seq_along(trace)) {
      if (trace[j] == -1L && !in_break) {
        in_break    <- TRUE
        break_start <- j
      } else if (trace[j] != -1L && in_break) {
        break_end <- j - 1L
        if (break_start > 1L && break_end < length(trace)) {
          pos_before <- trace[break_start - 1L]
          pos_after  <- trace[break_end + 1L]
          if (pos_before != -1L && pos_after != -1L &&
              (pos_before + k - 1L > pos_after)) {
            dup_len <- (pos_before + k - 1L) - pos_after + 1L
            if (dup_len >= min_size) {
              # --- Auto-PTD conversion for long duplications ---
              if (!ptd_mode && dup_len > max_itd_length) {
                if (convert_long_to_ptd) {
                  if (verbose) message("[TALOS] Long duplication (", dup_len, " bp) converted to PTD at breakpoint ", genomic_start + pos_after - 1L)
                  candidates[[cand_idx]] <- list(
                    read_name = read_qname, local_breakpoint = pos_after,
                    itd_seq = NA_character_, read_seq = read_seq,
                    length = 0L, type = "ptd_clip", mapq = read_mapq,
                    is_reverse = is_reverse, cigar = read_cigar
                  )
                  cand_idx <- cand_idx + 1L
                }
                next
              }
              dup_seq <- NA_character_
              if (!ptd_mode && dup_len <= 500)
                dup_seq <- substr(read_seq, break_start, break_start + dup_len - 1L)
              candidates[[cand_idx]] <- list(
                read_name = read_qname, local_breakpoint = pos_after,
                itd_seq = dup_seq, read_seq = read_seq,
                length = if (ptd_mode) 0L else as.integer(dup_len),
                type = "missing_kmer", mapq = read_mapq,
                is_reverse = is_reverse, cigar = read_cigar)
              cand_idx <- cand_idx + 1L
            }
          }
        } else if (ptd_mode && break_start == 1L) {
          pos_after <- trace[j]
          if (pos_after != -1L) {
            refined_bp <- pos_after
            if (use_cigar_bp && grepl("^\\d+S", read_cigar, perl = TRUE)) {
              refined_bp <- read_start - genomic_start + 1L
            } else if (refine_bp) {
              max_offset <- min(40L, max(0, pos_after - 1L))
              for (offset in 0:max_offset) {
                test_pos <- pos_after - offset
                if (test_pos + k - 1L <= ref_len) {
                  ref_kmer  <- substr(ref_dna, test_pos, test_pos + k - 1L)
                  read_kmer <- substr(read_seq, j, j + k - 1L)
                  if (ref_kmer == read_kmer) {
                    match_len <- 0L
                    while ((j - match_len > 0L) && (test_pos - match_len > 0L) &&
                           substr(read_seq, j - match_len, j - match_len) ==
                           substr(ref_dna, test_pos - match_len, test_pos - match_len)) {
                      match_len <- match_len + 1L
                    }
                    refined_bp <- test_pos - match_len + 1L
                    break
                  }
                }
              }
            }
            refined_bp <- max(1L, refined_bp)
            candidates[[cand_idx]] <- list(
              read_name = read_qname, local_breakpoint = refined_bp,
              itd_seq = NA_character_, read_seq = read_seq,
              length = 0L,
              type = "ptd_clip", mapq = read_mapq,
              is_reverse = is_reverse, cigar = read_cigar)
            cand_idx <- cand_idx + 1L
          }
        }
        in_break <- FALSE
      }
    }

    if (in_break && ptd_mode && break_start > 1L) {
      pos_before <- trace[break_start - 1L]
      if (pos_before != -1L) {
        refined_bp <- pos_before + k - 1L
        if (use_cigar_bp && grepl("\\d+S$", read_cigar, perl = TRUE))
          refined_bp <- read_end - genomic_start + 1L
        refined_bp <- max(1L, refined_bp)
        candidates[[cand_idx]] <- list(
          read_name = read_qname, local_breakpoint = refined_bp,
          itd_seq = NA_character_, read_seq = read_seq,
          length = 0L,
          type = "ptd_clip", mapq = read_mapq,
          is_reverse = is_reverse, cigar = read_cigar)
        cand_idx <- cand_idx + 1L
      }
    }
    reads_processed <- reads_processed + 1L
    if (verbose && reads_processed %% 1000 == 0) message("[TALOS] Processed ", reads_processed, " / ", n_reads, " reads for candidate extraction")
  }
  candidates <- candidates[seq_len(cand_idx - 1L)]
  if (verbose) message("[TALOS] Found ", length(candidates), " standard candidates")
  candidates
}


# ----------------------------------------------------------------------------
# Candidate list → data.frame
# ----------------------------------------------------------------------------
.candidates_to_df <- function(candidates, genomic_start) {
  n <- length(candidates)
  if (n == 0L) return(data.frame())
  
  data.frame(
    breakpoint = vapply(candidates, function(x)
      as.integer(as.numeric(genomic_start) + as.numeric(x$local_breakpoint) - 1L),
      integer(1L)),
    length     = vapply(candidates, function(x) {
      v <- x$length; if (is.null(v) || (length(v) == 1L && is.na(v))) NA_integer_ else as.integer(v[1L])
    }, integer(1L)),
    read_name  = vapply(candidates, function(x) x$read_name,  character(1L)),
    read_seq   = vapply(candidates, function(x) {
      s <- x$read_seq; if (is.null(s) || is.na(s)) NA_character_ else s
    }, character(1L)),
    type       = vapply(candidates, function(x) x$type,       character(1L)),
    mapq       = vapply(candidates, function(x) {
      v <- x$mapq; if (is.null(v) || (length(v) == 1L && is.na(v))) NA_integer_ else as.integer(v[1L])
    }, integer(1L)),
    is_reverse = vapply(candidates, function(x) x$is_reverse, logical(1L)),
    cigar      = vapply(candidates, function(x) x$cigar,      character(1L)),
    itd_seq    = vapply(candidates, function(x) {
      s <- x$itd_seq; if (is.null(s) || is.na(s)) NA_character_ else s
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}

.cluster_breakpoints <- function(bp_values, cluster_tolerance) {
  if (length(bp_values) == 0) return(list())
  bp_sorted <- sort(bp_values)
  if (length(bp_sorted) == 1L) return(list(bp_sorted))

  clusters          <- list()
  cluster_start_idx <- 1L
  for (j in seq(2L, length(bp_sorted))) {
    if (bp_sorted[j] - bp_sorted[j - 1L] > cluster_tolerance) {
      clusters[[length(clusters) + 1L]] <- bp_sorted[cluster_start_idx:(j - 1L)]
      cluster_start_idx <- j
    }
  }
  clusters[[length(clusters) + 1L]] <- bp_sorted[cluster_start_idx:length(bp_sorted)]
  clusters
}


.prepare_wildtype_info <- function(all_reads, genomic_start, genomic_end) {
  raw_qnames <- .safe_qnames(all_reads)
  all_flags  <- S4Vectors::mcols(all_reads)$flag
  is_primary <- bitwAnd(all_flags, 0x100L) == 0L & bitwAnd(all_flags, 0x800L) == 0L
  all_reads_gr <- GenomicRanges::granges(all_reads)
  list(qnames = raw_qnames, is_primary = is_primary, gr = all_reads_gr,
       region_start = genomic_start, region_end = genomic_end)
}


.compute_left_right_coverage <- function(wildtype_info, len_specific_bp,
                                          support_qnames, buffer = 10) {
  strict_idx <- which(wildtype_info$is_primary &
                      start(wildtype_info$gr) <= (len_specific_bp - buffer) &
                      end(wildtype_info$gr)   >= (len_specific_bp + buffer))
  strict_qnames <- unique(wildtype_info$qnames[strict_idx])
  pure_wt       <- setdiff(strict_qnames, support_qnames)
  wt_count      <- length(pure_wt)
  list(left = wt_count, right = wt_count, avg = wt_count)
}


.calculate_vaf_and_depth <- function(raw_support, wildtype_info, len_specific_bp,
                                      support_qnames, best_len,
                                      nominal_read_len, max_correction, buffer = 10) {
  lr        <- .compute_left_right_coverage(wildtype_info, len_specific_bp,
                                            support_qnames, buffer)
  strict_wt <- lr$avg
  len_used  <- if (is.na(best_len)) 0L else best_len
  size_scalar       <- min(1 + (0.25 * len_used / nominal_read_len), max_correction)
  corrected_support <- round(raw_support * size_scalar)
  depth_at_bp       <- corrected_support + strict_wt
  corrected_af      <- if (depth_at_bp == 0) 0 else corrected_support / depth_at_bp

  list(corrected_support = corrected_support, depth_at_bp = depth_at_bp,
       raw_af = corrected_af, raw_depth = raw_support + strict_wt,
       left_wt = strict_wt, right_wt = strict_wt, avg_wt = strict_wt)
}


# ----------------------------------------------------------------------------
# Compute variant metrics – ENHANCED with fallback coverage and corrected filter
# ----------------------------------------------------------------------------
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
                                      debug = FALSE, verbose = FALSE) {

  if (verbose) message("[TALOS] Computing metrics for cluster at ", paste(range(cluster_bps), collapse="-"), " length=", best_len)

  # ---- Initialize all metric variables to avoid missing values ----
  consistency_score    <- NA_real_
  itd_coverage_percent <- NA_real_
  itd_coverage_rle     <- NA_character_
  ref_match_observed   <- NA_real_
  ref_match_total      <- NA_real_
  alignment_score_total<- NA_real_
  total_support_bases  <- NA_integer_
  coverage_drop        <- NA_real_
  median_microhomology <- NA_real_
  repeat_entropy       <- NA_real_
  discordant_ratio     <- NA_real_
  orientation          <- NA_character_
  hgvs                 <- list(c_notation = NA_character_, p_notation = NA_character_)

  len_specific_bp <- round(median(cluster_bps))
  support_qnames  <- unique(support_rows$read_name)
  raw_support     <- length(support_qnames)

  # ---- Compute softclip side counts ----
  left_sc_count  <- 0L
  right_sc_count <- 0L
  both_sc_count  <- 0L
  if (nrow(support_rows) > 0 && !is.null(support_rows$cigar)) {
    for (cig in support_rows$cigar) {
      if (is.na(cig)) next
      has_left  <- grepl("^\\d+S", cig, perl = TRUE)
      has_right <- grepl("\\d+S$", cig, perl = TRUE)
      if (has_left && has_right) {
        both_sc_count <- both_sc_count + 1L
        left_sc_count <- left_sc_count + 1L
        right_sc_count <- right_sc_count + 1L
      } else if (has_left) {
        left_sc_count <- left_sc_count + 1L
      } else if (has_right) {
        right_sc_count <- right_sc_count + 1L
      }
    }
  }

  # ---- ITD sequence extraction (with observed part tracking) ----
  itd_seq          <- NA_character_
  imputed          <- TRUE
  sequence_partial <- FALSE
  observed_seq     <- NA_character_   # bases actually observed in reads
  observed_len     <- 0L

  valid_seqs <- support_rows$itd_seq[!is.na(support_rows$itd_seq) &
                                      nchar(support_rows$itd_seq) > 0]
  if (length(valid_seqs) > 0) {
    itd_seq <- names(sort(table(valid_seqs), decreasing = TRUE))[1L]
    imputed <- FALSE
    observed_seq <- itd_seq
    observed_len <- nchar(itd_seq)
  }

  if (is.na(itd_seq)) {
    extracted_seqs <- vector("character", 2L * nrow(support_rows))
    seq_idx <- 1L
    for (i in seq_len(nrow(support_rows))) {
      cig   <- support_rows$cigar[i]
      r_seq <- support_rows$read_seq[i]
      sc    <- .get_softclips(cig, r_seq)
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
    if (length(extracted_seqs) > 0) {
      itd_seq <- names(sort(table(extracted_seqs), decreasing = TRUE))[1L]
      imputed <- FALSE
      observed_seq <- itd_seq
      observed_len <- nchar(itd_seq)
    }
  }

  if (!is.na(itd_seq) && !is.na(best_len) && best_len > 0 && nchar(itd_seq) > best_len) {
    itd_seq <- substr(itd_seq, 1L, best_len)
    if (!is.na(observed_seq) && nchar(observed_seq) > best_len)
      observed_seq <- substr(observed_seq, 1L, best_len)
    observed_len <- min(observed_len, best_len)
  }

  if (!is.na(itd_seq) && !is.na(best_len) && best_len > 0 &&
      nchar(itd_seq) < best_len && !ptd_mode) {
    partial_obs_len <- nchar(itd_seq)
    ref_local_start <- len_specific_bp - genomic_start + partial_obs_len + 1L
    ref_local_end   <- ref_local_start + (best_len - partial_obs_len) - 1L
    if (ref_local_end <= nchar(ref_dna)) {
      itd_seq          <- paste0(itd_seq, substr(ref_dna, ref_local_start, ref_local_end))
      sequence_partial <- TRUE
      imputed          <- TRUE
      observed_len     <- partial_obs_len
    }
  }

  if (is.na(itd_seq) && !is.na(best_len) && best_len > 0 && !ptd_mode) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end   <- local_start + best_len - 1L
    if (local_end <= nchar(ref_dna)) {
      itd_seq <- substr(ref_dna, local_start, local_end)
      imputed <- TRUE
      sequence_partial <- FALSE
      observed_seq <- NA_character_
      observed_len <- 0L
    }
  }

  # ---- Consistency / ITD coverage (with reproducibility seed) ----
  if ((do_consistency || do_itd_coverage) &&
      !is.na(itd_seq) && nchar(itd_seq) > 0 && nrow(support_rows) > 0) {

    valid_mask <- !is.na(support_rows$read_seq) & nchar(support_rows$read_seq) > 0
    v_seqs     <- support_rows$read_seq[valid_mask]

    if (length(v_seqs) > 0) {
      if (length(v_seqs) > max_pairwise_alignments) {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .GlobalEnv$.Random.seed else NULL
        set.seed(42)
        v_seqs <- v_seqs[sample(length(v_seqs), max_pairwise_alignments)]
        if (!is.null(old_seed)) .GlobalEnv$.Random.seed <- old_seed else rm(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      }

      if (.has_pwalign && .has_biostrings) {
        aln  <- pwalign::pairwiseAlignment(
          pattern = Biostrings::DNAStringSet(v_seqs),
          subject = Biostrings::DNAString(itd_seq),
          type    = "local"
        )
        pids <- pwalign::pid(aln)

        if (do_consistency) {
          matches           <- sum(pids >= 90, na.rm = TRUE)
          consistency_score <- matches / length(v_seqs) * 100
        }

        if (do_itd_coverage) {
          good_aln_idx <- which(pids >= 90)
          if (length(good_aln_idx) > 0) {
            subj_ranges <- pwalign::subject(aln[good_aln_idx])
            starts <- start(subj_ranges); ends <- end(subj_ranges)
            cov_vec <- integer(nchar(itd_seq))
            for (ki in seq_along(starts))
              cov_vec[starts[ki]:ends[ki]] <- cov_vec[starts[ki]:ends[ki]] + 1L
            itd_coverage_percent <- sum(cov_vec > 0) / length(cov_vec) * 100
            rle_res              <- rle(cov_vec)
            itd_coverage_rle     <- paste(
              paste0(rle_res$values, ":", rle_res$lengths), collapse = ",")
          } else {
            itd_coverage_percent <- 0
            itd_coverage_rle     <- paste0("0:", nchar(itd_seq))
          }
        }
      } else {
        if (do_consistency) {
          matches           <- sum(grepl(itd_seq, v_seqs, fixed = TRUE))
          consistency_score <- matches / length(v_seqs) * 100
        }
        if (do_itd_coverage) {
          itd_coverage_percent <- if (any(grepl(itd_seq, v_seqs, fixed = TRUE))) 100 else 0
          itd_coverage_rle     <- paste0(as.integer(itd_coverage_percent > 0),
                                         ":", nchar(itd_seq))
        }
      }
    }
  }

  # ---- Alignment scores (observed and total) ----
  if (do_alignment_score && !is.na(itd_seq) && !is.na(best_len) && best_len > 0) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end   <- local_start + best_len - 1L
    if (local_end <= nchar(ref_dna)) {
      ref_segment <- substr(ref_dna, local_start, local_end)
      min_len_total <- min(nchar(itd_seq), nchar(ref_segment))
      if (min_len_total > 0) {
        if (itd_seq == ref_segment) {
          alignment_score_total <- 1.0
        } else {
          itd_raw <- charToRaw(substr(itd_seq, 1L, min_len_total))
          ref_raw <- charToRaw(substr(ref_segment, 1L, min_len_total))
          alignment_score_total <- sum(itd_raw == ref_raw) / min_len_total
        }
        ref_match_total <- round(alignment_score_total * 100, 1)

        if (!is.na(observed_seq) && observed_len > 0) {
          obs_len <- min(observed_len, min_len_total)
          if (obs_len > 0) {
            obs_sub <- substr(observed_seq, 1L, obs_len)
            ref_sub <- substr(ref_segment, 1L, obs_len)
            if (obs_sub == ref_sub) {
              alignment_score_obs <- 1.0
            } else {
              obs_raw <- charToRaw(obs_sub)
              ref_raw_obs <- charToRaw(ref_sub)
              alignment_score_obs <- sum(obs_raw == ref_raw_obs) / obs_len
            }
            ref_match_observed <- round(alignment_score_obs * 100, 1)
          }
        }
      }
    }
  }

  # ---- Support bases ----
  if (do_support_bases && nrow(support_rows) > 0) {
    cigars  <- support_rows$cigar; cigars[is.na(cigars)] <- ""
    lead_S  <- as.integer(regmatches(cigars, regexec("^(\\d+)S", cigars)) |>
                          sapply(function(m) if (length(m) > 1) m[2] else "0"))
    trail_S <- as.integer(regmatches(cigars, regexec("(\\d+)S$", cigars)) |>
                          sapply(function(m) if (length(m) > 1) m[2] else "0"))
    total_support_bases <- sum(lead_S, trail_S, na.rm = TRUE)
  }

  # ---- VAF / depth ----
  if (raw_support < min_support) return(NULL)

  vaf_metrics       <- .calculate_vaf_and_depth(
    raw_support, wildtype_info, len_specific_bp, support_qnames,
    best_len, nominal_read_len, max_correction, buffer = 10
  )
  corrected_support <- vaf_metrics$corrected_support
  depth_at_bp       <- vaf_metrics$depth_at_bp
  raw_af            <- vaf_metrics$raw_af
  wildtype_reads    <- round(vaf_metrics$avg_wt)

  if (wildtype_reads < min_wt_reads) return(NULL)
  if (raw_af < vaf_threshold)        return(NULL)

  # ---- Wildtype-relative softclip percentages ----
  left_sc_pct_wt  <- if (wildtype_reads > 0) (left_sc_count / wildtype_reads) * 100 else NA_real_
  right_sc_pct_wt <- if (wildtype_reads > 0) (right_sc_count / wildtype_reads) * 100 else NA_real_

  # ---- Diagnostic metrics ----
  strand_bias       <- if (nrow(support_rows) > 0)
    round(mean(support_rows$is_reverse, na.rm = TRUE), 4) else NA_real_
  mean_support_mapq <- if (nrow(support_rows) > 0)
    round(mean(support_rows$mapq, na.rm = TRUE), 1) else NA_real_
  bp_spread         <- if (nrow(support_rows) > 1)
    as.integer(max(support_rows$breakpoint) - min(support_rows$breakpoint)) else 0L
  softclip_frac     <- if (nrow(support_rows) > 0) {
    round(mean(grepl("S", support_rows$cigar, fixed = TRUE), na.rm = TRUE), 4)
  } else NA_real_
  unique_bps        <- length(unique(support_rows$breakpoint))

  # ---- Percentage of softclip reads relative to support ----
  left_sc_pct_support <- if (corrected_support > 0) (left_sc_count / corrected_support) * 100 else NA_real_
  right_sc_pct_support <- if (corrected_support > 0) (right_sc_count / corrected_support) * 100 else NA_real_

  # ---- Coverage drop ----
  if (do_coverage_drop && !is.na(best_len) && best_len > 0) {
    cov_ok <- FALSE
    if (!is.null(all_reads_cov)) {
      coverage_drop <- compute_coverage_drop(all_reads_cov, gene_config$chrom,
                                             len_specific_bp, flank = 200,
                                             verbose = debug)
      if (!is.na(coverage_drop)) cov_ok <- TRUE
    }
    if (!cov_ok && !is.null(wildtype_info$gr) && length(wildtype_info$gr) > 0) {
      if (debug) message("Coverage drop: recomputing from wildtype reads")
      reg_gr <- GenomicRanges::GRanges(
        seqnames = gene_config$chrom,
        ranges = IRanges::IRanges(start = gene_config$genomic_start,
                                  end   = gene_config$genomic_end)
      )
      ov <- GenomicRanges::findOverlaps(wildtype_info$gr, reg_gr)
      if (length(ov) > 0) {
        reads_sub <- wildtype_info$gr[unique(S4Vectors::queryHits(ov))]
        cov_fallback <- GenomicAlignments::coverage(reads_sub)[[gene_config$chrom]]
        if (!is.null(cov_fallback)) {
          coverage_drop <- compute_coverage_drop(cov_fallback, gene_config$chrom,
                                                 len_specific_bp, flank = 200,
                                                 verbose = debug)
        }
      }
    }
    if (debug && is.na(coverage_drop)) message("Coverage drop: still NA after fallback")
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

  # ---- Orientation ----
  if (do_detect_orientation && !is.na(itd_seq) && !is.na(best_len) && best_len > 0) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end   <- local_start + best_len - 1L
    if (local_end <= nchar(ref_dna))
      orientation <- detect_orientation(itd_seq, substr(ref_dna, local_start, local_end))
  }

  # ---- HGVS ----
  is_exonic <- FALSE
  if (!is.null(gene_config$exons)) {
    bp_gr <- GenomicRanges::GRanges(
      seqnames = gene_config$chrom,
      ranges   = IRanges::IRanges(start = len_specific_bp, end = len_specific_bp)
    )
    is_exonic <- length(GenomicRanges::findOverlaps(bp_gr, gene_config$exons,
                                                    type = "within")) > 0
  }

  if (do_hgvs) {
    dup_len_for_hgvs <- if (!is.na(itd_seq) && nchar(itd_seq) > 0)
      nchar(itd_seq)
    else if (!is.na(best_len)) best_len
    else NA_integer_

    if (!is.null(gene_config$all_exons) || !is.null(gene_config$exons)) {
      hgvs <- compute_hgvs_annotations(gene_config, len_specific_bp,
                                       itd_seq, dup_len_for_hgvs, debug)
      if (!is_exonic && !is.na(hgvs$p_notation)) {
        hgvs$p_notation <- "intronic"
        if (is.na(hgvs$c_notation)) hgvs$c_notation <- "c.?+?"
      }
    }
  }

  # ---- Return list (all variables now guaranteed to exist) ----
  result <- list(
    Sample = wildtype_info$sample_name, Gene = gene_config$gene,
    Genome = if (is.null(gene_config$build)) "custom" else gene_config$build,
    GenomicPosition = len_specific_bp, ITD_Sequence = itd_seq, Length = best_len,
    SupportingReads = corrected_support, WildtypeReads = wildtype_reads,
    DepthAtBreakpoint = depth_at_bp, AlleleFrequency = raw_af,
    HGVS_cDNA = hgvs$c_notation,
    HGVS_Protein = hgvs$p_notation,
    StrandBias = strand_bias, MeanSupportMAPQ = mean_support_mapq,
    BreakpointSpread = bp_spread, SoftclipFraction = softclip_frac,
    UniqueBreakpoints = unique_bps, CoverageDrop = coverage_drop,
    MedianMicrohomology = median_microhomology, DiscordantRatio = discordant_ratio,
    RepeatEntropy = repeat_entropy,
    SequenceImputed   = imputed, SequencePartial = sequence_partial,
    SequenceSource    = if (!imputed) "observed"
                        else if (sequence_partial) "partial_read+ref"
                        else "ref_imputed",
    RefMatch_Observed = ref_match_observed,
    RefMatch_Total    = ref_match_total,
    ITDReadCoverage   = itd_coverage_percent,
    ITDCoverageRLE    = itd_coverage_rle,
    SupportConsistency = consistency_score,
    AlignmentScore    = alignment_score_total,   
    TotalSupportBases = total_support_bases,
    Orientation       = orientation,
    LeftSoftclipCount  = left_sc_count,
    RightSoftclipCount = right_sc_count,
    BothSoftclipCount  = both_sc_count,
    LeftSoftclipPctSupport  = left_sc_pct_support,
    RightSoftclipPctSupport = right_sc_pct_support,
    LeftSoftclipPctWT       = left_sc_pct_wt,
    RightSoftclipPctWT      = right_sc_pct_wt
  )
  if (verbose) message("[TALOS] Metrics computed: support=", corrected_support, ", AF=", round(raw_af,4), ", length=", best_len)
  result
}


# ============================================================================
# MODIFIED .apply_filters – added min_length / max_length checks,
#               absolute minimum softclip counts for ITDs,
#               AND minimum wildtype softclip percentages for ITDs
# ============================================================================
.apply_filters <- function(metrics, thresholds, min_length = NULL, max_length = NULL) { 
  with(thresholds, {
    # ---- Length filters ----
    if (!is.null(min_length) && !is.na(metrics$Length) && metrics$Length < min_length) return(FALSE)
    if (!is.null(max_length) && !is.na(metrics$Length) && metrics$Length > max_length) return(FALSE)

    if (!is.na(metrics$CoverageDrop) &&
        metrics$CoverageDrop < min_coverage_drop)       return(FALSE)
    if (!is.na(metrics$MedianMicrohomology) &&
        metrics$MedianMicrohomology < min_microhomology) return(FALSE)
    if (!is.na(metrics$DiscordantRatio) &&
        metrics$DiscordantRatio < min_discordant_ratio)  return(FALSE)
    if (!is.na(metrics$RepeatEntropy) &&
        metrics$RepeatEntropy < min_entropy)             return(FALSE)
    if (!is.na(metrics$StrandBias) &&
        (metrics$StrandBias < min_strand_bias ||
         metrics$StrandBias > max_strand_bias))          return(FALSE)
    if (!is.na(metrics$MeanSupportMAPQ) &&
        metrics$MeanSupportMAPQ < min_mean_support_mapq) return(FALSE)
    if (!is.na(metrics$BreakpointSpread) &&
        metrics$BreakpointSpread > max_breakpoint_spread) return(FALSE)
    if (!is.na(metrics$SoftclipFraction) &&
        metrics$SoftclipFraction < min_softclip_fraction) return(FALSE)
    if (!is.na(metrics$UniqueBreakpoints) &&
        metrics$UniqueBreakpoints < min_unique_breakpoints) return(FALSE)
    if (!is.na(metrics$ITDReadCoverage) &&
        metrics$ITDReadCoverage < min_itd_read_coverage) return(FALSE)
    if (!is.na(metrics$SequenceImputed) && isTRUE(metrics$SequenceImputed) &&
        !isTRUE(metrics$SequencePartial) &&
        !is.na(metrics$Length) && metrics$Length > 40)   return(FALSE)

    # ---- Soft‑clip filters (different for ITDs vs PTDs) ----
    left  <- metrics$LeftSoftclipCount
    right <- metrics$RightSoftclipCount
    if (!is.na(left) && !is.na(right)) {
      # ITD case: length > 0
      if (!is.na(metrics$Length) && metrics$Length > 0L) {
        # 1) Absolute minimum count on each side
        if (left < min_side_softclip_reads || right < min_side_softclip_reads) {
          return(FALSE)
        }
        # 2) Minimum wildtype softclip percentage on each side
        left_wt_pct  <- metrics$LeftSoftclipPctWT
        right_wt_pct <- metrics$RightSoftclipPctWT
        if (!is.na(left_wt_pct) && left_wt_pct < min_left_softclip_pct_wt) return(FALSE)
        if (!is.na(right_wt_pct) && right_wt_pct < min_right_softclip_pct_wt) return(FALSE)
      } else {
        # PTD case (Length == 0): use ratio + absolute filter (unchanged)
        if (left < min_side_softclip_reads || right < min_side_softclip_reads) {
          if (max(left, right) / min(left, right) > max_side_ratio) {
            return(FALSE)
          }
        }
        # PTD absolute side filter (min_abs_side_softclip)
        if (left < min_abs_side_softclip || right < min_abs_side_softclip) {
          return(FALSE)
        }
        # PTD soft‑clip percentage filters (unchanged)
        left_pct  <- metrics$LeftSoftclipPctSupport
        right_pct <- metrics$RightSoftclipPctSupport
        if (!is.na(left_pct) && left_pct < min_softclip_pct_side) return(FALSE)
        if (!is.na(right_pct) && right_pct < min_softclip_pct_side) return(FALSE)
        left_wt_pct  <- metrics$LeftSoftclipPctWT
        right_wt_pct <- metrics$RightSoftclipPctWT
        if (!is.na(left_wt_pct) && left_wt_pct < min_left_softclip_pct_wt) return(FALSE)
        if (!is.na(right_wt_pct) && right_wt_pct < min_right_softclip_pct_wt) return(FALSE)
      }
    }

    TRUE
  })
}
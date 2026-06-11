# ============================================================================
# TALOS – Candidate extraction and breakpoint clustering
# ============================================================================

# ---------------------------------------------------------------------------
# Safe query-name extraction
# ---------------------------------------------------------------------------
.safe_qnames <- function(reads) {
  qn <- as.character(S4Vectors::mcols(reads)$qname)
  if (length(qn) == 0L || all(is.na(qn))) qn <- names(reads)
  na_idx <- which(is.na(qn))
  if (length(na_idx) > 0L) qn[na_idx] <- paste0("read_", na_idx)
  qn
}


# ---------------------------------------------------------------------------
# PTD fast path – pure CIGAR soft-clip extraction (no k-mer analysis)
# ---------------------------------------------------------------------------
.extract_candidates_ptd <- function(reads, genomic_start, ref_len,
                                     verbose = FALSE) {
  n_reads <- length(reads)
  if (n_reads == 0L) return(list())
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


# ---------------------------------------------------------------------------
# Standard k-mer tracing (both backward-jump and missing-block strategies)
# with auto-PTD conversion for duplications longer than max_itd_length.
# ---------------------------------------------------------------------------
.extract_candidates_standard <- function(reads, ref_kmers, ptd_mode, min_size,
                                          max_missing_kmers, refine_bp,
                                          use_cigar_bp, genomic_start, ref_dna,
                                          max_itd_length = 1000L,
                                          convert_long_to_ptd = TRUE,
                                          verbose = FALSE) {
  n_reads <- length(reads)
  k       <- nchar(ref_kmers[1L])
  if (verbose) message(
    "[TALOS] Extracting standard candidates from ", n_reads, " reads",
    " (k=", k, ", ptd_mode=", ptd_mode, ", max_itd_length=", max_itd_length, ")"
  )

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
  reads_processed <- 0L

  for (i in seq_len(n_reads)) {
    read_seq <- all_seqs[i]
    if (is.na(read_seq)) next
    read_len <- nchar(read_seq)
    if (read_len < k || read_len < (k + min_size)) next

    is_reverse <- bitwAnd(all_flags[i], 0x10L) != 0L

    if (.has_biostrings) {
      dna_obj   <- Biostrings::DNAString(read_seq)
      v         <- Biostrings::Views(dna_obj, start = seq_len(read_len - k + 1L), width = k)
      all_kmers <- as.character(v)
    } else {
      kmer_starts <- seq_len(read_len - k + 1L)
      all_kmers   <- substring(read_seq, kmer_starts, kmer_starts + k - 1L)
    }

    trace              <- .kmer_match(all_kmers, ref_kmers)
    trace[is.na(trace)] <- -1L

    missing <- sum(trace == -1L)
    if (!ptd_mode && (missing / length(trace) > max_missing_kmers)) next

    read_qname <- all_qnames[i]; read_mapq  <- all_mapqs[i]
    read_cigar <- all_cigars[i]; read_start <- all_starts[i]; read_end <- all_ends[i]

    # --- Strategy 1: backward jump (direct adjacency) ---
    if (length(trace) >= 2L) {
      for (j in seq(2L, length(trace))) {
        if (trace[j] == -1L || trace[j - 1L] == -1L) next
        if (trace[j - 1L] + k - 1L > trace[j]) {
          pos_after  <- trace[j]
          pos_before <- trace[j - 1L]
          dup_len    <- (pos_before + k - 1L) - pos_after + 1L
          if (dup_len >= min_size) {
            if (!ptd_mode && dup_len > max_itd_length) {
              if (convert_long_to_ptd) {
                if (verbose) message("[TALOS] Long dup (", dup_len, " bp) → PTD at ",
                                     genomic_start + pos_after - 1L)
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
            if (!ptd_mode && dup_len <= 500L)
              dup_seq <- substr(read_seq, j, j + dup_len - 1L)
            candidates[[cand_idx]] <- list(
              read_name = read_qname, local_breakpoint = pos_after,
              itd_seq = dup_seq, read_seq = read_seq,
              length = if (ptd_mode) 0L else as.integer(dup_len),
              type = "bwd_jump", mapq = read_mapq,
              is_reverse = is_reverse, cigar = read_cigar
            )
            cand_idx <- cand_idx + 1L
          }
        }
      }
    }

    # --- Strategy 2: missing k-mer block ---
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
              if (!ptd_mode && dup_len > max_itd_length) {
                if (convert_long_to_ptd) {
                  if (verbose) message("[TALOS] Long dup (", dup_len, " bp) → PTD at ",
                                       genomic_start + pos_after - 1L)
                  candidates[[cand_idx]] <- list(
                    read_name = read_qname, local_breakpoint = pos_after,
                    itd_seq = NA_character_, read_seq = read_seq,
                    length = 0L, type = "ptd_clip", mapq = read_mapq,
                    is_reverse = is_reverse, cigar = read_cigar
                  )
                  cand_idx <- cand_idx + 1L
                }
                in_break <- FALSE
                next
              }
              dup_seq <- NA_character_
              if (!ptd_mode && dup_len <= 500L)
                dup_seq <- substr(read_seq, break_start, break_start + dup_len - 1L)
              candidates[[cand_idx]] <- list(
                read_name = read_qname, local_breakpoint = pos_after,
                itd_seq = dup_seq, read_seq = read_seq,
                length = if (ptd_mode) 0L else as.integer(dup_len),
                type = "missing_kmer", mapq = read_mapq,
                is_reverse = is_reverse, cigar = read_cigar
              )
              cand_idx <- cand_idx + 1L
            }
          }
        } else if (ptd_mode && break_start == 1L) {
          pos_after <- trace[j]
          if (pos_after != -1L) {
            refined_bp <- pos_after
            if (use_cigar_bp && grepl("^\\d+S", read_cigar, perl = TRUE))
              refined_bp <- read_start - genomic_start + 1L
            else if (refine_bp) {
              max_offset <- min(40L, max(0L, pos_after - 1L))
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
              length = 0L, type = "ptd_clip", mapq = read_mapq,
              is_reverse = is_reverse, cigar = read_cigar
            )
            cand_idx <- cand_idx + 1L
          }
        }
        in_break <- FALSE
      }
    }

    # Trailing missing block (PTD trailing soft-clip)
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
          length = 0L, type = "ptd_clip", mapq = read_mapq,
          is_reverse = is_reverse, cigar = read_cigar
        )
        cand_idx <- cand_idx + 1L
      }
    }

    reads_processed <- reads_processed + 1L
    if (verbose && reads_processed %% 1000L == 0L)
      message("[TALOS] Processed ", reads_processed, " / ", n_reads,
              " reads for candidate extraction")
  }

  candidates <- candidates[seq_len(cand_idx - 1L)]
  if (verbose) message("[TALOS] Found ", length(candidates), " standard candidates")
  candidates
}


# ---------------------------------------------------------------------------
# Candidate list → data.frame
# ---------------------------------------------------------------------------
.candidates_to_df <- function(candidates, genomic_start) {
  n <- length(candidates)
  if (n == 0L) return(data.frame())
  data.frame(
    breakpoint = vapply(candidates, function(x)
      as.integer(as.numeric(genomic_start) + as.numeric(x$local_breakpoint) - 1L),
      integer(1L)),
    length     = vapply(candidates, function(x) {
      v <- x$length
      if (is.null(v) || (length(v) == 1L && is.na(v))) NA_integer_ else as.integer(v[1L])
    }, integer(1L)),
    read_name  = vapply(candidates, `[[`, character(1L), "read_name"),
    read_seq   = vapply(candidates, function(x) {
      s <- x$read_seq; if (is.null(s) || is.na(s)) NA_character_ else s
    }, character(1L)),
    type       = vapply(candidates, `[[`, character(1L), "type"),
    mapq       = vapply(candidates, function(x) {
      v <- x$mapq
      if (is.null(v) || (length(v) == 1L && is.na(v))) NA_integer_ else as.integer(v[1L])
    }, integer(1L)),
    is_reverse = vapply(candidates, `[[`, logical(1L), "is_reverse"),
    cigar      = vapply(candidates, `[[`, character(1L), "cigar"),
    itd_seq    = vapply(candidates, function(x) {
      s <- x$itd_seq; if (is.null(s) || is.na(s)) NA_character_ else s
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Cluster nearby breakpoints (within cluster_tolerance bp)
# ---------------------------------------------------------------------------
.cluster_breakpoints <- function(bp_values, cluster_tolerance) {
  if (length(bp_values) == 0L) return(list())
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


# ---------------------------------------------------------------------------
# Prepare wildtype-read info struct (used downstream for VAF / depth)
# ---------------------------------------------------------------------------
.prepare_wildtype_info <- function(all_reads, genomic_start, genomic_end) {
  raw_qnames <- .safe_qnames(all_reads)
  all_flags  <- S4Vectors::mcols(all_reads)$flag
  is_primary <- bitwAnd(all_flags, 0x100L) == 0L & bitwAnd(all_flags, 0x800L) == 0L
  list(
    qnames       = raw_qnames,
    is_primary   = is_primary,
    gr           = GenomicRanges::granges(all_reads),
    region_start = genomic_start,
    region_end   = genomic_end
  )
}
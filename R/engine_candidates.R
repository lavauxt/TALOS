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
# Assemble PTD consensus from soft-clipped sequences (new for PTD improvement)
# ---------------------------------------------------------------------------
.assemble_ptd_consensus <- function(support_rows, min_reads = 3L, min_len = 15L) {
  # Collect all soft-clip sequences (lead and trail)
  soft_seqs <- character()
  for (i in seq_len(nrow(support_rows))) {
    sc <- .get_softclips(support_rows$cigar[i], support_rows$read_seq[i])
    if (!is.na(sc$lead)) soft_seqs <- c(soft_seqs, sc$lead)
    if (!is.na(sc$trail)) soft_seqs <- c(soft_seqs, sc$trail)
  }
  if (length(soft_seqs) < min_reads) return(NA_character_)

  # Keep only sequences with length >= min_len
  soft_seqs <- soft_seqs[nchar(soft_seqs) >= min_len]
  if (length(soft_seqs) < min_reads) return(NA_character_)

  # Build consensus by taking the most frequent base at each position
  max_len <- max(nchar(soft_seqs))
  # Pad all to max_len with NAs
  mat <- matrix(NA_character_, nrow = length(soft_seqs), ncol = max_len)
  for (i in seq_along(soft_seqs)) {
    s <- soft_seqs[i]
    mat[i, seq_len(nchar(s))] <- strsplit(s, "")[[1]]
  }
  consensus <- vapply(seq_len(max_len), function(pos) {
    col <- mat[, pos]
    col <- col[!is.na(col)]
    if (length(col) == 0) return("N")
    names(sort(table(col), decreasing = TRUE))[1L]
  }, character(1L))
  paste(consensus, collapse = "")
}


# ---------------------------------------------------------------------------
# CIGAR-only candidate extraction (for use_kmers = FALSE)
# ---------------------------------------------------------------------------
.extract_candidates_cigar <- function(reads, genomic_start, ref_len,
                                       min_size = 10L, verbose = FALSE) {
  n_reads <- length(reads)
  if (n_reads == 0L) return(list())
  if (verbose) message("[TALOS] Extracting candidates from CIGAR ops only")

  all_qnames  <- .safe_qnames(reads)
  all_mapqs   <- S4Vectors::mcols(reads)$mapq; all_mapqs[is.na(all_mapqs)] <- 0L
  all_flags   <- S4Vectors::mcols(reads)$flag
  is_reverse  <- bitwAnd(all_flags, 0x10L) != 0L
  read_starts <- BiocGenerics::start(reads)
  read_ends   <- BiocGenerics::end(reads)
  read_cigars <- GenomicAlignments::cigar(reads)
  read_seqs   <- as.character(S4Vectors::mcols(reads)$seq)

  candidates <- vector("list", 2L * n_reads)
  cand_idx   <- 1L

  op_table <- if (.has_cigarillo)
    cigarillo::tabulate_cigar_ops(read_cigars)
  else
    suppressWarnings(GenomicAlignments::cigarOpTable(read_cigars))

  has_ins <- if ("I" %in% colnames(op_table)) op_table[, "I"] > 0L else rep(FALSE, n_reads)
  has_soft <- if ("S" %in% colnames(op_table)) op_table[, "S"] > 0L else rep(FALSE, n_reads)

  for (i in seq_len(n_reads)) {
    # Insertions (I) – breakpoint at the insertion start
    if (has_ins[i]) {
      local_bp <- read_starts[i] - genomic_start + 1L
      if (local_bp >= 1L && local_bp <= ref_len) {
        candidates[[cand_idx]] <- list(
          read_name = all_qnames[i], local_breakpoint = local_bp,
          length = as.integer(op_table[i, "I"]), type = "cigar_ins",
          mapq = all_mapqs[i], is_reverse = is_reverse[i],
          cigar = read_cigars[i], read_seq = read_seqs[i], itd_seq = NA_character_
        )
        cand_idx <- cand_idx + 1L
      }
    }
    # Soft‑clips – same as PTD mode
    if (has_soft[i]) {
      lead_clip <- grepl("^\\d+S", read_cigars[i], perl = TRUE)
      trail_clip <- grepl("\\d+S$", read_cigars[i], perl = TRUE)
      if (lead_clip) {
        local_bp <- read_starts[i] - genomic_start + 1L
        if (local_bp >= 1L && local_bp <= ref_len) {
          candidates[[cand_idx]] <- list(
            read_name = all_qnames[i], local_breakpoint = local_bp, length = 0L,
            type = "ptd_clip", mapq = all_mapqs[i], is_reverse = is_reverse[i],
            cigar = read_cigars[i], read_seq = read_seqs[i], itd_seq = NA_character_
          )
          cand_idx <- cand_idx + 1L
        }
      }
      if (trail_clip) {
        local_bp <- read_ends[i] - genomic_start + 1L
        if (local_bp >= 1L && local_bp <= ref_len) {
          candidates[[cand_idx]] <- list(
            read_name = all_qnames[i], local_breakpoint = local_bp, length = 0L,
            type = "ptd_clip", mapq = all_mapqs[i], is_reverse = is_reverse[i],
            cigar = read_cigars[i], read_seq = read_seqs[i], itd_seq = NA_character_
          )
          cand_idx <- cand_idx + 1L
        }
      }
    }
  }

  candidates <- candidates[seq_len(cand_idx - 1L)]
  if (verbose) message("[TALOS] Found ", length(candidates), " CIGAR‑based candidates")
  candidates
}


# ---------------------------------------------------------------------------
# Anchor‑based ITD size estimation (primary geometry source)
# Uses sliding seed walking on reference, not k‑mer gaps.
# Returns a single size estimate (or NA).
# ---------------------------------------------------------------------------
.estimate_itd_from_anchors <- function(lead_clip, trail_clip, ref_dna,
                                      orig_bp_local, search_window = 5000L,
                                      min_anchor_len = 20L,
                                      kmer_len = NA_integer_,
                                      kmer_band = 40L,
                                      debug_label = NULL,
                                      verbose = FALSE) {
  ref_len <- nchar(ref_dna)
  candidate_sizes <- integer(0)

  .collect_mode <- function(sizes) {
    if (length(sizes) == 0L) return(NA_integer_)
    tb <- sort(table(sizes), decreasing = TRUE)
    as.integer(names(tb)[1L])
  }

  .top5_str <- function(sizes) {
    if (length(sizes) == 0L) return("none")
    tb <- sort(table(sizes), decreasing = TRUE)
    top <- head(tb, 5)
    paste0(names(top), "x", as.integer(top), collapse = ", ")
  }

  # Left anchor: extend lead soft-clip downstream
  if (!is.na(lead_clip) && nchar(lead_clip) >= min_anchor_len) {
    start_pos <- orig_bp_local + 1L
    end_pos <- min(ref_len, start_pos + search_window - 1L)
    if (start_pos <= end_pos) {
      win <- substr(ref_dna, start_pos, end_pos)
      max_prefix <- min(50L, nchar(lead_clip))
      left_sizes <- integer(0)
      for (plen in seq(min_anchor_len, max_prefix, by = 2L)) {
        q <- substr(lead_clip, 1L, plen)
        hits <- gregexpr(q, win, fixed = TRUE)[[1L]]
        if (length(hits) == 0L || hits[1L] == -1L) next
        offsets <- start_pos + hits - 1L - orig_bp_local
        offsets <- offsets[offsets > 0L]
        if (!is.na(kmer_len))
          offsets <- offsets[abs(offsets - as.integer(kmer_len)) <= as.integer(kmer_band)]
        if (length(offsets) > 0L)
          left_sizes <- c(left_sizes, as.integer(offsets))
      }
      candidate_sizes <- c(candidate_sizes, left_sizes)
    }
  }

  # Right anchor: extend trail soft-clip upstream
  if (!is.na(trail_clip) && nchar(trail_clip) >= min_anchor_len) {
    end_pos <- orig_bp_local - 1L
    start_pos <- max(1L, end_pos - search_window + 1L)
    if (start_pos <= end_pos) {
      win <- substr(ref_dna, start_pos, end_pos)
      max_prefix <- min(50L, nchar(trail_clip))
      right_sizes <- integer(0)
      for (plen in seq(min_anchor_len, max_prefix, by = 2L)) {
        q <- substr(trail_clip, 1L, plen)
        hits <- gregexpr(q, win, fixed = TRUE)[[1L]]
        if (length(hits) == 0L || hits[1L] == -1L) next
        offsets <- orig_bp_local - (start_pos + hits - 1L)
        offsets <- offsets[offsets > 0L]
        if (!is.na(kmer_len))
          offsets <- offsets[abs(offsets - as.integer(kmer_len)) <= as.integer(kmer_band)]
        if (length(offsets) > 0L)
          right_sizes <- c(right_sizes, as.integer(offsets))
      }
      candidate_sizes <- c(candidate_sizes, right_sizes)
    }
  }

  candidate_sizes <- candidate_sizes[!is.na(candidate_sizes) & candidate_sizes > 0L]
  if (verbose && !is.null(debug_label)) {
    message("[TALOS] Anchor sizes ", debug_label, ": ", .top5_str(candidate_sizes))
  }
  .collect_mode(candidate_sizes)
}

# ---------------------------------------------------------------------------
# k‑mer trace validation (consistency check, not geometry)
# Returns TRUE if the read's k‑mer pattern is compatible with a duplication.
# ---------------------------------------------------------------------------
.validate_kmer_trace <- function(trace, max_missing_kmers = 0.5,
                                 require_backward_jump = FALSE,
                                 k = 11L) {   # k added as parameter
  if (length(trace) == 0L) return(FALSE)
  missing_rate <- sum(trace == -1L) / length(trace)
  if (missing_rate > max_missing_kmers) return(FALSE)

  if (require_backward_jump) {
    has_jump <- FALSE
    for (j in seq(2L, length(trace))) {
      if (trace[j] == -1L || trace[j-1L] == -1L) next
      if (trace[j-1L] + k - 1L > trace[j]) {
        has_jump <- TRUE
        break
      }
    }
    if (!has_jump) return(FALSE)
  }
  TRUE
}


# ---------------------------------------------------------------------------
# Standard k-mer tracing (backward-jump and missing-block strategies)
# Uses anchor‑based size when soft‑clips are available, falls back to k‑mer.
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

    trace               <- .kmer_match(all_kmers, ref_kmers)
    trace[is.na(trace)] <- -1L

    missing <- sum(trace == -1L)
    if (!ptd_mode && (missing / length(trace) > max_missing_kmers)) next

    read_qname <- all_qnames[i]; read_mapq  <- all_mapqs[i]
    read_cigar <- all_cigars[i]; read_start <- all_starts[i]; read_end <- all_ends[i]

    # --- Extract soft‑clips for anchor‑based estimation ---
    sc <- .get_softclips(read_cigar, read_seq)

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
            # Use anchor‑based size if available (primary), else k‑mer size
            anchor_size <- .estimate_itd_from_anchors(
              sc$lead, sc$trail, ref_dna,
              orig_bp_local = pos_after,
              search_window = 5000L,
              kmer_len = dup_len,
              debug_label = paste0(read_qname, " @", pos_after),
              verbose = verbose
            )
            if (!is.na(anchor_size) && anchor_size >= min_size) {
              final_len <- anchor_size
              dup_seq <- NA_character_
              dup_type <- "anchor_ext"
            } else {
              final_len <- dup_len
              dup_seq <- if (!ptd_mode && dup_len <= 500L)
                substr(read_seq, j, j + dup_len - 1L) else NA_character_
              dup_type <- "bwd_jump"
            }
            candidates[[cand_idx]] <- list(
              read_name = read_qname, local_breakpoint = pos_after,
              itd_seq = dup_seq, read_seq = read_seq,
              length = if (ptd_mode) 0L else as.integer(final_len),
              type = dup_type, mapq = read_mapq,
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
            dup_len <- pos_before - pos_after + (break_end - break_start) + 2L
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
              # Anchor‑based size (primary) or fallback
              anchor_size <- .estimate_itd_from_anchors(
                sc$lead, sc$trail, ref_dna,
                orig_bp_local = pos_after,
                search_window = 5000L,
                kmer_len = dup_len,
                debug_label = paste0(read_qname, " @", pos_after),
                verbose = verbose
              )
              if (!is.na(anchor_size) && anchor_size >= min_size) {
                final_len <- anchor_size
                dup_seq <- NA_character_
                dup_type <- "anchor_ext"
              } else {
                final_len <- dup_len
                dup_seq <- if (!ptd_mode && dup_len <= 500L)
                  substr(read_seq, break_start, break_start + dup_len - 1L) else NA_character_
                dup_type <- "missing_kmer"
              }
              candidates[[cand_idx]] <- list(
                read_name = read_qname, local_breakpoint = pos_after,
                itd_seq = dup_seq, read_seq = read_seq,
                length = if (ptd_mode) 0L else as.integer(final_len),
                type = dup_type, mapq = read_mapq,
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


# ---------------------------------------------------------------------------
# In-silico read extension (uses anchor‑based size and fallback)
# ---------------------------------------------------------------------------
.extend_candidates <- function(candidates_df, ref_dna, genomic_start,
                               min_size = 15L, search_window = 5000L,
                               verbose = FALSE) {
  if (is.null(candidates_df) || nrow(candidates_df) == 0L) return(data.frame())
  results <- vector("list", nrow(candidates_df))
  n_out <- 0L

  for (i in seq_len(nrow(candidates_df))) {
    row <- candidates_df[i, ]
    read_seq <- row$read_seq
    cigar <- row$cigar
    if (is.na(read_seq) || nchar(read_seq) == 0L || is.na(cigar)) next

    sc <- .get_softclips(cigar, read_seq)
    if (is.na(sc$lead) && is.na(sc$trail)) next

    orig_bp <- as.integer(row$breakpoint)
    refined_bp <- orig_bp
    orig_bp_local <- orig_bp - genomic_start + 1L

    # Anchor‑based size (primary)
    anchor_size <- .estimate_itd_from_anchors(
      lead_clip = sc$lead,
      trail_clip = sc$trail,
      ref_dna = ref_dna,
      orig_bp_local = orig_bp_local,
      search_window = search_window,
      kmer_len = row$length,
      debug_label = as.character(row$read_name),
      verbose = verbose
    )

    # Use anchor size if valid, otherwise fall back to k‑mer length from candidates_df
    if (!is.na(anchor_size) && anchor_size >= min_size) {
      itdsize <- anchor_size
    } else if (!is.na(row$length) && row$length >= min_size) {
      itdsize <- as.integer(row$length)
    } else {
      next
    }

    itdsize <- as.integer(itdsize)
    if (itdsize < min_size) next

    n_out <- n_out + 1L
    results[[n_out]] <- list(
      breakpoint_original = orig_bp,
      breakpoint_refined  = refined_bp,
      itdsize             = itdsize,
      read_name           = as.character(row$read_name),
      extended            = NA_character_
    )
  }

  if (n_out == 0L) return(data.frame())
  results <- results[seq_len(n_out)]
  data.frame(
    breakpoint = vapply(results, function(x) {
      v <- x[["breakpoint_original"]]
      if (is.null(v) || (length(v) == 1L && is.na(v))) NA_integer_ else as.integer(v[1L])
    }, integer(1L)),
    breakpoint_refined = vapply(results, function(x) {
      v <- x[["breakpoint_refined"]]
      if (is.null(v) || (length(v) == 1L && is.na(v))) NA_integer_ else as.integer(v[1L])
    }, integer(1L)),
    itdsize = vapply(results, function(x) {
      v <- x[["itdsize"]]
      if (is.null(v) || (length(v) == 1L && is.na(v))) NA_integer_ else as.integer(v[1L])
    }, integer(1L)),
    read_name = vapply(results, function(x) {
      v <- x[["read_name"]]
      if (is.null(v) || length(v) == 0L || is.na(v[1L])) NA_character_ else as.character(v[1L])
    }, character(1L)),
    extended = vapply(results, function(x) {
      v <- x[["extended"]]
      if (is.null(v) || length(v) == 0L || is.na(v[1L])) NA_character_ else as.character(v[1L])
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}
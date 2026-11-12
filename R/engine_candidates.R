.safe_qnames <- function(reads) {
  qn <- as.character(S4Vectors::mcols(reads)$qname)
  if (length(qn) == 0L || all(is.na(qn))) {
    qn <- names(reads)
  }
  na_idx <- which(is.na(qn))
  if (length(na_idx) > 0L) {
    qn[na_idx] <- paste0("read_", na_idx)
  }
  qn
}


.unique_fixed_match <- function(query, subject) {
  hits <- integer(0)
  search_from <- 1L
  subject_len <- nchar(subject)

  while (search_from <= subject_len && length(hits) < 2L) {
    hit <- regexpr(query, substring(subject, search_from), fixed = TRUE)[1L]
    if (hit < 0L) {
      break
    }
    absolute_hit <- search_from + hit - 1L
    hits <- c(hits, absolute_hit)
    search_from <- absolute_hit + 1L
  }

  hits
}


.apply_itd_length_limit <- function(
  length,
  type,
  max_itd_length,
  convert_long_to_ptd
) {
  if (length > max_itd_length) {
    if (!convert_long_to_ptd) {
      return(NULL)
    }
    return(list(length = 0L, type = "ptd_clip", converted = TRUE))
  }

  list(length = as.integer(length), type = type, converted = FALSE)
}


# One constructor for every candidate record (DRY: the record shape was
# spelled out by hand in nine places). read_start / read_end are the read's
# aligned genomic span; the anchor search needs them (see
# .estimate_itd_from_anchors()). read_mean_bq carries the read's mean Phred
# base quality (v0.4.20) so base-quality metrics report real measurements
# instead of a hard-coded "unavailable" status.
.new_candidate <- function(read_name, local_breakpoint, length, type, mapq,
                           is_reverse, cigar, read_seq, itd_seq = NA_character_,
                           read_start = NA_integer_, read_end = NA_integer_,
                           read_mean_bq = NA_real_) {
  list(read_name = read_name, local_breakpoint = as.integer(local_breakpoint),
       length = length, type = type, mapq = mapq, is_reverse = is_reverse,
       cigar = cigar, read_seq = read_seq, itd_seq = itd_seq,
       read_start = as.integer(read_start), read_end = as.integer(read_end),
       read_mean_bq = as.numeric(read_mean_bq))
}

# Per-read vectors used by every candidate extractor (DRY).
.read_fields <- function(reads) {
  mapqs <- S4Vectors::mcols(reads)$mapq
  mapqs[is.na(mapqs)] <- 0L
  cigars <- GenomicAlignments::cigar(reads)
  # Mean Phred base quality per read; NA when the BAM QUAL column was not
  # loaded (e.g. precomputed candidate frames) so downstream status fields
  # can distinguish "unavailable" from "computed".
  mean_bqs <- .mean_base_quality(S4Vectors::mcols(reads)$qual)
  list(
    qnames = .safe_qnames(reads),
    mapqs = mapqs,
    flags = S4Vectors::mcols(reads)$flag,
    is_reverse = bitwAnd(S4Vectors::mcols(reads)$flag, 0x10L) != 0L,
    starts = BiocGenerics::start(reads),
    ends = BiocGenerics::end(reads),
    cigars = cigars,
    seqs = as.character(S4Vectors::mcols(reads)$seq),
    mean_bqs = mean_bqs,
    clips = lapply(cigars, .get_softclip_lengths)
  )
}

# Zero-length clip candidates at the aligned start (leading clip) and/or end
# (trailing clip) of read i. Shared by the PTD and CIGAR extractors.
.clip_candidates <- function(rf, i, genomic_start, ref_len) {
  out <- list()
  sides <- c(lead = rf$starts[i], trail = rf$ends[i])
  for (side in names(sides)) {
    if (rf$clips[[i]][[side]] <= 0L) next
    local_bp <- sides[[side]] - genomic_start + 1L
    if (local_bp < 1L || local_bp > ref_len) next
    out[[length(out) + 1L]] <- .new_candidate(
      rf$qnames[i], local_bp, 0L, "ptd_clip", rf$mapqs[i], rf$is_reverse[i],
      rf$cigars[i], rf$seqs[i], read_start = rf$starts[i], read_end = rf$ends[i],
      read_mean_bq = rf$mean_bqs[i]
    )
  }
  out
}

.extract_candidates_ptd <- function(
  reads,
  genomic_start,
  ref_len,
  verbose = FALSE
) {
  if (length(reads) == 0L) return(list())
  if (verbose) message("[TALOS] Extracting PTD candidates from ", length(reads), " reads")
  rf <- .read_fields(reads)
  candidates <- unlist(lapply(seq_along(reads), function(i) {
    .clip_candidates(rf, i, genomic_start, ref_len)
  }), recursive = FALSE)
  if (is.null(candidates)) candidates <- list()
  if (verbose) message("[TALOS] Found ", length(candidates), " PTD candidates")
  candidates
}

.assemble_ptd_consensus <- function(
  support_rows,
  min_reads = 3L,
  min_len = 15L
) {
  soft_seqs <- character()
  for (i in seq_len(nrow(support_rows))) {
    sc <- .get_softclips(support_rows$cigar[i], support_rows$read_seq[i])
    if (!is.na(sc$lead)) {
      soft_seqs <- c(soft_seqs, sc$lead)
    }
    if (!is.na(sc$trail)) soft_seqs <- c(soft_seqs, sc$trail)
  }
  if (length(soft_seqs) < min_reads) {
    return(NA_character_)
  }

  soft_seqs <- soft_seqs[nchar(soft_seqs) >= min_len]
  if (length(soft_seqs) < min_reads) {
    return(NA_character_)
  }

  .assemble_weighted_consensus(soft_seqs)
}


.extract_candidates_cigar <- function(
  reads,
  genomic_start,
  ref_len,
  min_size = 10L,
  verbose = FALSE
) {
  n_reads <- length(reads)
  if (n_reads == 0L) return(list())
  if (verbose) message("[TALOS] Extracting candidates from CIGAR ops only")
  rf <- .read_fields(reads)
  op_table <- if (.has_cigarillo) {
    cigarillo::tabulate_cigar_ops(rf$cigars)
  } else {
    suppressWarnings(GenomicAlignments::cigarOpTable(rf$cigars))
  }
  has_ins <- if ("I" %in% colnames(op_table)) op_table[, "I"] > 0L else rep(FALSE, n_reads)

  candidates <- list()
  for (i in seq_len(n_reads)) {
    if (has_ins[i]) {
      insertions <- .parse_cigar_insertions(rf$cigars[i], rf$starts[i], genomic_start,
                                            ref_len, min_size = min_size)
      for (j in seq_len(nrow(insertions))) {
        candidates[[length(candidates) + 1L]] <- .new_candidate(
          rf$qnames[i], insertions$local_breakpoint[j], insertions$length[j], "cigar_ins",
          rf$mapqs[i], rf$is_reverse[i], rf$cigars[i], rf$seqs[i],
          read_start = rf$starts[i], read_end = rf$ends[i],
          read_mean_bq = rf$mean_bqs[i]
        )
      }
    }
    candidates <- c(candidates, .clip_candidates(rf, i, genomic_start, ref_len))
  }
  if (verbose) message("[TALOS] Found ", length(candidates), " CIGAR-based candidates")
  candidates
}

#' Estimate a duplication length from soft-clip re-anchoring
#'
#' For a tandem duplication [S, E] (length L = E - S + 1):
#' \itemize{
#'   \item a read aligned from S with a leading clip C (length c) carries the
#'     end of the first copy, C = ref[E - c + 1 .. E]. The prefix of C is found
#'     downstream at h = E - c + 1, hence L = h - S + c;
#'   \item a read aligned up to E with a trailing clip C carries the start of
#'     the second copy, C = ref[S .. S + c - 1]. The prefix of C is found
#'     upstream of E at h = S, hence L = E - h + 1.
#' }
#' Bug fixed: the lead side returned h - S (short by the clip length; the
#' 40 bp k-mer band is why HORIZON FLT3 300 bp came out at ~260 bp) and the
#' trail side searched upstream of the breakpoint S instead of the read end E
#' (it could only hit repeats).
#'
#' [BUG 8] Two-pass search: pass 1 keeps hits within kmer_len +/- kmer_band;
#' when a side has none, pass 2 retries without the band.
#'
#' @param lead_clip,trail_clip Leading / trailing soft-clip sequences (or NA).
#' @param ref_dna Local reference sequence.
#' @param orig_bp_local Local breakpoint (default anchor for both sides).
#' @param lead_bp_local Local aligned start of the read (S side).
#' @param trail_bp_local Local aligned end of the read (E side); the trail
#'   side is skipped when NA.
#' @return Mode of the candidate sizes, or NA.
.estimate_itd_from_anchors <- function(
  lead_clip,
  trail_clip,
  ref_dna,
  orig_bp_local,
  search_window = 5000L,
  min_anchor_len = 20L,
  kmer_len = NA_integer_,
  kmer_band = 40L,
  require_unique_anchor = FALSE,
  debug_label = NULL,
  verbose = FALSE,
  lead_bp_local = orig_bp_local,
  trail_bp_local = NA_integer_
) {
  ref_len <- nchar(ref_dna)
  .top5_str <- function(sizes) {
    if (length(sizes) == 0L) return("none")
    top <- utils::head(sort(table(sizes), decreasing = TRUE), 5)
    paste0(names(top), "x", as.integer(top), collapse = ", ")
  }
  .collect_side <- function(clip_seq, side, use_band) {
    bp <- if (side == "lead") lead_bp_local else trail_bp_local
    if (is.na(clip_seq) || is.na(bp) || nchar(clip_seq) < min_anchor_len) return(integer(0))
    if (side == "lead") {
      start_pos <- bp + 1L
      end_pos <- min(ref_len, start_pos + search_window - 1L)
    } else {
      end_pos <- bp - 1L
      start_pos <- max(1L, end_pos - search_window + 1L)
    }
    if (start_pos > end_pos) return(integer(0))
    win <- substr(ref_dna, start_pos, end_pos)
    clip_len <- nchar(clip_seq)
    sizes <- integer(0)
    for (plen in seq(min_anchor_len, min(50L, clip_len), by = 2L)) {
      q <- substr(clip_seq, 1L, plen)
      if (require_unique_anchor) {
        hits <- .unique_fixed_match(q, win)
        if (length(hits) != 1L) next
      } else {
        hits <- gregexpr(q, win, fixed = TRUE)[[1L]]
        if (length(hits) == 0L || hits[1L] == -1L) next
      }
      abs_hits <- start_pos + hits - 1L
      offsets <- if (side == "lead") abs_hits - bp + clip_len else bp - abs_hits + 1L
      offsets <- offsets[offsets > 0L]
      if (use_band && !is.na(kmer_len)) {
        offsets <- offsets[abs(offsets - as.integer(kmer_len)) <= as.integer(kmer_band)]
      }
      sizes <- c(sizes, as.integer(offsets))
    }
    sizes
  }
  candidate_sizes <- integer(0)
  recovered_unbanded <- FALSE
  for (side in c("lead", "trail")) {
    clip <- if (side == "lead") lead_clip else trail_clip
    sz <- .collect_side(clip, side, use_band = TRUE)
    if (length(sz) == 0L && !is.na(kmer_len)) {
      sz <- .collect_side(clip, side, use_band = FALSE)
      if (length(sz) > 0L) recovered_unbanded <- TRUE
    }
    candidate_sizes <- c(candidate_sizes, sz)
  }
  candidate_sizes <- candidate_sizes[!is.na(candidate_sizes) & candidate_sizes > 0L]
  if (verbose && !is.null(debug_label)) {
    if (recovered_unbanded) {
      message("[TALOS] Anchor ", debug_label, ": kmer band (", kmer_len, " +/- ", kmer_band,
              ") excluded all hits; recovered un-banded")
    }
    message("[TALOS] Anchor sizes ", debug_label, ": ", .top5_str(candidate_sizes))
  }
  if (length(candidate_sizes) == 0L) return(NA_integer_)
  as.integer(names(sort(table(candidate_sizes), decreasing = TRUE))[1L])
}

# ---------------------------------------------------------------------------
# Standard k-mer tracing (backward-jump and missing-block strategies)
# Uses anchor-based size when soft-clips are available, falls back to k-mer.
.extract_candidates_standard <- function(
  reads,
  ref_kmers,
  ptd_mode,
  min_size,
  max_missing_kmers,
  refine_bp,
  use_cigar_bp,
  genomic_start,
  ref_dna,
  max_itd_length = 1000L,
  convert_long_to_ptd = TRUE,
  search_window = 5000L,
  require_unique_anchor = FALSE,
  verbose = FALSE
) {
  n_reads <- length(reads)
  k <- nchar(ref_kmers[1L])
  if (verbose) {
    message("[TALOS] Extracting standard candidates from ", n_reads, " reads (k=", k,
            ", ptd_mode=", ptd_mode, ", max_itd_length=", max_itd_length, ")")
  }
  rf <- .read_fields(reads)
  candidates <- vector("list", 2L * n_reads)
  cand_idx <- 1L
  ref_len <- nchar(ref_dna)
  reads_processed <- 0L
  add <- function(cand) {
    candidates[[cand_idx]] <<- cand
    cand_idx <<- cand_idx + 1L
  }

  for (i in seq_len(n_reads)) {
    read_seq <- rf$seqs[i]
    if (is.na(read_seq)) next
    read_len <- nchar(read_seq)
    if (read_len < k || read_len < (k + min_size)) next

    all_kmers <- if (.has_biostrings) {
      as.character(Biostrings::Views(Biostrings::DNAString(read_seq),
                                     start = seq_len(read_len - k + 1L), width = k))
    } else {
      kmer_starts <- seq_len(read_len - k + 1L)
      substring(read_seq, kmer_starts, kmer_starts + k - 1L)
    }
    trace <- .kmer_match(all_kmers, ref_kmers)
    trace[is.na(trace)] <- -1L
    if (!ptd_mode && (sum(trace == -1L) / length(trace) > max_missing_kmers)) next

    read_start <- rf$starts[i]
    read_end <- rf$ends[i]
    sc <- .get_softclips(rf$cigars[i], read_seq)
    mk <- function(local_bp, len, type, itd_seq = NA_character_) {
      .new_candidate(rf$qnames[i], local_bp, len, type, rf$mapqs[i], rf$is_reverse[i],
                     rf$cigars[i], read_seq, itd_seq, read_start, read_end,
                     read_mean_bq = rf$mean_bqs[i])
    }

    # Shared by both k-mer strategies: anchor refinement, ITD length limit,
    # candidate creation. Returns FALSE when the length limit drops it.
    emit_dup <- function(pos_after, dup_len, seq_from, kmer_type) {
      anchor_size <- .estimate_itd_from_anchors(
        sc$lead, sc$trail, ref_dna,
        orig_bp_local = pos_after,
        search_window = search_window,
        require_unique_anchor = require_unique_anchor,
        kmer_len = dup_len,
        debug_label = paste0(rf$qnames[i], " @", pos_after),
        verbose = verbose,
        lead_bp_local = read_start - genomic_start + 1L,
        trail_bp_local = read_end - genomic_start + 1L
      )
      if (!is.na(anchor_size) && anchor_size >= min_size) {
        final_len <- anchor_size
        dup_seq <- NA_character_
        dup_type <- "anchor_ext"
      } else {
        final_len <- dup_len
        dup_seq <- if (!ptd_mode && dup_len <= 500L) substr(read_seq, seq_from, seq_from + dup_len - 1L) else NA_character_
        dup_type <- kmer_type
      }
      if (!ptd_mode) {
        length_result <- .apply_itd_length_limit(final_len, dup_type, max_itd_length, convert_long_to_ptd)
        if (is.null(length_result)) return(FALSE)
        if (length_result$converted) {
          dup_seq <- NA_character_
          if (verbose) {
            message("[TALOS] Anchor-refined long dup (", final_len, " bp) -> PTD at ",
                    genomic_start + pos_after - 1L)
          }
        }
        final_len <- length_result$length
        dup_type <- length_result$type
      }
      add(mk(pos_after, if (ptd_mode) 0L else as.integer(final_len), dup_type, dup_seq))
      TRUE
    }

    # Strategy 1: backward jump between consecutive matched k-mers.
    if (length(trace) >= 2L) {
      for (j in seq(2L, length(trace))) {
        if (trace[j] == -1L || trace[j - 1L] == -1L) next
        if (trace[j - 1L] + k - 1L > trace[j]) {
          dup_len <- (trace[j - 1L] + k - 1L) - trace[j] + 1L
          if (dup_len >= min_size) emit_dup(trace[j], dup_len, j, "bwd_jump")
        }
      }
    }

    # Strategy 2: block of missing k-mers bracketed by a backward jump.
    in_break <- FALSE
    break_start <- NULL
    for (j in seq_along(trace)) {
      if (trace[j] == -1L && !in_break) {
        in_break <- TRUE
        break_start <- j
      } else if (trace[j] != -1L && in_break) {
        break_end <- j - 1L
        if (break_start > 1L && break_end < length(trace)) {
          pos_before <- trace[break_start - 1L]
          pos_after <- trace[break_end + 1L]
          if (pos_before != -1L && pos_after != -1L && (pos_before + k - 1L > pos_after)) {
            dup_len <- pos_before - pos_after + (break_end - break_start) + 2L
            if (dup_len >= min_size) emit_dup(pos_after, dup_len, break_start, "missing_kmer")
          }
        } else if (ptd_mode && break_start == 1L) {
          pos_after <- trace[j]
          refined_bp <- pos_after
          if (use_cigar_bp && !is.na(sc$lead)) {
            refined_bp <- read_start - genomic_start + 1L
          } else if (refine_bp) {
            max_offset <- min(40L, max(0L, pos_after - 1L))
            for (offset in 0:max_offset) {
              test_pos <- pos_after - offset
              if (test_pos + k - 1L <= ref_len &&
                  substr(ref_dna, test_pos, test_pos + k - 1L) == substr(read_seq, j, j + k - 1L)) {
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
          add(mk(max(1L, refined_bp), 0L, "ptd_clip"))
        }
        in_break <- FALSE
      }
    }

    if (in_break && ptd_mode && break_start > 1L) {
      pos_before <- trace[break_start - 1L]
      if (pos_before != -1L) {
        refined_bp <- pos_before + k - 1L
        if (use_cigar_bp && !is.na(sc$trail)) refined_bp <- read_end - genomic_start + 1L
        add(mk(max(1L, refined_bp), 0L, "ptd_clip"))
      }
    }

    reads_processed <- reads_processed + 1L
    if (verbose && reads_processed %% 1000L == 0L) {
      message("[TALOS] Processed ", reads_processed, " / ", n_reads, " reads for candidate extraction")
    }
  }

  candidates <- candidates[seq_len(cand_idx - 1L)]
  if (verbose) message("[TALOS] Found ", length(candidates), " standard candidates")
  candidates
}

# Typed field extraction from a list of records, NULL/NA safe (DRY helper for
# the data.frame builders below).
.pluck <- function(records, field, type = c("integer", "character", "logical", "numeric")) {
  type <- match.arg(type)
  na <- switch(type, integer = NA_integer_, character = NA_character_, logical = NA, numeric = NA_real_)
  cast <- switch(type, integer = as.integer, character = as.character, logical = as.logical, numeric = as.numeric)
  vapply(records, function(x) {
    v <- x[[field]]
    if (is.null(v) || length(v) == 0L || is.na(v[1L])) na else cast(v[1L])
  }, switch(type, integer = integer(1L), character = character(1L), logical = logical(1L), numeric = numeric(1L)),
  USE.NAMES = FALSE)
}

.candidates_to_df <- function(candidates, genomic_start) {
  if (length(candidates) == 0L) return(data.frame())
  data.frame(
    breakpoint = as.integer(as.numeric(genomic_start) + .pluck(candidates, "local_breakpoint") - 1L),
    length = .pluck(candidates, "length"),
    read_name = .pluck(candidates, "read_name", "character"),
    read_seq = .pluck(candidates, "read_seq", "character"),
    type = .pluck(candidates, "type", "character"),
    mapq = .pluck(candidates, "mapq"),
    is_reverse = .pluck(candidates, "is_reverse", "logical"),
    cigar = .pluck(candidates, "cigar", "character"),
    itd_seq = .pluck(candidates, "itd_seq", "character"),
    read_start = .pluck(candidates, "read_start"),
    read_end = .pluck(candidates, "read_end"),
    dup_exons = .pluck(candidates, "dup_exons", "character"),
    read_mean_bq = .pluck(candidates, "read_mean_bq", "numeric"),
    stringsAsFactors = FALSE
  )
}

# Breakpoint clustering.
#
# v0.4.20 adds length-scaled tolerance: breakpoint placement uncertainty
# grows with the duplication length, so the effective tolerance for a
# breakpoint with duplication length L is
#   max(cluster_tolerance, length_scale * L)
# capped at max_scaled_tolerance. The anti-chaining rule below still bounds
# the TOTAL cluster width: every member must lie within the (scaled)
# tolerance of the distance from the cluster seed, so nearby distinct
# events can no longer chain into one arbitrarily wide cluster.
#
# length_scale = 0 (or lengths = NULL) restores the legacy scalar tolerance.
.cluster_breakpoints <- function(bp_values, cluster_tolerance, lengths = NULL,
                                 length_scale = 0, max_scaled_tolerance = Inf) {
  tolerance <- suppressWarnings(as.numeric(cluster_tolerance))
  if (length(tolerance) != 1L || is.na(tolerance) || !is.finite(tolerance) || tolerance < 0) {
    stop("cluster_tolerance must be a single finite non-negative number.", call. = FALSE)
  }
  length_scale <- suppressWarnings(as.numeric(length_scale))
  if (length(length_scale) != 1L || is.na(length_scale) || length_scale < 0 || length_scale > 1) {
    stop("length_scale must be a single number between 0 and 1.", call. = FALSE)
  }
  max_scaled_tolerance <- suppressWarnings(as.numeric(max_scaled_tolerance))
  if (length(max_scaled_tolerance) != 1L || is.na(max_scaled_tolerance) ||
      (is.finite(max_scaled_tolerance) && max_scaled_tolerance < tolerance)) {
    stop("max_scaled_tolerance must be a single number (>= cluster_tolerance, or Inf for no cap).", call. = FALSE)
  }
  bp_values <- suppressWarnings(as.numeric(bp_values))
  keep <- is.finite(bp_values)
  bp_values <- bp_values[keep]
  if (!is.null(lengths) && length(lengths) == length(keep)) {
    lengths <- suppressWarnings(as.numeric(lengths[keep]))
  } else {
    lengths <- NULL
  }
  if (length(bp_values) == 0L) return(list())
  ord <- order(bp_values)
  bp_sorted <- bp_values[ord]
  lengths_sorted <- if (!is.null(lengths) && length_scale > 0) {
    ls <- lengths[ord]
    ls[is.na(ls)] <- 0 # NA length (unknown size) keeps the base tolerance
    ls
  } else NULL
  if (length(bp_sorted) == 1L) return(list(bp_sorted))

  # Per-breakpoint effective tolerance (length-scaled, capped).
  tol_sorted <- if (!is.null(lengths_sorted)) {
    pmin(max_scaled_tolerance, pmax(tolerance, length_scale * lengths_sorted))
  } else {
    rep(tolerance, length(bp_sorted))
  }

  # Bound the TOTAL cluster width, not merely each adjacent gap. The old
  # single-link rule chained a run of nearby breakpoints into one arbitrarily
  # wide cluster (e.g. 100,108,116 with tolerance=10). That can pool separate
  # nearby events and inflate their support.
  clusters <- list()
  cluster_start_idx <- 1L
  for (j in seq.int(2L, length(bp_sorted))) {
    # Width bound uses the wider of the seed's and the candidate member's
    # effective tolerance: a long duplication anywhere in the cluster widens
    # its own tolerance without un-bounding the cluster.
    tol_pair <- max(tol_sorted[cluster_start_idx], tol_sorted[j])
    if (bp_sorted[j] - bp_sorted[cluster_start_idx] > tol_pair) {
      clusters[[length(clusters) + 1L]] <- bp_sorted[cluster_start_idx:(j - 1L)]
      cluster_start_idx <- j
    }
  }
  clusters[[length(clusters) + 1L]] <- bp_sorted[cluster_start_idx:length(bp_sorted)]
  clusters
}

# Select a length mode without allowing input row order to resolve ties.
# Anchor-backed estimates take precedence over k-mer fallbacks; if several
# modal sizes remain at the same evidence tier, choose the smaller observed
# mode and retain the disagreement for the PE concordance diagnostic.
.mode_extension_length <- function(x, source = NULL) {
  x <- suppressWarnings(as.numeric(x))
  keep <- is.finite(x) & x > 0
  if (!any(keep)) return(NA_integer_)
  x <- as.integer(x[keep])
  tab <- table(x)
  modal_sizes <- as.integer(names(tab)[tab == max(tab)])
  if (!is.null(source) && length(source) == length(keep)) {
    source <- as.character(source[keep])
    anchored <- unique(x[source == "anchor" & !is.na(source)])
    anchored <- anchored[anchored %in% modal_sizes]
    if (length(anchored)) return(as.integer(min(anchored)))
  }
  as.integer(min(modal_sizes))
}

# Reduce multiple candidate rows from one QNAME to one extension hypothesis.
# A paired fragment or a read with several CIGAR/k-mer candidates must not
# cast multiple votes for the same event length.
.collapse_extension_candidates_by_read <- function(ext_df) {
  empty <- data.frame(
    read_name = character(), breakpoint_refined = numeric(),
    itdsize = integer(), itdsize_source = character(),
    stringsAsFactors = FALSE
  )
  if (is.null(ext_df) || nrow(ext_df) == 0L ||
      !all(c("read_name", "breakpoint_refined", "itdsize") %in% names(ext_df))) {
    return(empty)
  }
  keep <- !is.na(ext_df$read_name) & nzchar(as.character(ext_df$read_name)) &
    !is.na(ext_df$breakpoint_refined) & is.finite(ext_df$breakpoint_refined) &
    !is.na(ext_df$itdsize) & is.finite(ext_df$itdsize) & ext_df$itdsize > 0L
  ext_df <- ext_df[keep, , drop = FALSE]
  if (nrow(ext_df) == 0L) return(empty)

  groups <- split(seq_len(nrow(ext_df)), as.character(ext_df$read_name))
  rows <- lapply(names(groups), function(read_name) {
    idx <- groups[[read_name]]
    sizes <- as.integer(ext_df$itdsize[idx])
    tab <- table(sizes)
    modal_sizes <- as.integer(names(tab)[tab == max(tab)])
    # Prefer a modal size backed by an explicit sequence anchor when there
    # is a tie; otherwise choose the smaller observed size. Never invent an
    # unobserved midpoint between competing length hypotheses.
    anchor_sizes <- integer(0)
    if ("itdsize_source" %in% names(ext_df)) {
      is_anchor <- !is.na(ext_df$itdsize_source[idx]) &
        ext_df$itdsize_source[idx] == "anchor"
      anchor_sizes <- unique(sizes[sizes %in% modal_sizes & is_anchor])
      anchor_sizes <- anchor_sizes[!is.na(anchor_sizes)]
    }
    selected_size <- if (length(anchor_sizes)) min(anchor_sizes) else min(modal_sizes)
    selected <- idx[sizes == selected_size]
    source <- if ("itdsize_source" %in% names(ext_df) &&
                  any(ext_df$itdsize_source[selected] == "anchor", na.rm = TRUE)) {
      "anchor"
    } else if ("itdsize_source" %in% names(ext_df)) {
      as.character(ext_df$itdsize_source[selected[1L]])
    } else {
      NA_character_
    }
    data.frame(
      read_name = read_name,
      breakpoint_refined = stats::median(ext_df$breakpoint_refined[selected], na.rm = TRUE),
      itdsize = selected_size,
      itdsize_source = source,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


.prepare_wildtype_info <- function(all_reads, genomic_start, genomic_end) {
  raw_qnames <- .safe_qnames(all_reads)
  all_flags <- S4Vectors::mcols(all_reads)$flag
  is_primary <- bitwAnd(all_flags, 0x100L) == 0L &
    bitwAnd(all_flags, 0x800L) == 0L
  list(
    qnames = raw_qnames,
    is_primary = is_primary,
    gr = GenomicRanges::granges(all_reads),
    region_start = genomic_start,
    region_end = genomic_end
  )
}


.count_wildtype_at_position <- function(
  wildtype_info,
  position,
  support_qnames = character()
) {
  if (
    is.null(wildtype_info$gr) ||
      length(wildtype_info$gr) == 0L ||
      is.na(position)
  ) {
    return(0L)
  }

  keep <- wildtype_info$is_primary &
    GenomicRanges::start(wildtype_info$gr) <= position &
    GenomicRanges::end(wildtype_info$gr) >= position
  keep[is.na(keep)] <- FALSE
  qnames <- wildtype_info$qnames[keep]
  qnames <- unique(qnames[!is.na(qnames) & !qnames %in% support_qnames])
  as.integer(length(qnames))
}


.extend_candidates <- function(
  candidates_df,
  ref_dna,
  genomic_start,
  min_size = 15L,
  search_window = 5000L,
  require_unique_anchor = FALSE,
  verbose = FALSE
) {
  if (is.null(candidates_df) || nrow(candidates_df) == 0L) return(data.frame())
  has_span <- all(c("read_start", "read_end") %in% names(candidates_df))
  results <- list()
  for (i in seq_len(nrow(candidates_df))) {
    row <- candidates_df[i, ]
    if (is.na(row$read_seq) || nchar(row$read_seq) == 0L || is.na(row$cigar)) next
    sc <- .get_softclips(row$cigar, row$read_seq)
    if (is.na(sc$lead) && is.na(sc$trail)) next

    orig_bp <- as.integer(row$breakpoint)
    orig_bp_local <- orig_bp - genomic_start + 1L
    # Lead anchor at the read's aligned start, trail anchor at its aligned
    # end (falls back to breakpoint + length - 1 for records without a span).
    lead_bp <- if (has_span && !is.na(row$read_start)) row$read_start - genomic_start + 1L else orig_bp_local
    trail_bp <- if (has_span && !is.na(row$read_end)) {
      row$read_end - genomic_start + 1L
    } else if (!is.na(row$length) && row$length > 0L) {
      orig_bp_local + row$length - 1L
    } else {
      NA_integer_
    }
    anchor_size <- .estimate_itd_from_anchors(
      lead_clip = sc$lead, trail_clip = sc$trail, ref_dna = ref_dna,
      orig_bp_local = orig_bp_local, search_window = search_window,
      require_unique_anchor = require_unique_anchor, kmer_len = row$length,
      debug_label = as.character(row$read_name), verbose = verbose,
      lead_bp_local = lead_bp, trail_bp_local = trail_bp
    )
    if (!is.na(anchor_size) && anchor_size >= min_size) {
      itdsize <- anchor_size
      itdsize_source <- "anchor"
    } else if (!is.na(row$length) && row$length >= min_size) {
      itdsize <- as.integer(row$length)
      itdsize_source <- "kmer_fallback"
    } else {
      next
    }
    results[[length(results) + 1L]] <- list(
      breakpoint = orig_bp, breakpoint_refined = orig_bp, itdsize = as.integer(itdsize),
      itdsize_source = itdsize_source, read_name = as.character(row$read_name),
      extended = NA_character_
    )
  }
  if (length(results) == 0L) return(data.frame())
  data.frame(
    breakpoint = .pluck(results, "breakpoint"),
    breakpoint_refined = .pluck(results, "breakpoint_refined"),
    itdsize = .pluck(results, "itdsize"),
    itdsize_source = .pluck(results, "itdsize_source", "character"),
    read_name = .pluck(results, "read_name", "character"),
    extended = .pluck(results, "extended", "character"),
    stringsAsFactors = FALSE
  )
}

# Filter empty soft-clips while keeping per-read weights in lockstep.
.filter_weighted_sequences <- function(seqs, weights) {
  if (length(seqs) != length(weights)) {
    stop("seqs and weights must have the same length.", call. = FALSE)
  }
  keep <- !is.na(seqs) & nzchar(seqs)
  list(seqs = seqs[keep], weights = weights[keep])
}

#' Turn exon-graph PTD cycles into read-level candidates
#'
#' One candidate per supporting read present in \code{reads}, placed at the
#' duplicated interval found by the graph (dup_start, length dup_length).
#' Previously the breakpoint was a KDE of the reads' own start/end positions,
#' the length was 0 and every read's sequence was replaced by the junction
#' consensus, so soft-clip metrics read the wrong sequence.
#'
#' @param cycles Output of extract_graph_candidates().
#' @param reads GAlignments to look the supporting reads up in (use the
#'   MAPQ-filtered set, not the CIGAR-prefiltered one: discordant mates have
#'   no clip and would be lost).
#' @param genomic_start Window start.
#' @return list(candidates, cycles_used)
.graph_cycles_to_candidates <- function(cycles, reads, genomic_start, ref_len = NA_integer_) {
  if (!length(cycles) || length(reads) == 0L) return(list(candidates = list(), cycles_used = 0L))
  rf <- .read_fields(reads)
  out <- list()
  used <- 0L
  for (cand in cycles) {
    qn <- unique(unlist(strsplit(cand$support_qnames, ";", fixed = TRUE)))
    idx <- which(rf$qnames %in% qn)
    idx <- idx[!duplicated(rf$qnames[idx])]
    if (!length(idx) || is.na(cand$dup_start) || is.na(cand$dup_length)) next
    sc_seqs <- vapply(idx, function(i) {
      sc <- .get_softclips(rf$cigars[i], rf$seqs[i])
      paste(stats::na.omit(c(sc$lead, sc$trail)), collapse = "")
    }, character(1L))
    # Keep the same mask for the sequence and MAPQ vectors. Re-applying
    # nzchar() after filtering returned an all-TRUE mask over the original
    # read indices, causing a length mismatch and silently losing weights.
    filtered_sc <- .filter_weighted_sequences(sc_seqs, rf$mapqs[idx])
    sc_seqs <- filtered_sc$seqs
    sc_weights <- filtered_sc$weights
    junction_seq <- if (length(sc_seqs)) .build_local_debruijn(sc_seqs) else NA_character_
    if (is.na(junction_seq) && length(sc_seqs)) {
      junction_seq <- .assemble_weighted_consensus(sc_seqs, sc_weights)
    }
    local_bp <- cand$dup_start - genomic_start + 1L
    dup_len <- as.integer(cand$dup_length)
    if (!is.na(ref_len) && (local_bp < 1L || local_bp > ref_len)) next
    if (!is.na(ref_len)) dup_len <- min(dup_len, ref_len - local_bp + 1L)
    if (is.na(dup_len) || dup_len < 1L) next
    n_split <- if (!is.null(cand$n_split)) cand$n_split else NA_integer_
    cand_type <- if (!is.na(n_split) && n_split > 0L) "ptd_graph" else "ptd_graph_imprecise"
    for (i in idx) {
      c1 <- .new_candidate(rf$qnames[i], local_bp, dup_len, cand_type,
                           rf$mapqs[i], rf$is_reverse[i], rf$cigars[i], rf$seqs[i],
                           itd_seq = junction_seq, read_start = rf$starts[i], read_end = rf$ends[i])
      c1$dup_exons <- paste0(sub("Exon_", "", cand$dup_start_exon), "-", sub("Exon_", "", cand$dup_end_exon))
      out[[length(out) + 1L]] <- c1
    }
    used <- used + 1L
  }
  list(candidates = out, cycles_used = used)
}

#' Kernel Density Estimation for Breakpoint Refinement
#' @param positions Numeric vector of breakpoint positions
#' @return Refined breakpoint (mode of density)
.kde_breakpoint <- function(positions) {
  positions <- positions[!is.na(positions)]
  if (length(positions) == 0) return(NA_integer_)
  # density() needs a non-zero bandwidth: identical positions are the answer.
  if (length(unique(positions)) == 1L) return(as.integer(positions[1]))
  d <- stats::density(positions, bw = "nrd")
  as.integer(d$x[which.max(d$y)])
}

#' Build Local de Bruijn Graph from Soft-Clip Sequences
#'
#' @param seqs Character vector of soft-clip sequences
#' @param k K-mer size (default 31)
#' @param min_coverage Minimum k-mer coverage to keep
#' @param max_steps Hard cap on walk length (safety bound against cycles)
#' @return Assembled sequence, or NA if assembly isn't possible
.build_local_debruijn <- function(seqs, k = 31, min_coverage = 2, max_steps = 2000L) {
  if (length(seqs) == 0) return(NA_character_)
  all_seqs <- seqs[!is.na(seqs)]
  # Fallbacks get every input sequence: the k-filtered set was passed before,
  # which is empty on this branch and always returned NA.
  fallback <- function() .simple_overlap_assembly(all_seqs, min_overlap = min(k, 10L))
  seqs <- all_seqs[nchar(all_seqs) >= k]
  if (length(seqs) == 0) return(fallback())
  kmer_env <- new.env(hash = TRUE, parent = emptyenv())
  for (s in seqs) {
    n_k <- nchar(s) - k + 1L
    if (n_k < 1L) next
    for (i in seq_len(n_k)) {
      kmer <- substr(s, i, i + k - 1L)
      cur  <- kmer_env[[kmer]]
      kmer_env[[kmer]] <- if (is.null(cur)) 1L else cur + 1L
    }
  }
  kmer_names <- ls(kmer_env)
  if (length(kmer_names) == 0L) return(fallback())
  coverage <- vapply(kmer_names, function(kn) kmer_env[[kn]], integer(1L))
  names(coverage) <- kmer_names
  kmer_names <- kmer_names[coverage >= min_coverage]
  if (length(kmer_names) == 0L) return(fallback())
  coverage <- coverage[kmer_names]
  prefix_idx <- new.env(hash = TRUE, parent = emptyenv())
  for (kn in kmer_names) {
    px <- substr(kn, 1L, k - 1L)
    bucket <- prefix_idx[[px]]
    prefix_idx[[px]] <- if (is.null(bucket)) kn else c(bucket, kn)
  }
  successors_of <- function(kmer) prefix_idx[[substr(kmer, 2L, k)]]
  indeg_env <- new.env(hash = TRUE, parent = emptyenv())
  for (kn in kmer_names) {
    for (nx in successors_of(kn)) {
      if (identical(nx, kn)) next
      cur <- indeg_env[[nx]]
      indeg_env[[nx]] <- if (is.null(cur)) 1L else cur + 1L
    }
  }
  indeg <- vapply(kmer_names, function(kn) {
    v <- indeg_env[[kn]]; if (is.null(v)) 0L else v
  }, integer(1L))
  names(indeg) <- kmer_names
  zero_indeg <- kmer_names[indeg == 0L]
  seed <- if (length(zero_indeg) > 0L) zero_indeg[which.max(coverage[zero_indeg])]
  else kmer_names[which.max(coverage)]
  path <- seed; visited <- seed; current <- seed
  for (step in seq_len(max_steps)) {
    cand <- setdiff(successors_of(current), visited)
    if (length(cand) == 0L) break
    nxt <- cand[which.max(coverage[cand])]
    path <- c(path, nxt); visited <- c(visited, nxt); current <- nxt
  }
  if (length(path) < 2L) return(path[1])
  paste0(path[1], paste(substr(path[-1], k, k), collapse = ""))
}

#' Simple Overlap Assembly (fallback for short/sparse soft-clip sets)
#' @param seqs Character vector of sequences
#' @param min_overlap Minimum overlap length to merge two sequences
#' @return Assembled contig, or NA
.simple_overlap_assembly <- function(seqs, min_overlap = 10) {
  if (length(seqs) == 0) return(NA_character_)
  seqs <- seqs[!is.na(seqs) & nchar(seqs) >= min_overlap]
  if (length(seqs) == 0) return(NA_character_)
  contig <- seqs[which.max(nchar(seqs))]
  seqs <- seqs[-which.max(nchar(seqs))]
  changed <- TRUE
  while (changed && length(seqs) > 0) {
    changed <- FALSE
    for (i in seq_along(seqs)) {
      s <- seqs[i]
      if (grepl(s, contig, fixed = TRUE)) {
        seqs <- seqs[-i]; changed <- TRUE; break
      }
      for (ov in seq(min_overlap, nchar(s), by = 1)) {
        if (substr(s, 1, ov) == substr(contig, nchar(contig) - ov + 1, nchar(contig))) {
          contig <- paste0(contig, substr(s, ov + 1, nchar(s)))
          seqs <- seqs[-i]; changed <- TRUE; break
        }
        if (substr(contig, 1, ov) == substr(s, nchar(s) - ov + 1, nchar(s))) {
          contig <- paste0(substr(s, 1, nchar(s) - ov), contig)
          seqs <- seqs[-i]; changed <- TRUE; break
        }
      }
      if (changed) break
    }
  }
  contig
}

#' [BUG 6 FIX] Weighted Consensus from Soft-Clip Sequences
#'
#' Per-position weighted majority vote. Weights (e.g. MAPQ) are now actually
#' used: each base vote is weighted by its read's weight. Falls back to
#' unweighted voting when weights are NULL, malformed, or all zero.
#'
#' @param seqs Character vector of sequences
#' @param weights Numeric vector of weights, one per sequence (optional)
#' @return Consensus sequence
.assemble_weighted_consensus <- function(seqs, weights = NULL) {
  if (length(seqs) == 0L) return(NA_character_)
  if (is.null(weights) || length(weights) != length(seqs)) {
    weights <- rep(1, length(seqs))
  }
  weights <- suppressWarnings(as.numeric(weights))
  weights[is.na(weights) | !is.finite(weights) | weights < 0] <- 0
  if (sum(weights) <= 0) weights <- rep(1, length(seqs))
  max_len <- max(nchar(seqs), na.rm = TRUE)
  if (max_len == 0L) return(NA_character_)
  mat <- matrix(NA_character_, nrow = length(seqs), ncol = max_len)
  for (i in seq_along(seqs)) {
    chars <- strsplit(seqs[i], "")[[1]]
    mat[i, seq_len(length(chars))] <- chars
  }
  consensus <- vapply(seq_len(max_len), function(pos) {
    ok  <- !is.na(mat[, pos])
    col <- mat[ok, pos]
    if (length(col) == 0L) return("N")
    w <- weights[ok]
    score <- tapply(w, col, sum)
    names(score)[which.max(score)]
  }, character(1L))
  paste(consensus, collapse = "")
}

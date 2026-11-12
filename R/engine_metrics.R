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
  len_used <- if (length(best_len) != 1L || is.na(best_len)) 0L else best_len
  nominal_read_len <- suppressWarnings(as.numeric(nominal_read_len))
  max_correction <- suppressWarnings(as.numeric(max_correction))
  if (length(nominal_read_len) != 1L || !is.finite(nominal_read_len) || nominal_read_len <= 0) {
    stop("nominal_read_len must be one positive finite number.", call. = FALSE)
  }
  if (length(max_correction) != 1L || !is.finite(max_correction) || max_correction < 1) {
    stop("max_correction must be one finite number >= 1; use 1 to disable size correction.", call. = FALSE)
  }
  size_scalar <- min(1 + (0.25 * len_used / nominal_read_len), max_correction)
  corrected_support <- round(raw_support * size_scalar)
  raw_depth <- raw_support + strict_wt
  corrected_depth <- corrected_support + strict_wt
  raw_af <- if (raw_depth == 0L) 0 else raw_support / raw_depth
  corrected_af <- if (corrected_depth == 0L) 0 else corrected_support / corrected_depth
  list(
    raw_support = as.integer(raw_support),
    corrected_support = as.integer(corrected_support),
    size_correction_factor = as.numeric(size_scalar),
    raw_depth = as.integer(raw_depth),
    corrected_depth = as.integer(corrected_depth),
    raw_af = as.numeric(raw_af),
    corrected_af = as.numeric(corrected_af),
    depth_at_bp = as.integer(raw_depth),
    left_wt = strict_wt,
    right_wt = strict_wt,
    avg_wt = strict_wt
  )
}


# Longest duplication on which sequence-level work (reference imputation,
# pairwise ref / orientation alignment) is attempted.
.MAX_SEQUENCE_OPS_LEN <- 2000L

# Serialize a QNAME set for the result row (SupportQNAMEs / WildtypeQNAMEs).
# .merge_ptd_intervals() splits these strings back into sets and takes the
# exact union, so reads supporting several merged fragments are counted once.
# Oversized sets are truncated with a marker so the merge can fall back to
# the legacy sum instead of computing a wrong union from partial sets.
.qname_set_to_string <- function(qnames, max_chars = 50000L) {
  qnames <- sort(unique(stats::na.omit(as.character(qnames))))
  if (!length(qnames)) return(NA_character_)
  s <- paste(qnames, collapse = ";")
  if (nchar(s) > max_chars) {
    keep <- qnames[cumsum(nchar(qnames) + 1L) <= max_chars]
    s <- paste0(paste(keep, collapse = ";"), ";...TRUNCATED")
  }
  s
}

# Keep one representative alignment row per read before calculating
# row-count-based metrics (soft-clip counts, strand bias, MAPQ, sequence
# consensus, etc.). Support is defined by distinct QNAMEs; duplicate
# evidence rows from one read must not inflate those metrics.
# Build read-level rows for both mates of supporting QNAMEs from the BAM.
# Candidate extraction often retains only the breakpoint-spanning/soft-clipped
# mate. Its paired mate can provide real sequence coverage elsewhere in the
# duplicated segment, but must be aligned as a read sequence (not inferred
# from the pair's outer span or projected genomic coordinates).
.extract_support_mate_coverage_rows <- function(all_reads, support_qnames) {
  empty <- data.frame(
    read_name = character(), read_seq = character(), read_start = integer(),
    read_end = integer(), cigar = character(), is_reverse = logical(),
    stringsAsFactors = FALSE
  )
  if (is.null(all_reads) || length(all_reads) == 0L ||
      is.null(support_qnames) || !length(support_qnames)) return(empty)

  mc <- S4Vectors::mcols(all_reads)
  seqs <- if ("seq" %in% names(mc)) as.character(mc$seq) else character()
  if (length(seqs) != length(all_reads)) return(empty)

  qnames <- names(all_reads)
  if (is.null(qnames) || length(qnames) != length(all_reads) ||
      all(is.na(qnames) | !nzchar(qnames))) {
    qnames <- if ("qname" %in% names(mc)) as.character(mc$qname) else rep(NA_character_, length(all_reads))
  }
  qnames <- .normalize_read_qname(as.character(qnames))
  support_norm <- unique(stats::na.omit(.normalize_read_qname(as.character(support_qnames))))
  keep <- !is.na(qnames) & qnames %in% support_norm &
    !is.na(seqs) & nzchar(seqs)
  if (!any(keep)) return(empty)

  flags <- if ("flag" %in% names(mc)) as.integer(mc$flag) else rep(0L, length(all_reads))
  flags[is.na(flags)] <- 0L
  cigars <- if ("cigar" %in% names(mc)) as.character(mc$cigar) else rep(NA_character_, length(all_reads))
  data.frame(
    read_name = qnames[keep], read_seq = toupper(seqs[keep]),
    read_start = as.integer(GenomicRanges::start(all_reads)[keep]),
    read_end = as.integer(GenomicRanges::end(all_reads)[keep]),
    cigar = cigars[keep], is_reverse = bitwAnd(flags[keep], 0x10L) != 0L,
    stringsAsFactors = FALSE
  )
}

# Keep one representative alignment row per read before calculating
# row-count-based metrics. Distinct mates remain separate observations.
.deduplicate_coverage_read_rows <- function(read_rows) {
  if (is.null(read_rows) || !is.data.frame(read_rows) || nrow(read_rows) == 0L ||
      !all(c("read_name", "read_seq") %in% names(read_rows))) {
    return(read_rows)
  }
  qname <- .normalize_read_qname(as.character(read_rows$read_name))
  read_seq <- toupper(as.character(read_rows$read_seq))
  keep <- !is.na(qname) & nzchar(qname) & !is.na(read_seq) & nzchar(read_seq)
  read_rows <- read_rows[keep, , drop = FALSE]
  qname <- qname[keep]
  read_seq <- read_seq[keep]
  if (nrow(read_rows) < 2L) return(read_rows)

  # Candidate extraction can emit multiple rows for the same alignment when
  # one read supports more than one nearby seed/anchor. Collapse only identical
  # read-end alignments, not all rows sharing a QNAME: R1 and R2 are distinct
  # observations for per-base read coverage.
  signature_cols <- intersect(
    c("read_start", "read_end", "cigar", "is_reverse"),
    names(read_rows)
  )
  parts <- c(list(qname), lapply(signature_cols, function(nm) {
    x <- as.character(read_rows[[nm]])
    x[is.na(x)] <- "<NA>"
    x
  }), list(read_seq))
  key <- do.call(paste, c(parts, sep = "\034"))
  read_rows[!duplicated(key), , drop = FALSE]
}

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

#' PE reference-coordinate projection diagnostic for the ITD interval
#'
#' For every valid pair whose QNAME is among the supporting reads, each mate
#' that falls entirely inside the duplicated reference segment
#' [breakpoint, breakpoint + span_len - 1] contributes coverage over its own
#' extent, mapped into ITD-sequence coordinates (position i of the ITD
#' sequence corresponds to reference position breakpoint + i - 1). Mate
#' intervals are unioned per fragment, so an overlapping base is counted
#' once while non-overlapping bases from both mates are retained. This is
#' strictly less specific than the direct PID >= 90 alignment coverage: a
#' wild-type fragment whose mate happens to land inside the segment also
#' counts. It is therefore reported in separate columns
#' (ITDCoverageRLE_PE / ITDReadCoverage_PE / ITDPE_SpanningPairs); PE coordinate
#' projection is diagnostic and is excluded from combined allele coverage.
#' and the min_itd_read_coverage filter keeps operating on the direct
#' metric.
#'
#' FF / RR pairs are skipped: they are not normal spanning fragments.
#' RF (everted) pairs are kept - they are the tandem-duplication junction
#' signature and their inside-segment mate is still a real fragment.
#'
#' Pure (plain vectors in, list out) so it can be unit-tested without any
#' BAM / GAlignmentPairs machinery.
#'
#' @param qname,s1,e1,s2,e2,orientation,valid Per-pair vectors from
#'   \code{.pair_info()} (starts/ends of both mates, pair orientation,
#'   primary + MAPQ validity).
#' @param support_qnames Post-dedup QNAMEs of the supporting reads (one row
#'   per pair in GAlignmentPairs, so no read is counted twice).
#' @param breakpoint Leftmost genomic base of the duplicated segment.
#' @param span_len Duplication length (bp).
#' @param seq_len Length of the reconstructed ITD sequence (the vector the
#'   coverage is expressed over; may be shorter than span_len, in which
#'   case contributions are clamped).
#' @param direct_cov Optional direct-alignment coverage vector (length
#'   seq_len) used only to compute a legacy diagnostic union percentage; this
#'   value is not used for combined allele-specific coverage.
#' @return list(cov_vec, rle_str, percent, percent_combined,
#'   n_spanning_pairs); percent is PE-only. percent_combined is retained as a
#'   legacy geometric-union diagnostic and is not used for allele coverage.
#' @keywords internal
.pe_spanning_coverage <- function(
    qname,
    s1, e1, s2, e2,
    orientation,
    valid,
    support_qnames,
    breakpoint,
    span_len,
    seq_len,
    direct_cov = integer(0)
) {
  empty <- list(
    cov_vec = integer(0),
    rle_str = NA_character_,
    percent = NA_real_,
    percent_combined = NA_real_,
    n_spanning_pairs = 0L
  )
  if (is.null(qname) || length(qname) == 0L) {
    return(empty)
  }
  if (length(breakpoint) != 1L || is.na(breakpoint) ||
      length(span_len) != 1L || is.na(span_len) || span_len <= 0L ||
      length(seq_len) != 1L || is.na(seq_len) || seq_len <= 0L) {
    return(empty)
  }
  if (!is.null(direct_cov) && length(direct_cov) > 0L &&
      length(direct_cov) != seq_len) {
    stop("direct_cov must have length seq_len (or be empty).", call. = FALSE)
  }
  n <- length(qname)
  stopifnot(
    length(s1) == n, length(e1) == n, length(s2) == n, length(e2) == n,
    length(orientation) == n, length(valid) == n
  )

  span_lo <- as.integer(breakpoint)
  span_hi <- span_lo + as.integer(span_len) - 1L

  # Normalize both sources before matching. BAM query names and read names
  # may differ only by the conventional /1 or /2 suffix or whitespace tags.
  qname <- .normalize_read_qname(qname)
  support_qnames <- unique(stats::na.omit(.normalize_read_qname(support_qnames)))
  keep <- valid & !is.na(orientation) & !is.na(qname) & qname %in% support_qnames
  keep[is.na(keep)] <- FALSE
  keep <- which(keep)
  keep <- keep[orientation[keep] %in% c("FR", "RF")]

  interval_starts <- integer(0)
  interval_ends <- integer(0)
  n_used <- 0L
  for (ki in keep) {
    # Convert each fragment's two mate intervals to a disjoint union. The
    # per-fragment loop is retained for pair geometry, but depth accumulation
    # is deferred to one difference-array/tabulation pass below.
    lo <- c(s1[ki], s2[ki])
    hi <- c(e1[ki], e2[ki])
    inside <- is.finite(lo) & is.finite(hi) & lo >= span_lo &
      hi <= span_hi & hi >= lo
    lo <- pmax(1L, as.integer(lo[inside] - span_lo + 1L))
    hi <- pmin(as.integer(seq_len), as.integer(hi[inside] - span_lo + 1L))
    keep_interval <- hi >= lo
    lo <- lo[keep_interval]
    hi <- hi[keep_interval]
    if (!length(lo)) next
    ord <- order(lo, hi)
    lo <- lo[ord]; hi <- hi[ord]
    if (length(lo) == 2L && lo[2L] <= hi[1L] + 1L) {
      lo <- lo[1L]
      hi <- max(hi)
    }
    interval_starts <- c(interval_starts, lo)
    interval_ends <- c(interval_ends, hi)
    n_used <- n_used + 1L
  }
  # Difference-array accumulation counts each fragment once per covered base,
  # including overlapping mates, without allocating a template-length logical
  # vector or updating every covered base inside the pair loop.
  events <- tabulate(interval_starts, nbins = as.integer(seq_len) + 1L) -
    tabulate(interval_ends + 1L, nbins = as.integer(seq_len) + 1L)
  cov_vec <- as.integer(cumsum(events)[seq_len])
  rle_res <- rle(cov_vec)
  rle_str <- paste(paste0(rle_res$values, ":", rle_res$lengths), collapse = ",")
  pe_covered <- sum(cov_vec > 0L)
  union_covered <- if (length(direct_cov) == seq_len) {
    sum(cov_vec > 0L | direct_cov > 0L)
  } else {
    pe_covered
  }
  list(
    cov_vec = cov_vec,
    rle_str = rle_str,
    percent = pe_covered / seq_len * 100,
    percent_combined = union_covered / seq_len * 100,
    n_spanning_pairs = n_used
  )
}

#' Fold direct read-alignment / k-mer coverage into one display vector
#'
#' The displayed vector takes, per base, the per-base maximum across direct
#' read-sequence alignments and sequence-anchored k-mer evidence. PE
#' reference-coordinate projections are deliberately excluded because they
#' cannot distinguish the original WT copy from the duplicated ITD copy. This is a
#' non-additive depth envelope: methods can observe the same fragment, so
#' summing them would double-count evidence. Taking the maximum avoids that
#' duplication while ensuring the combined total is not lower than any one
#' contributing source. All inputs are optional (pass \code{integer(0)});
#' vectors of the wrong length are ignored.
#'
#' @param cov_direct Integer vector - direct-alignment per-base counts.
#' @param cov_pe     Legacy optional vector; caller passes an empty vector for allele-specific coverage.
#' @param cov_kmer   Integer vector - k-mer per-base counts.
#' @param direct_is_binary Logical - TRUE only when \code{cov_direct} comes
#'   from the no-pwalign fallback (0/1 "covered" flags, not read counts).
#'   The flag is supplied by the caller; binary direct flags are normalized
#'   to 0/1 before taking the per-base maximum.
#' @return Integer vector (length = the common vector length), or
#'   \code{integer(0)} when no input has a valid (non-zero, mutually
#'   consistent) length.
#' @keywords internal
.fold_display_coverage <- function(cov_direct, cov_pe, cov_kmer,
                                   direct_is_binary = FALSE) {
  # The combined vector is a per-position depth envelope, not a sum of
  # method-specific depths: the same fragment/read can contribute to more
  # than one method, so summing would double-count evidence. Conversely,
  # filling only uncovered positions under-reports depth whenever a less
  # specific method has greater depth at an already-covered position.
  lens <- c(
    if (is.null(cov_direct) || length(cov_direct) == 0L) 0L else length(cov_direct),
    if (is.null(cov_pe) || length(cov_pe) == 0L) 0L else length(cov_pe),
    if (is.null(cov_kmer) || length(cov_kmer) == 0L) 0L else length(cov_kmer)
  )
  n <- max(lens)
  if (n == 0L) return(integer(0))
  keep <- lens == n
  sources <- list(cov_direct, cov_pe, cov_kmer)
  out <- integer(n)
  for (i in seq_along(sources)) {
    if (!keep[i]) next
    x <- as.integer(sources[[i]])
    if (i == 1L && isTRUE(direct_is_binary)) x <- as.integer(x > 0L)
    out <- pmax(out, x)
  }
  out
}

#' Position-resolved k-mer coverage of the ITD sequence
#'
#' K-mers are only allowed to contribute when they identify one position in
#' the reconstructed ITD sequence. A repeated k-mer cannot tell which copy
#' of a tandem repeat was observed; assigning it to every occurrence
#' overstates per-position coverage. Hits are also required to agree on a
#' common read-to-ITD offset (at least two overlapping k-mers by default).
#' Reads are evaluated in both orientations, but each read contributes at
#' most one count to each ITD base.
#'
#' @param itd_seq Character sequence of the reconstructed ITD.
#' @param read_seqs Character vector of supporting read sequences.
#' @param k K-mer size.
#' @param min_offset_hits Minimum number of k-mer hits supporting one
#'   collinear offset before that read contributes coverage.
#' @return Per-base coverage and an RLE summary.
#' @keywords internal
.kmer_itd_coverage <- function(itd_seq, read_seqs, k = 11L,
                               min_offset_hits = 2L) {
  empty <- list(
    cov_vec = integer(0),
    rle_str = NA_character_,
    percent = 0,
    n_reads = 0L,
    matched_reads = logical(0),
    per_read_hits = list(),
    ambiguous_cov_vec = integer(0),
    ambiguous_percent = 0
  )
  if (!is.numeric(k) || !is.numeric(min_offset_hits) ||
      length(k) != 1L || length(min_offset_hits) != 1L ||
      is.na(k) || is.na(min_offset_hits) ||
      !is.finite(k) || !is.finite(min_offset_hits) ||
      k < 1L || min_offset_hits < 1L ||
      k != as.integer(k) || min_offset_hits != as.integer(min_offset_hits)) {
    stop("k and min_offset_hits must be positive integers.", call. = FALSE)
  }
  k <- as.integer(k)
  min_offset_hits <- as.integer(min_offset_hits)
  if (length(itd_seq) != 1L || is.na(itd_seq) || !nzchar(itd_seq)) {
    return(empty)
  }
  itd_length <- nchar(itd_seq)
  if (itd_length < k || length(read_seqs) == 0L) {
    empty$matched_reads <- rep(FALSE, length(read_seqs))
    empty$per_read_hits <- rep(list(logical(max(0L, itd_length))), length(read_seqs))
    if (itd_length < k) empty$cov_vec <- integer(itd_length)
    return(empty)
  }

  itd_seq <- toupper(itd_seq)
  cache_key <- paste0("kmer|", k, "|", itd_length, "|", itd_seq)
  cached_index <- .talos_cache_get(cache_key)
  if (is.null(cached_index)) {
    n_km <- itd_length - k + 1L
    itd_kmers <- substring(itd_seq, seq_len(n_km), seq_len(n_km) + k - 1L)
    kmer_positions <- split(seq_len(n_km), itd_kmers)
    # A repeated k-mer is not position-specific and cannot establish coverage
    # at every copy. Cache only the unique index, not read-derived evidence.
    unique_names <- names(kmer_positions)[vapply(
      kmer_positions, function(pos) length(pos) == 1L, logical(1L))
    ]
    unique_positions <- if (length(unique_names)) {
      vapply(kmer_positions[unique_names], function(x) as.integer(x[1L]), integer(1L))
    } else integer(0)
    names(unique_positions) <- unique_names
    cached_index <- unique_positions
    .talos_cache_set(cache_key, cached_index)
  }
  unique_positions <- cached_index
  if (length(unique_positions) == 0L) {
    cov_vec <- integer(itd_length)
    return(list(cov_vec = cov_vec, rle_str = paste0("0:", itd_length),
                percent = 0, n_reads = 0L,
                matched_reads = rep(FALSE, length(read_seqs)),
                per_read_hits = rep(list(logical(itd_length)), length(read_seqs))))
  }
  per_read_ambiguous_hits <- rep(list(logical(itd_length)), length(read_seqs))
  matched_reads <- rep(FALSE, length(read_seqs))
  per_read_hits <- rep(list(logical(itd_length)), length(read_seqs))
  n_reads_with_matches <- 0L
  for (read_idx in seq_along(read_seqs)) {
    read_seq <- read_seqs[[read_idx]]
    if (length(read_seq) != 1L || is.na(read_seq) || nchar(read_seq) < k) next
    read_seq <- toupper(read_seq)
    oriented_reads <- read_seq
    if (requireNamespace("Biostrings", quietly = TRUE)) {
      rc <- tryCatch(
        as.character(Biostrings::reverseComplement(Biostrings::DNAString(read_seq))),
        error = function(e) NA_character_
      )
      if (!is.na(rc) && !identical(rc, read_seq)) oriented_reads <- c(oriented_reads, rc)
    }

    read_hit <- logical(itd_length)
    read_ambiguous_hit <- logical(itd_length)
    for (oriented_seq in oriented_reads) {
      read_len <- nchar(oriented_seq)
      read_starts <- seq_len(read_len - k + 1L)
      read_kmers <- substring(oriented_seq, read_starts, read_starts + k - 1L)
      hit_idx <- match(read_kmers, names(unique_positions))
      keep <- which(!is.na(hit_idx))
      if (length(keep) < min_offset_hits) next

      target_starts <- unname(unique_positions[hit_idx[keep]])
      offsets <- target_starts - read_starts[keep]
      offset_counts <- table(offsets)
      accepted_offsets <- as.integer(names(offset_counts)[offset_counts >= min_offset_hits])
      if (!length(accepted_offsets)) next

      # A read with multiple qualifying offsets is positionally ambiguous.
      # Never union those placements into the certain k-mer coverage vector.
      # Retain the possible union separately for diagnostics/reporting.
      if (length(accepted_offsets) != 1L) {
        for (offset in accepted_offsets) {
          starts <- target_starts[offsets == offset]
          for (pos in starts) {
            end_pos <- pos + k - 1L
            if (pos >= 1L && end_pos <= itd_length) read_ambiguous_hit[pos:end_pos] <- TRUE
          }
        }
        next
      }
      starts <- target_starts[offsets == accepted_offsets[[1L]]]
      for (pos in starts) {
        end_pos <- pos + k - 1L
        if (pos >= 1L && end_pos <= itd_length) read_hit[pos:end_pos] <- TRUE
      }
    }

    per_read_hits[[read_idx]] <- read_hit
    per_read_ambiguous_hits[[read_idx]] <- read_ambiguous_hit
    if (any(read_hit)) {
      matched_reads[read_idx] <- TRUE
      n_reads_with_matches <- n_reads_with_matches + 1L
    }
  }

  # Preserve each read's positional mask, then aggregate all positions in one
  # compiled tabulation pass rather than updating depth in the read loop.
  flatten_mask_positions <- function(x) {
    unlist(lapply(x, which), use.names = FALSE)
  }
  cov_vec <- tabulate(flatten_mask_positions(per_read_hits), nbins = itd_length)
  ambiguous_cov_vec <- tabulate(flatten_mask_positions(per_read_ambiguous_hits),
                                nbins = itd_length)
  rle_res <- rle(cov_vec)
  rle_str <- paste(paste0(rle_res$values, ":", rle_res$lengths), collapse = ",")
  percent <- sum(cov_vec > 0L) / itd_length * 100
  list(cov_vec = cov_vec, rle_str = rle_str, percent = percent,
       n_reads = n_reads_with_matches, matched_reads = matched_reads,
       per_read_hits = per_read_hits,
       ambiguous_cov_vec = ambiguous_cov_vec,
       ambiguous_percent = sum(ambiguous_cov_vec > 0L) / itd_length * 100)
}

#' Multi-k, all-offset k-mer coverage for tandem-duplication templates.
#'
#' Each k is scanned at every possible nucleotide offset (not in a fixed
#' reading frame). Only template-unique k-mers are used, and anchors within a
#' read must agree on a shared read-to-template offset. The combined vector
#' unions each read end across all k values before adding depth, so one read
#' can contribute at most +1 to a base even when several k sizes support it.
#' Shorter k values require more independent collinear anchors to control
#' chance matches in repetitive DNA.
#'
#' @param itd_seq ITD coordinate template.
#' @param read_seqs Supporting read-end sequences.
#' @param k_values K-mer sizes to evaluate; each scans every possible offset.
#' @return Combined and per-k coverage summaries.
#' @keywords internal
.multi_kmer_itd_coverage <- function(
    itd_seq,
    read_seqs,
    k_values = c(7L, 9L, 11L, 13L, 15L),
    min_offset_hits = c(3L, 2L, 2L, 2L, 2L)
) {
  if (!is.numeric(k_values) || length(k_values) == 0L || anyNA(k_values) ||
      any(!is.finite(k_values)) || any(k_values < 1L) ||
      any(k_values != as.integer(k_values))) {
    stop("k_values must contain positive integers.", call. = FALSE)
  }
  k_values <- as.integer(k_values)
  if (anyDuplicated(k_values)) {
    stop("k_values must not contain duplicates.", call. = FALSE)
  }
  if (!is.numeric(min_offset_hits)) {
    stop("min_offset_hits must be a positive integer per k value.", call. = FALSE)
  }
  if (length(min_offset_hits) == 1L) min_offset_hits <- rep(min_offset_hits, length(k_values))
  if (length(min_offset_hits) != length(k_values) || anyNA(min_offset_hits) ||
      any(!is.finite(min_offset_hits)) || any(min_offset_hits < 1L) ||
      any(min_offset_hits != as.integer(min_offset_hits))) {
    stop("min_offset_hits must be a positive integer per k value.", call. = FALSE)
  }
  min_offset_hits <- as.integer(min_offset_hits)
  n <- if (length(itd_seq) == 1L && !is.na(itd_seq)) nchar(itd_seq) else 0L
  if (n == 0L) {
    return(list(cov_vec = integer(0), rle_str = NA_character_, percent = NA_real_,
                n_reads = 0L, matched_reads = rep(FALSE, length(read_seqs)),
                per_k = list(), k_values = k_values,
                ambiguous_cov_vec = integer(0), ambiguous_percent = NA_real_))
  }
  per_k <- lapply(seq_along(k_values), function(i) {
    z <- .kmer_itd_coverage(
      itd_seq, read_seqs, k = k_values[i],
      min_offset_hits = min_offset_hits[i]
    )
    if (length(z$cov_vec) != n) {
      z$cov_vec <- integer(n)
      z$rle_str <- paste0("0:", n)
      z$percent <- 0
      z$n_reads <- 0L
      z$matched_reads <- rep(FALSE, length(read_seqs))
      z$per_read_hits <- rep(list(logical(n)), length(read_seqs))
      z$ambiguous_cov_vec <- integer(n)
      z$ambiguous_percent <- 0
    }
    z
  })
  names(per_k) <- as.character(k_values)

  combined <- integer(n)
  ambiguous_by_k <- lapply(per_k, function(z) {
    if (length(z$ambiguous_cov_vec) == n) as.integer(z$ambiguous_cov_vec) else integer(n)
  })
  ambiguous_combined <- if (length(ambiguous_by_k)) Reduce(pmax, ambiguous_by_k) else integer(n)
  matched_reads <- rep(FALSE, length(read_seqs))
  for (ri in seq_along(read_seqs)) {
    per_read <- lapply(per_k, function(z) {
      if (length(z$per_read_hits) >= ri) z$per_read_hits[[ri]] else logical(n)
    })
    if (length(per_read)) {
      hit <- Reduce(`|`, per_read)
      if (length(hit) == n && any(hit)) {
        combined[hit] <- combined[hit] + 1L
        matched_reads[ri] <- TRUE
      }
    }
  }
  rr <- rle(combined)
  list(
    cov_vec = combined,
    rle_str = paste(paste0(rr$values, ":", rr$lengths), collapse = ","),
    percent = sum(combined > 0L) / n * 100,
    n_reads = as.integer(sum(matched_reads)),
    matched_reads = matched_reads,
    per_k = per_k,
    k_values = k_values,
    ambiguous_cov_vec = ambiguous_combined,
    ambiguous_percent = sum(ambiguous_combined > 0L) / n * 100
  )
}

# Align every distinct supporting read-end sequence to an ITD coordinate
# template in bounded chunks. Coverage must never be calculated from the small
# subsample used by older versions for optional consistency scoring: that made
# the percentage depend on max_pairwise_alignments (30 by default) rather than
# on the sample. Both sequence orientations are tested, and each input read end
# contributes at most one selected alignment. The caller separately tracks
# distinct QNAME fragments so R1/R2 can both contribute coverage without
# inflating fragment-level support counts. Long, high-identity alignments are
# preferred over short incidental exact matches.
# Small bounded cache for immutable template-derived structures. This avoids
# rebuilding unique k-mer indexes and tandem-junction windows for every sample
# or candidate that shares the same reference template. The FIFO bound prevents
# long cohort runs from retaining every candidate sequence indefinitely.
.TALOS_TEMPLATE_CACHE <- new.env(parent = emptyenv())
.TALOS_TEMPLATE_CACHE$keys <- character(0)
.TALOS_TEMPLATE_CACHE$max_entries <- 24L

.talos_cache_get <- function(key) {
  if (exists(key, envir = .TALOS_TEMPLATE_CACHE, inherits = FALSE)) {
    .TALOS_TEMPLATE_CACHE$keys <- c(key, .TALOS_TEMPLATE_CACHE$keys[
      .TALOS_TEMPLATE_CACHE$keys != key])
    return(get(key, envir = .TALOS_TEMPLATE_CACHE, inherits = FALSE))
  }
  NULL
}

.talos_cache_set <- function(key, value) {
  if (exists(key, envir = .TALOS_TEMPLATE_CACHE, inherits = FALSE)) {
    rm(list = key, envir = .TALOS_TEMPLATE_CACHE)
  }
  assign(key, value, envir = .TALOS_TEMPLATE_CACHE)
  .TALOS_TEMPLATE_CACHE$keys <- c(key, .TALOS_TEMPLATE_CACHE$keys[
    .TALOS_TEMPLATE_CACHE$keys != key])
  limit <- max(1L, as.integer(.TALOS_TEMPLATE_CACHE$max_entries)[1L])
  while (length(.TALOS_TEMPLATE_CACHE$keys) > limit) {
    old <- tail(.TALOS_TEMPLATE_CACHE$keys, 1L)
    if (exists(old, envir = .TALOS_TEMPLATE_CACHE, inherits = FALSE))
      rm(list = old, envir = .TALOS_TEMPLATE_CACHE)
    .TALOS_TEMPLATE_CACHE$keys <- head(.TALOS_TEMPLATE_CACHE$keys, -1L)
  }
  invisible(value)
}

.template_unique_seed_map <- function(template, seed_k) {
  template <- toupper(as.character(template)[1L])
  seed_k <- max(1L, as.integer(seed_k)[1L])
  key <- paste0("seed|", seed_k, "|", nchar(template), "|", template)
  cached <- .talos_cache_get(key)
  if (!is.null(cached)) return(cached)
  n <- nchar(template)
  out <- integer(0)
  if (n >= seed_k) {
    starts <- seq_len(n - seed_k + 1L)
    seeds <- substring(template, starts, starts + seed_k - 1L)
    counts <- table(seeds)
    unique_names <- names(counts)[counts == 1L]
    out <- setNames(starts[match(unique_names, seeds)], unique_names)
  }
  .talos_cache_set(key, out)
  out
}

.template_junction_window <- function(template, anchor_width) {
  template <- toupper(as.character(template)[1L])
  n <- nchar(template)
  width <- max(1L, min(as.integer(anchor_width)[1L], n))
  key <- paste0("junction|", width, "|", n, "|", template)
  cached <- .talos_cache_get(key)
  if (!is.null(cached)) return(cached)
  out <- paste0(substring(template, n - width + 1L, n),
                substring(template, 1L, width))
  .talos_cache_set(key, out)
  out
}

.align_reads_to_itd_template <- function(
    read_seqs,
    template,
    chunk_size = 100L,
    min_pid = 90,
    min_aligned = 20L
) {
  if (length(template) != 1L || is.na(template) || !nzchar(template)) {
    stop("'template' must be one non-empty ITD sequence.", call. = FALSE)
  }
  template <- toupper(as.character(template))
  n <- nchar(template)
  empty <- list(
    cov_vec = integer(n),
    rle_str = paste0("0:", n),
    percent = 0,
    reads_processed = 0L,
    reads_matched = 0L,
    pid_by_read = numeric(0),
    score_by_read = numeric(0),
    aligned_query_width = integer(0),
    aligned_query_start = integer(0),
    alignment_orientation = character(0),
    aligned_subject_start = integer(0),
    aligned_subject_end = integer(0),
    query_coverage_pct = numeric(0),
    per_read_positions = list()
  )
  if (!length(read_seqs)) return(empty)
  if (!requireNamespace("pwalign", quietly = TRUE) ||
      !requireNamespace("Biostrings", quietly = TRUE)) {
    stop(
      "ITD coverage requires the imported pwalign and Biostrings packages.",
      call. = FALSE
    )
  }
  if (length(chunk_size) != 1L || is.na(chunk_size) ||
      !is.finite(chunk_size) || chunk_size < 1L) {
    stop("'chunk_size' must be a positive integer.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)
  min_aligned <- as.integer(min_aligned)

  read_seqs <- toupper(as.character(read_seqs))
  valid <- !is.na(read_seqs) & nzchar(read_seqs)
  read_seqs <- read_seqs[valid]
  if (!length(read_seqs)) return(empty)

  pid_by_read <- rep(NA_real_, length(read_seqs))
  score_by_read <- rep(NA_real_, length(read_seqs))
  aligned_query_width <- integer(length(read_seqs))
  aligned_query_start <- rep(NA_integer_, length(read_seqs))
  alignment_orientation <- rep(NA_character_, length(read_seqs))
  aligned_subject_start <- rep(NA_integer_, length(read_seqs))
  aligned_subject_end <- rep(NA_integer_, length(read_seqs))
  query_coverage_pct <- rep(NA_real_, length(read_seqs))
  per_read_positions <- rep(list(integer(0)), length(read_seqs))
  read_starts <- seq.int(1L, length(read_seqs), by = chunk_size)

  for (chunk_start in read_starts) {
    chunk_end <- min(length(read_seqs), chunk_start + chunk_size - 1L)
    chunk <- read_seqs[chunk_start:chunk_end]
    n_chunk <- length(chunk)

    rc <- tryCatch(
      as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(chunk))),
      error = function(e) chunk
    )
    oriented <- c(chunk, rc)
    read_index <- rep(seq_len(n_chunk), 2L)
    aln <- pwalign::pairwiseAlignment(
      pattern = Biostrings::DNAStringSet(oriented),
      subject = Biostrings::DNAString(template),
      type = "local"
    )

    pid_all <- suppressWarnings(as.numeric(pwalign::pid(aln)))
    pattern_ranges <- pwalign::pattern(aln)
    subject_ranges <- pwalign::subject(aln)
    aligned_pattern_all <- as.character(pwalign::alignedPattern(aln))
    aligned_subject_all <- as.character(pwalign::alignedSubject(aln))
    query_width_all <- suppressWarnings(as.integer(
      BiocGenerics::width(pattern_ranges@range)
    ))
    query_start_all <- suppressWarnings(as.integer(
      BiocGenerics::start(pattern_ranges@range)
    ))
    subject_start_all <- suppressWarnings(as.integer(
      BiocGenerics::start(subject_ranges)
    ))
    subject_end_all <- suppressWarnings(as.integer(
      BiocGenerics::end(subject_ranges)
    ))
    subject_width_all <- subject_end_all - subject_start_all + 1L

    for (read_i in seq_len(n_chunk)) {
      candidates <- which(
        read_index == read_i & is.finite(pid_all) &
          !is.na(query_width_all) & query_width_all >= min_aligned &
          !is.na(subject_width_all) & subject_width_all >= min_aligned &
          !is.na(subject_start_all) & !is.na(subject_end_all) &
          subject_start_all >= 1L & subject_end_all <= n &
          subject_end_all >= subject_start_all
      )
      if (!length(candidates)) next

      # Prefer the longest query span first. PID breaks ties, preventing a
      # short exact repeat from replacing a much more informative alignment.
      # A longer alignment below the identity threshold must not suppress a
      # shorter, genuinely qualifying alignment from the other read orientation.
      # Within each class, retain the longest aligned query/subject span, then PID.
      best <- candidates[order(
        -(pid_all[candidates] >= min_pid),
        -query_width_all[candidates],
        -subject_width_all[candidates],
        -pid_all[candidates],
        candidates
      )][1L]
      global_i <- chunk_start + read_i - 1L
      pid_by_read[global_i] <- pid_all[best]
      score_by_read[global_i] <- suppressWarnings(as.numeric(BiocGenerics::score(aln))[best])
      aligned_query_width[global_i] <- query_width_all[best]
      aligned_query_start[global_i] <- query_start_all[best]
      alignment_orientation[global_i] <- if (best <= n_chunk) "forward" else "reverse"
      aligned_subject_start[global_i] <- subject_start_all[best]
      aligned_subject_end[global_i] <- subject_end_all[best]
      query_coverage_pct[global_i] <- 100 * query_width_all[best] / nchar(chunk[read_i])

      if (pid_all[best] >= min_pid) {
        # Mark only reference bases aligned to an actual read base. A deletion
        # in the read (pattern gap) is not read coverage and must remain zero.
        ap <- strsplit(aligned_pattern_all[best], "", fixed = TRUE)[[1L]]
        asub <- strsplit(aligned_subject_all[best], "", fixed = TRUE)[[1L]]
        if (length(ap) == length(asub) && length(asub)) {
          subject_has_base <- asub != "-"
          pattern_has_base <- ap != "-"
          subject_positions <- subject_start_all[best] +
            cumsum(subject_has_base) - 1L
          covered_positions <- subject_positions[
            subject_has_base & pattern_has_base
          ]
          covered_positions <- covered_positions[
            !is.na(covered_positions) & covered_positions >= 1L &
              covered_positions <= n
          ]
          if (length(covered_positions)) {
            covered_positions <- unique(as.integer(covered_positions))
            per_read_positions[[global_i]] <- covered_positions
          }
        }
      }
    }
  }

  # The read-end placement list is the source of truth. Aggregate all unique
  # read-end positions at once after alignment, avoiding an R-level depth
  # update for every aligned read while retaining the complete read-level map.
  flat_positions <- unlist(lapply(per_read_positions, unique), use.names = FALSE)
  cov_vec <- tabulate(flat_positions[!is.na(flat_positions) &
    flat_positions >= 1L & flat_positions <= n], nbins = n)
  rle_res <- rle(cov_vec)
  list(
    cov_vec = cov_vec,
    rle_str = paste(paste0(rle_res$values, ":", rle_res$lengths), collapse = ","),
    percent = sum(cov_vec > 0L) / n * 100,
    reads_processed = as.integer(length(read_seqs)),
    reads_matched = as.integer(sum(is.finite(pid_by_read) & pid_by_read >= min_pid)),
    pid_by_read = pid_by_read,
    score_by_read = score_by_read,
    aligned_query_width = aligned_query_width,
    aligned_query_start = aligned_query_start,
    alignment_orientation = alignment_orientation,
    aligned_subject_start = aligned_subject_start,
    aligned_subject_end = aligned_subject_end,
    query_coverage_pct = query_coverage_pct,
    per_read_positions = per_read_positions
  )
}


# Enumerate co-optimal local placements by masking each discovered template
# interval and realigning. This is a bounded, repeat-aware alternative search:
# it does not claim to enumerate every suboptimal alignment in the genome, but
# it identifies equivalent placements within the ITD template. Per read-end,
# the intersection across plausible placements is the conservative coverage
# and the union is the possible-coverage upper bound.
.align_reads_to_itd_template_ambiguity <- function(
    read_seqs,
    template,
    min_pid = 90,
    min_aligned = 20L,
    max_alternatives = 8L,
    score_delta = 1,
    best_result = NULL,
    seed_k = 9L,
    progress_every = 25L,
    try_reverse = TRUE
) {
  template <- toupper(as.character(template)[1L])
  n <- nchar(template)
  empty <- list(
    certain_cov = integer(n), possible_cov = integer(n),
    ambiguous_cov = integer(n), certain_percent = 0, possible_percent = 0,
    ambiguous_percent = 0, ambiguous_read_ends = 0L, multi_placement_read_ends = 0L,
    alternative_search_limit_read_ends = 0L,
    reads_with_qualifying_alignment = 0L, mean_placements = NA_real_,
    max_placements = 0L, status = "computed_no_qualifying_alignments"
  )
  if (!length(read_seqs) || n < min_aligned) return(empty)
  if (!requireNamespace("pwalign", quietly = TRUE) ||
      !requireNamespace("Biostrings", quietly = TRUE)) {
    empty$status <- "unavailable_alignment_packages"
    return(empty)
  }
  max_alternatives <- max(1L, as.integer(max_alternatives)[1L])
  min_aligned <- max(1L, as.integer(min_aligned)[1L])
  score_delta <- max(0, as.numeric(score_delta)[1L])
  seed_k <- max(5L, as.integer(seed_k)[1L])
  has_best_positions <- is.list(best_result) &&
    length(best_result$per_read_positions) == length(read_seqs)

  # Fast-path reads whose sequence has no exact seed shared by multiple
  # template positions. Their primary local alignment is an unambiguous
  # positional placement; expensive masked re-alignment is reserved for reads
  # that could genuinely have multiple repeat placements.
  unique_seed_map <- .template_unique_seed_map(template, seed_k)
  collinear_unique_seed_offsets <- function(seq) {
    if (!length(unique_seed_map) || nchar(seq) < seed_k) return(character(0))
    q_starts <- seq_len(nchar(seq) - seed_k + 1L)
    q_seeds <- substring(seq, q_starts, q_starts + seed_k - 1L)
    target_starts <- unname(unique_seed_map[q_seeds])
    keep <- which(!is.na(target_starts))
    if (length(keep) < 2L) return(character(0))
    offsets <- target_starts[keep] - q_starts[keep]
    names(which(table(offsets) >= 2L))
  }
  read_seqs <- toupper(as.character(read_seqs))
  good_read <- !is.na(read_seqs) & nzchar(read_seqs)
  # Retain read-end-level positions so later evidence sources can be merged
  # without either double-counting the same read or undercounting distinct reads.
  per_read_certain_positions <- vector("list", length(read_seqs))
  per_read_possible_positions <- vector("list", length(read_seqs))
  n_ambiguous <- 0L
  n_multi <- 0L
  n_search_limit <- 0L
  n_qualifying <- 0L
  placement_counts <- integer(0)
  progress_every <- suppressWarnings(as.integer(progress_every)[1L])
  if (is.na(progress_every) || progress_every < 1L) progress_every <- 25L
  read_indices <- which(good_read)
  n_good_reads <- length(read_indices)

  alignment_positions <- function(aln, subject_start) {
    ap <- as.character(pwalign::alignedPattern(aln))[1L]
    asub <- as.character(pwalign::alignedSubject(aln))[1L]
    if (is.na(ap) || is.na(asub) || !nzchar(ap) || nchar(ap) != nchar(asub)) {
      return(integer(0))
    }
    ap <- strsplit(ap, "", fixed = TRUE)[[1L]]
    asub <- strsplit(asub, "", fixed = TRUE)[[1L]]
    subject_has_base <- asub != "-"
    pattern_has_base <- ap != "-"
    subject_positions <- as.integer(subject_start) + cumsum(subject_has_base) - 1L
    positions <- subject_positions[subject_has_base & pattern_has_base]
    unique(as.integer(positions[!is.na(positions) & positions >= 1L & positions <= n]))
  }

  for (read_counter in seq_along(read_indices)) {
    ri <- read_indices[read_counter]
    if (read_counter == 1L || read_counter %% progress_every == 0L ||
        read_counter == n_good_reads) {
      message("[TALOS] Ambiguity alignment progress: ", read_counter, "/",
              n_good_reads, " read ends")
    }
    read_seq <- read_seqs[ri]
    # Reverse-complement once per read end and reuse it for seed screening and
    # the bounded alternative search; the previous code computed it twice for
    # ambiguous reads without changing the candidate set.
    rc_read <- if (isTRUE(try_reverse) || has_best_positions) tryCatch(
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(read_seq))),
      error = function(e) NA_character_
    ) else NA_character_
    best_positions <- if (has_best_positions) {
      unique(as.integer(best_result$per_read_positions[[ri]]))
    } else integer(0)
    best_positions <- best_positions[!is.na(best_positions) &
      best_positions >= 1L & best_positions <= n]

    if (has_best_positions && !length(best_positions)) next
    if (has_best_positions && length(best_positions)) {
      oriented_for_seed <- read_seq
      rc_for_seed <- rc_read
      if (!is.na(rc_for_seed) && !identical(rc_for_seed, read_seq)) {
        oriented_for_seed <- c(oriented_for_seed, rc_for_seed)
      }
      # A read is fast-pathed only when at least two exact, template-unique
      # k-mer anchors agree with the offset/orientation of the primary local
      # alignment. Merely lacking a repeated exact k-mer is not enough to
      # claim unique placement.
      selected_orientation <- if (length(best_result$alignment_orientation) >= ri) {
        best_result$alignment_orientation[ri]
      } else NA_character_
      selected_seq <- if (identical(selected_orientation, "reverse") && !is.na(rc_for_seed)) {
        rc_for_seed
      } else read_seq
      seed_offsets <- collinear_unique_seed_offsets(selected_seq)
      direct_offset <- if (length(best_result$aligned_subject_start) >= ri &&
                           length(best_result$aligned_query_start) >= ri) {
        best_result$aligned_subject_start[ri] - best_result$aligned_query_start[ri]
      } else NA_integer_
      unique_collinear <- length(seed_offsets) > 0L && is.finite(direct_offset) &&
        as.character(direct_offset) %in% seed_offsets
      if (unique_collinear) {
        per_read_certain_positions[[ri]] <- best_positions
        per_read_possible_positions[[ri]] <- best_positions
        n_qualifying <- n_qualifying + 1L
        placement_counts <- c(placement_counts, 1L)
        next
      }
    }

    oriented <- read_seq
    if (isTRUE(try_reverse)) {
      rc <- rc_read
      if (!is.na(rc) && !identical(rc, read_seq)) oriented <- c(oriented, rc)
    }

    hits <- list()
    primary_valid <- FALSE
    primary_orientation <- NA_character_
    primary_start <- primary_end <- NA_integer_
    if (has_best_positions && length(best_positions)) {
      base_score <- suppressWarnings(as.numeric(best_result$score_by_read[ri]))
      base_pid <- suppressWarnings(as.numeric(best_result$pid_by_read[ri]))
      base_qwidth <- suppressWarnings(as.integer(best_result$aligned_query_width[ri]))
      primary_orientation <- if (length(best_result$alignment_orientation) >= ri) {
        as.character(best_result$alignment_orientation[ri])
      } else NA_character_
      primary_start <- if (length(best_result$aligned_subject_start) >= ri) {
        suppressWarnings(as.integer(best_result$aligned_subject_start[ri]))
      } else NA_integer_
      primary_end <- if (length(best_result$aligned_subject_end) >= ri) {
        suppressWarnings(as.integer(best_result$aligned_subject_end[ri]))
      } else NA_integer_
      primary_valid <- is.finite(base_score) && is.finite(base_pid) &&
        base_pid >= min_pid && !is.na(base_qwidth) && base_qwidth >= min_aligned &&
        primary_orientation %in% c("forward", "reverse") &&
        !is.na(primary_start) && !is.na(primary_end) &&
        primary_start >= 1L && primary_end <= n && primary_end >= primary_start
      if (primary_valid) {
        hits[[1L]] <- list(
          positions = best_positions, score = base_score, pid = base_pid,
          qwidth = base_qwidth, orientation = primary_orientation
        )
      }
    }
    read_reached_search_limit <- FALSE
    for (orientation_i in seq_along(oriented)) {
      orientation_label <- if (orientation_i == 1L) "forward" else "reverse"
      masked <- template
      # The primary local alignment is already available in best_result and
      # has been retained above. Mask its subject span before searching this
      # orientation so pairwiseAlignment does not spend another DP pass
      # rediscovering the same placement.
      if (primary_valid && identical(primary_orientation, orientation_label)) {
        chars <- strsplit(masked, "", fixed = TRUE)[[1L]]
        chars[primary_start:primary_end] <- "N"
        masked <- paste(chars, collapse = "")
      }
      pattern_dna <- Biostrings::DNAString(oriented[orientation_i])
      for (attempt in seq_len(max_alternatives)) {
        aln <- tryCatch(
          pwalign::pairwiseAlignment(
            pattern = pattern_dna,
            subject = Biostrings::DNAString(masked),
            type = "local"
          ),
          error = function(e) NULL
        )
        if (is.null(aln)) break
        if (attempt == max_alternatives) read_reached_search_limit <- TRUE
        pid <- suppressWarnings(as.numeric(pwalign::pid(aln))[1L])
        score <- suppressWarnings(as.numeric(BiocGenerics::score(aln))[1L])
        pr <- pwalign::pattern(aln)@range
        sr <- pwalign::subject(aln)@range
        qwidth <- tryCatch(as.integer(BiocGenerics::width(pr)[1L]), error = function(e) NA_integer_)
        sstart <- tryCatch(as.integer(BiocGenerics::start(sr)[1L]), error = function(e) NA_integer_)
        send <- tryCatch(as.integer(BiocGenerics::end(sr)[1L]), error = function(e) NA_integer_)
        if (is.na(sstart) || is.na(send) || send < sstart || sstart < 1L || send > n) break

        positions <- alignment_positions(aln, sstart)
        if (is.finite(pid) && pid >= min_pid && !is.na(qwidth) &&
            qwidth >= min_aligned && length(positions) >= min_aligned && is.finite(score)) {
          hits[[length(hits) + 1L]] <- list(
            positions = positions, score = score, pid = pid,
            qwidth = qwidth, orientation = if (orientation_i == 1L) "forward" else "reverse"
          )
        }

        # Mask the full aligned subject span to force discovery of a different
        # placement on the next iteration. N is a mismatch to ordinary A/C/G/T
        # read bases, so the masked locus cannot win again.
        chars <- strsplit(masked, "", fixed = TRUE)[[1L]]
        chars[sstart:send] <- "N"
        next_masked <- paste(chars, collapse = "")
        if (identical(next_masked, masked)) break
        masked <- next_masked
      }
    }
    if (!length(hits)) next
    if (read_reached_search_limit) n_search_limit <- n_search_limit + 1L

    # Collapse the same coordinate placement found in both orientations, then
    # retain only co-optimal placements. A small absolute score tolerance keeps
    # near-ties visible without allowing weak incidental local matches.
    sig <- vapply(hits, function(h) paste(h$positions, collapse = ","), character(1L))
    ord <- order(-vapply(hits, `[[`, numeric(1L), "score"),
                 -vapply(hits, `[[`, integer(1L), "qwidth"),
                 -vapply(hits, `[[`, numeric(1L), "pid"))
    hits <- hits[ord[!duplicated(sig[ord])]]
    top_score <- hits[[1L]]$score
    hits <- Filter(function(h) h$score >= top_score - score_delta, hits)
    if (!length(hits)) next

    placement_counts <- c(placement_counts, length(hits))
    n_qualifying <- n_qualifying + 1L
    hit_masks <- lapply(hits, function(h) {
      z <- logical(n)
      z[h$positions] <- TRUE
      z
    })
    certain <- Reduce(`&`, hit_masks)
    possible <- Reduce(`|`, hit_masks)
    if (length(hits) > 1L) n_multi <- n_multi + 1L
    if (any(possible & !certain)) n_ambiguous <- n_ambiguous + 1L
    certain_positions <- which(certain)
    possible_positions <- which(possible)
    per_read_certain_positions[[ri]] <- as.integer(certain_positions)
    per_read_possible_positions[[ri]] <- as.integer(possible_positions)
  }

  # Aggregate only after all read-end placements have been resolved. Each
  # read end contributes its unique positional set once; tabulate() performs
  # the depth accumulation in compiled base R while the lists preserve the
  # original read-level evidence for downstream junction unions.
  flatten_unique_positions <- function(x) {
    if (!length(x)) return(integer(0))
    unlist(lapply(x, function(z) {
      z <- unique(as.integer(z))
      z[!is.na(z) & z >= 1L & z <= n]
    }), use.names = FALSE)
  }
  certain_flat <- flatten_unique_positions(per_read_certain_positions)
  possible_flat <- flatten_unique_positions(per_read_possible_positions)
  certain_cov <- tabulate(certain_flat, nbins = n)
  possible_cov <- tabulate(possible_flat, nbins = n)
  ambiguous_cov <- pmax(0L, possible_cov - certain_cov)
  rr_c <- rle(certain_cov)
  rr_p <- rle(possible_cov)
  list(
    certain_cov = certain_cov,
    possible_cov = possible_cov,
    ambiguous_cov = ambiguous_cov,
    certain_rle = paste(paste0(rr_c$values, ":", rr_c$lengths), collapse = ","),
    possible_rle = paste(paste0(rr_p$values, ":", rr_p$lengths), collapse = ","),
    ambiguous_rle = .rle_coverage_string(ambiguous_cov),
    certain_percent = sum(certain_cov > 0L) / n * 100,
    possible_percent = sum(possible_cov > 0L) / n * 100,
    ambiguous_percent = sum(ambiguous_cov > 0L) / n * 100,
    ambiguous_read_ends = as.integer(n_ambiguous),
    multi_placement_read_ends = as.integer(n_multi),
    alternative_search_limit_read_ends = as.integer(n_search_limit),
    reads_with_qualifying_alignment = as.integer(n_qualifying),
    mean_placements = if (length(placement_counts)) mean(placement_counts) else 0,
    max_placements = if (length(placement_counts)) max(placement_counts) else 0L,
    per_read_certain_positions = per_read_certain_positions,
    per_read_possible_positions = per_read_possible_positions,
    status = if (n_qualifying) "computed_masked_subject_cooptimal_placements" else "computed_no_qualifying_alignments"
  )
}

# Selective tandem-junction coverage helpers. These are positional diagnostics;
# they do not assert copy identity or independently establish an ITD call.
# Find reads that genuinely span the novel tandem junction. Unlike a single
# short junction k-mer, the candidate is required to match a contiguous window
# containing substantial sequence from BOTH sides of the junction. pwalign
# performs batched local alignment in compiled code and supports small indels
# as well as substitutions.
.find_tandem_junction_candidates <- function(read_seqs, template,
    anchor_width = 12L, max_mismatch = 2L, chunk_size = 256L) {
  template <- toupper(as.character(template)[1L])
  read_seqs <- toupper(as.character(read_seqs))
  n <- nchar(template)
  anchor_width <- max(6L, min(as.integer(anchor_width)[1L], n))
  max_mismatch <- max(0L, as.integer(max_mismatch)[1L])
  chunk_size <- max(1L, as.integer(chunk_size)[1L])
  if (n < anchor_width || !length(read_seqs) ||
      !requireNamespace("Biostrings", quietly = TRUE) ||
      !requireNamespace("pwalign", quietly = TRUE)) return(integer(0))

  junction_window <- .template_junction_window(template, anchor_width)
  window_width <- nchar(junction_window)
  min_query_width <- max(1L, window_width - max_mismatch)
  min_pid <- 100 * (1 - max_mismatch / window_width)
  valid <- which(!is.na(read_seqs) & nzchar(read_seqs) &
    nchar(read_seqs) >= min_query_width)
  if (!length(valid)) return(integer(0))

  # Filter invalid alphabets per read, not per batch, so one malformed read
  # cannot abort screening for the remainder of the cohort.
  valid_dna <- vapply(read_seqs[valid], function(x) {
    !is.null(tryCatch(Biostrings::DNAString(x), error = function(e) NULL))
  }, logical(1L))
  valid <- valid[valid_dna]
  if (!length(valid)) return(integer(0))

  matched_forward <- logical(length(valid))
  matched_reverse <- logical(length(valid))
  starts <- seq.int(1L, length(valid), by = chunk_size)

  align_batch <- function(subject_set) {
    # pwalign supports parallel XStringSet inputs. A repeated query set paired
    # with the read set performs indel-aware local alignments in compiled code.
    pattern_set <- Biostrings::DNAStringSet(
      rep(junction_window, length(subject_set))
    )
    aln <- pwalign::pairwiseAlignment(
      pattern = pattern_set, subject = subject_set, type = "local"
    )
    pid <- suppressWarnings(as.numeric(pwalign::pid(aln)))
    pr <- pwalign::pattern(aln)
    qwidth <- suppressWarnings(as.integer(BiocGenerics::width(pr@range)))
    if (length(pid) != length(subject_set) || length(qwidth) != length(subject_set)) {
      stop("Vectorized junction alignment returned an unexpected number of results")
    }
    is.finite(pid) & pid + 1e-8 >= min_pid &
      !is.na(qwidth) & qwidth >= min_query_width
  }

  for (chunk_start in starts) {
    ii <- seq.int(chunk_start, min(length(valid), chunk_start + chunk_size - 1L))
    dna <- Biostrings::DNAStringSet(read_seqs[valid[ii]])
    reverse_dna <- Biostrings::reverseComplement(dna)
    matched_forward[ii] <- align_batch(dna)
    matched_reverse[ii] <- align_batch(reverse_dna)
  }

  orientation <- ifelse(matched_forward, "forward",
                        ifelse(matched_reverse, "reverse", NA_character_))
  matched <- which(!is.na(orientation))
  out <- valid[matched]
  attr(out, "orientation") <- orientation[matched]
  out
}

# Merge primary-template and junction-window evidence by read end. The
# junction window's first half maps to the ITD unit's terminal flank and its
# second half maps to the initial flank. Rebuild depth from per-read unions so
# distinct read ends add depth, while the same read end cannot be counted twice.
.merge_junction_read_coverage <- function(base, junction_result, junction_idx,
    template_length, anchor_width, n_read_ends = NULL) {
  n <- as.integer(template_length)[1L]
  f <- as.integer(anchor_width)[1L]
  nr <- if (is.null(n_read_ends)) length(base$per_read_certain_positions) else as.integer(n_read_ends)[1L]
  if (is.null(base$per_read_certain_positions) || length(base$per_read_certain_positions) != nr)
    base$per_read_certain_positions <- rep(list(integer(0)), nr)
  if (is.null(base$per_read_possible_positions) || length(base$per_read_possible_positions) != nr)
    base$per_read_possible_positions <- rep(list(integer(0)), nr)
  if (length(base$per_read_certain_positions) != length(base$per_read_possible_positions))
    stop("Per-read certain/possible coverage vectors have inconsistent lengths")

  project_positions <- function(pos) {
    pos <- unique(as.integer(pos))
    pos <- pos[!is.na(pos) & pos >= 1L & pos <= 2L * f]
    unique(as.integer(c(
      if (any(pos <= f)) n - f + pos[pos <= f],
      if (any(pos > f)) pos[pos > f] - f
    )))
  }
  for (j in seq_along(junction_idx)) {
    original_i <- junction_idx[j]
    j_possible <- project_positions(junction_result$per_read_possible_positions[[j]])
    j_certain <- project_positions(junction_result$per_read_certain_positions[[j]])
    if (length(intersect(j_possible, seq.int(n - f + 1L, n))) < 4L ||
        length(intersect(j_possible, seq_len(f))) < 4L) next
    base$per_read_possible_positions[[original_i]] <- unique(c(
      base$per_read_possible_positions[[original_i]], j_possible
    ))
    base$per_read_certain_positions[[original_i]] <- unique(c(
      base$per_read_certain_positions[[original_i]], j_certain
    ))
  }

  # Vectorized depth rebuild: deduplicate positions within each read end first,
  # then tabulate all read-end positions in compiled base code. This preserves
  # exact per-read union semantics while avoiding an R-level update loop.
  clean_positions <- function(x) {
    z <- unique(as.integer(x))
    z[!is.na(z) & z >= 1L & z <= n]
  }
  certain_flat <- unlist(lapply(base$per_read_certain_positions, clean_positions),
                         use.names = FALSE)
  possible_flat <- unlist(lapply(base$per_read_possible_positions, clean_positions),
                          use.names = FALSE)
  base$certain_cov <- tabulate(certain_flat, nbins = n)
  base$possible_cov <- tabulate(possible_flat, nbins = n)
  base$ambiguous_cov <- pmax(0L, base$possible_cov - base$certain_cov)
  base$certain_percent <- if (n) 100 * mean(base$certain_cov > 0L) else 0
  base$possible_percent <- if (n) 100 * mean(base$possible_cov > 0L) else 0
  base$ambiguous_percent <- if (n) 100 * mean(base$ambiguous_cov > 0L) else 0
  base$certain_rle <- .rle_coverage_string(base$certain_cov)
  base$possible_rle <- .rle_coverage_string(base$possible_cov)
  base$ambiguous_rle <- .rle_coverage_string(base$ambiguous_cov)
  base
}

# Evaluate the single-copy template once, then selectively analyse reads that
# span the actual copy-to-copy junction. The second alignment is restricted to
# a short junction window (default 24 bp), not the entire doubled template.
# Read-end-level masks are merged by union before depth is rebuilt, avoiding
# both pmax undercounting and double-counting evidence from the same read end.
.align_reads_to_tandem_expanded_template <- function(read_seqs, template,
    min_pid = 90, min_aligned = 20L, max_alternatives = 8L, score_delta = 1,
    best_result = NULL, seed_k = 9L, progress_every = 25L,
    junction_anchor_width = 12L, junction_max_mismatch = 2L) {
  template <- toupper(as.character(template)[1L])
  n <- nchar(template)
  read_seqs <- toupper(as.character(read_seqs))
  n_reads <- length(read_seqs)
  if (!n) return(.align_reads_to_itd_template_ambiguity(read_seqs, template))
  anchor_width <- max(6L, min(as.integer(junction_anchor_width)[1L], n))
  junction_window <- .template_junction_window(template, anchor_width)

  message("[TALOS] Ambiguity-aware ITD coverage: ", n_reads,
          " read ends; primary placements plus two-sided junction screening")
  base <- .align_reads_to_itd_template_ambiguity(
    read_seqs = read_seqs, template = template,
    min_pid = min_pid, min_aligned = min_aligned,
    max_alternatives = max_alternatives, score_delta = score_delta,
    best_result = best_result, seed_k = seed_k,
    progress_every = progress_every
  )

  junction_idx <- .find_tandem_junction_candidates(
    read_seqs, template, anchor_width = anchor_width,
    max_mismatch = junction_max_mismatch
  )
  message("[TALOS] Two-sided tandem-junction screen: ", length(junction_idx),
          " candidates of ", n_reads, " read ends (", 2L * anchor_width,
          " bp junction window; <=", junction_max_mismatch, " edits)")

  if (length(junction_idx)) {
    # This short synthetic subject contains only the two flanks. Restricting
    # alternative enumeration to two placements is sufficient for this local
    # junction diagnostic and avoids repeated DP over a 600 bp subject.
    junction_orientation <- attr(junction_idx, "orientation")
    junction_reads <- read_seqs[junction_idx]
    if (length(junction_orientation) == length(junction_reads)) {
      reverse_idx <- which(junction_orientation == "reverse")
      if (length(reverse_idx)) {
        junction_reads[reverse_idx] <- vapply(junction_reads[reverse_idx], function(x) {
          as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
        }, character(1L), USE.NAMES = FALSE)
      }
    }
    # Candidate screening has already selected the matching strand. Avoid
    # repeating reverse-complement alignment; retain two alternative local
    # placements so repeat ambiguity is still represented.
    junction_result <- .align_reads_to_itd_template_ambiguity(
      read_seqs = junction_reads, template = junction_window,
      min_pid = min_pid, min_aligned = min(min_aligned, nchar(junction_window)),
      max_alternatives = min(2L, max_alternatives), score_delta = score_delta,
      best_result = NULL, seed_k = min(seed_k, nchar(junction_window)),
      progress_every = progress_every, try_reverse = FALSE
    )

    project_positions <- function(pos) {
      pos <- unique(as.integer(pos))
      pos <- pos[!is.na(pos) & pos >= 1L & pos <= 2L * anchor_width]
      unique(as.integer(c(
        if (any(pos <= anchor_width)) n - anchor_width + pos[pos <= anchor_width],
        if (any(pos > anchor_width)) pos[pos > anchor_width] - anchor_width
      )))
    }

    base <- .merge_junction_read_coverage(
      base, junction_result, junction_idx,
      template_length = n, anchor_width = anchor_width,
      n_read_ends = n_reads
    )
    base$junction_candidate_read_ends <- as.integer(length(junction_idx))
    base$junction_aligned_read_ends <- as.integer(sum(vapply(
      junction_result$per_read_possible_positions, length, integer(1L)) > 0L))
    base$junction_search_status <- junction_result$status
    base$template_mode <- "single_copy_plus_selective_two_sided_junction_window"
  } else {
    base$junction_candidate_read_ends <- 0L
    base$junction_aligned_read_ends <- 0L
    base$junction_search_status <- "no_two_sided_junction_candidates"
    base$template_mode <- "single_copy_plus_two_sided_junction_screen"
  }
  message("[TALOS] Ambiguity-aware ITD coverage complete: ", n_reads,
          " read ends; ", length(junction_idx), " junction candidates; certain=",
          sprintf("%.1f%%", base$certain_percent), ", possible=",
          sprintf("%.1f%%", base$possible_percent))
  base
}

# Best local identity between an observed read-derived fragment and the
# expected reference-derived duplicated segment. This deliberately compares
# a partial fragment wherever it aligns, not only against the first bases of
# the reference segment (the latter silently mis-scored internal fragments).
.local_sequence_match <- function(query, subject) {
  empty <- list(
    pid = NA_real_, query_coverage = 0, aligned_width = 0L,
    orientation = NA_character_, subject_start = NA_integer_,
    subject_end = NA_integer_
  )
  if (length(query) != 1L || is.na(query) || !nzchar(query) ||
      length(subject) != 1L || is.na(subject) || !nzchar(subject)) return(empty)
  query <- toupper(as.character(query))
  subject <- toupper(as.character(subject))

  rc <- if (requireNamespace("Biostrings", quietly = TRUE)) {
    tryCatch(
      as.character(Biostrings::reverseComplement(Biostrings::DNAString(query))),
      error = function(e) NA_character_
    )
  } else {
    NA_character_
  }

  # Exact-match fallback remains useful for minimal installations, but it is
  # deliberately not allowed to invent an approximate identity score.
  if (!requireNamespace("pwalign", quietly = TRUE) ||
      !requireNamespace("Biostrings", quietly = TRUE)) {
    hit <- regexpr(query, subject, fixed = TRUE)[1L]
    if (hit > 0L) {
      return(list(pid = 100, query_coverage = 1, aligned_width = nchar(query),
                  orientation = "tandem", subject_start = as.integer(hit),
                  subject_end = as.integer(hit + nchar(query) - 1L)))
    }
    if (!is.na(rc) && !identical(rc, query)) {
      hit_rc <- regexpr(rc, subject, fixed = TRUE)[1L]
      if (hit_rc > 0L) {
        return(list(pid = 100, query_coverage = 1, aligned_width = nchar(query),
                    orientation = "inverted", subject_start = as.integer(hit_rc),
                    subject_end = as.integer(hit_rc + nchar(query) - 1L)))
      }
    }
    return(empty)
  }

  queries <- query
  orientations <- "tandem"
  if (!is.na(rc) && !identical(rc, query)) {
    queries <- c(queries, rc)
    orientations <- c(orientations, "inverted")
  }
  candidates <- lapply(seq_along(queries), function(i) {
    aln <- tryCatch(
      pwalign::pairwiseAlignment(
        pattern = Biostrings::DNAString(queries[i]),
        subject = Biostrings::DNAString(subject),
        type = "local"
      ),
      error = function(e) NULL
    )
    if (is.null(aln)) return(NULL)
    pid <- suppressWarnings(as.numeric(pwalign::pid(aln))[1L])
    query_range <- pwalign::pattern(aln)@range
    subject_range <- pwalign::subject(aln)@range
    aligned_width <- tryCatch(
      as.integer(BiocGenerics::width(query_range)[1L]),
      error = function(e) NA_integer_
    )
    subject_start <- tryCatch(
      as.integer(BiocGenerics::start(subject_range)[1L]),
      error = function(e) NA_integer_
    )
    subject_end <- tryCatch(
      as.integer(BiocGenerics::end(subject_range)[1L]),
      error = function(e) NA_integer_
    )
    if (!is.finite(pid) || is.na(aligned_width) || aligned_width <= 0L ||
        is.na(subject_start) || is.na(subject_end) || subject_end < subject_start) {
      return(NULL)
    }
    list(
      pid = pid,
      query_coverage = min(1, aligned_width / nchar(queries[i])),
      aligned_width = aligned_width,
      orientation = orientations[i],
      subject_start = subject_start,
      subject_end = subject_end
    )
  })
  candidates <- Filter(Negate(is.null), candidates)
  if (!length(candidates)) return(empty)
  ord <- order(
    -vapply(candidates, `[[`, integer(1L), "aligned_width"),
    -vapply(candidates, `[[`, numeric(1L), "query_coverage"),
    -vapply(candidates, `[[`, numeric(1L), "pid"),
    seq_along(candidates)
  )
  candidates[[ord[1L]]]
}

# Summarise uncovered positions in a zero-based coverage vector using
# 1-based inclusive coordinates for human-readable / TSV output.
.coverage_uncovered_summary <- function(cov_vec) {
  if (length(cov_vec) == 0L) {
    return(list(bases = NA_integer_, intervals = NA_character_))
  }
  missing <- which(is.na(cov_vec) | cov_vec <= 0L)
  if (!length(missing)) {
    return(list(bases = 0L, intervals = ""))
  }
  group <- cumsum(c(TRUE, diff(missing) > 1L))
  starts <- vapply(split(missing, group), `[`, integer(1L), 1L)
  ends <- vapply(split(missing, group), function(x) x[length(x)], integer(1L))
  labels <- ifelse(starts == ends, as.character(starts), paste0(starts, "-", ends))
  list(bases = as.integer(length(missing)), intervals = paste(labels, collapse = "|"))
}

# Per-source depth summary. Empty vectors mean unavailable; a correctly sized
# all-zero vector is a measured zero-depth result.
.coverage_depth_summary <- function(cov_vec) {
  if (is.null(cov_vec) || !length(cov_vec)) {
    return(list(covered_bases = NA_integer_, total_depth = NA_real_,
                mean_depth = NA_real_, max_depth = NA_integer_))
  }
  x <- suppressWarnings(as.numeric(cov_vec))
  x[!is.finite(x) | x < 0] <- 0
  list(
    covered_bases = as.integer(sum(x > 0)),
    total_depth = sum(x),
    mean_depth = mean(x),
    max_depth = as.integer(max(x))
  )
}

.rle_coverage_string <- function(cov_vec) {
  if (is.null(cov_vec) || !length(cov_vec)) return(NA_character_)
  rr <- rle(as.integer(cov_vec))
  paste(paste0(rr$values, ":", rr$lengths), collapse = ",")
}


# Select a full-length coordinate template for ITD coverage without changing
# the clinical sequence call. When the read-derived sequence is unresolved,
# the reference segment is used ONLY to define ITD-relative coordinates for
# mapping independent read / PE / k-mer evidence; ITD_Sequence and
# SequenceSource remain unresolved.
.resolve_itd_coverage_template <- function(
    itd_seq,
    best_len,
    ref_dna,
    local_seg_start,
    is_ptd_scale = FALSE,
    seg_in_ref = FALSE,
    sequence_imputed = FALSE
) {
  has_sequence <- !is.null(itd_seq) && length(itd_seq) == 1L &&
    !is.na(itd_seq) && nzchar(as.character(itd_seq))
  if (has_sequence) {
    source <- if (isTRUE(sequence_imputed)) {
      "reference_imputed_sequence"
    } else {
      "read_derived_sequence"
    }
    return(list(sequence = as.character(itd_seq), source = source))
  }

  len <- suppressWarnings(as.integer(best_len)[1L])
  start <- suppressWarnings(as.integer(local_seg_start)[1L])
  ref_len <- if (is.null(ref_dna) || length(ref_dna) == 0L ||
                is.na(ref_dna[1L])) 0L else nchar(as.character(ref_dna[1L]))

  if (!isTRUE(is_ptd_scale) && isTRUE(seg_in_ref) &&
      length(len) == 1L && !is.na(len) && len > 0L &&
      length(start) == 1L && !is.na(start) && start >= 1L &&
      ref_len >= start + len - 1L) {
    template <- substr(as.character(ref_dna[1L]), start, start + len - 1L)
    if (nchar(template) == len) {
      return(list(
        sequence = template,
        source = "reference_segment_for_coverage_only"
      ))
    }
  }

  list(sequence = NA_character_, source = "unavailable")
}

.compute_variant_metrics <- function(
  cluster_bps,
  best_len,
  support_rows,
  genomic_start,
  ref_dna,
  gene_config,
  all_reads_cov,
  all_reads_for_coverage,
  all_pairs,
  wildtype_info,
  ptd_mode,
  min_support,
  min_wt_reads,
  nominal_read_len,
  max_correction,
  min_mapq,
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
  anchor_confirmed_qnames = NA_integer_,
  # v0.4.20: adapter screening + base-quality context, threaded from detect_itd().
  adapter_positive_qnames = NULL,
  adapter_screened_reads = NA_integer_,
  adapter_positive_reads = NA_integer_,
  adapter_n_configured = 0L,
  adapter_min_overlap = NA_integer_,
  base_quality_loaded = NA,
  min_base_quality = 0,
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
  support_consistency_status <- if (isTRUE(do_consistency)) "not_evaluated" else "disabled"
  ref_match_total_status <- if (isTRUE(do_alignment_score)) "not_evaluated" else "disabled"
  alignment_score_status <- if (isTRUE(do_alignment_score)) "not_evaluated" else "disabled"
  orientation_status <- if (isTRUE(do_detect_orientation)) "not_evaluated" else "disabled"
  sequence_dependent_metrics_status <- "not_evaluated"
  itd_coverage_percent <- NA_real_
  itd_coverage_rle <- NA_character_
  itd_coverage_rle_pe <- NA_character_
  itd_coverage_percent_combined <- NA_real_
  itd_pe_spanning_pairs <- NA_integer_
  itd_coverage_rle_kmer <- NA_character_
  itd_coverage_rle_kmer_by_k <- list()
  itd_kmer_percent_by_k <- list()
  itd_kmer_reads_by_k <- list()
  itd_kmer_fragments_by_k <- list()
  itd_kmer_reads <- NA_integer_
  itd_coverage_percent_kmer <- NA_real_
  itd_coverage_direct_certain_pct <- NA_real_
  itd_coverage_direct_possible_pct <- NA_real_
  itd_coverage_direct_ambiguous_pct <- NA_real_
  itd_coverage_direct_ambiguous_read_ends <- NA_integer_
  itd_coverage_direct_multi_placement_read_ends <- NA_integer_
  itd_coverage_direct_search_limit_read_ends <- NA_integer_
  itd_coverage_direct_mean_placements <- NA_real_
  itd_coverage_direct_max_placements <- NA_integer_
  itd_coverage_rle_direct_certain <- NA_character_
  itd_coverage_rle_direct_possible <- NA_character_
  itd_coverage_rle_direct_ambiguous <- NA_character_
  itd_coverage_combined_certain_pct <- NA_real_
  itd_coverage_combined_possible_pct <- NA_real_
  itd_coverage_rle_combined_certain <- NA_character_
  itd_coverage_rle_combined_possible <- NA_character_
  itd_kmer_k_values <- "7,9,11,13,15"
  itd_coverage_reads_processed <- NA_integer_
  itd_coverage_reads_matched <- NA_integer_
  itd_coverage_mate_read_ends_added <- 0L
  itd_coverage_fragments_processed <- NA_integer_
  itd_coverage_fragments_matched <- NA_integer_
  itd_coverage_direct_alignment_rate_pct <- NA_real_
  itd_coverage_direct_median_pid <- NA_real_
  itd_coverage_direct_median_aligned_bases <- NA_real_
  itd_coverage_direct_mean_pid <- NA_real_
  itd_coverage_direct_median_query_coverage_pct <- NA_real_
  itd_coverage_direct_n_any_alignment <- NA_integer_
  itd_coverage_direct_n_qualifying_alignment <- NA_integer_
  itd_coverage_percent_pe <- NA_real_
  itd_coverage_direct_status <- if (isTRUE(do_itd_coverage)) "not_evaluated" else "disabled"
  itd_coverage_pe_status <- if (isTRUE(do_itd_coverage)) "not_evaluated" else "disabled"
  itd_coverage_kmer_status <- if (isTRUE(do_itd_coverage)) "not_evaluated" else "disabled"
  itd_coverage_status <- if (isTRUE(do_itd_coverage)) "not_evaluated" else "disabled"
  itd_pe_observation_pairs <- NA_integer_
  itd_pe_valid_pairs <- NA_integer_
  itd_pe_support_pairs <- NA_integer_
  itd_pe_qname_usable_pairs <- NA_integer_
  itd_pe_qname_matched_pairs <- NA_integer_
  itd_pe_qname_status <- if (isTRUE(do_itd_coverage)) "not_evaluated" else "disabled"
  itd_kmer_tested_reads <- NA_integer_
  itd_kmer_matched_rate_pct <- NA_real_
  itd_kmer_tested_fragments <- NA_integer_
  itd_kmer_matched_fragments <- NA_integer_
  itd_kmer_matched_fragment_rate_pct <- NA_real_
  itd_uncovered_bases_pe <- NA_integer_
  itd_uncovered_intervals_pe <- NA_character_
  itd_uncovered_bases_kmer <- NA_integer_
  itd_uncovered_intervals_kmer <- NA_character_
  length_pe_nbackground <- NA_integer_
  length_pe_status <- if (is.null(all_pairs)) "no_pair_data" else "not_evaluated"
  coverage_template_length <- NA_integer_
  cov_kmer <- integer(0)
  cov_kmer_by_k <- list()
  cov_direct <- integer(0)
  cov_direct_certain <- integer(0)
  cov_direct_possible <- integer(0)
  cov_pe <- integer(0)
  cov_display <- integer(0)
  cov_display_certain <- integer(0)
  cov_display_possible <- integer(0)
  # [FIX v9] Declares whether cov_direct is the no-pwalign fallback's 0/1
  # covered-flags (not read counts); consumed by .fold_display_coverage()
  # so the fallback layer is ordered least specific WITHOUT value-based
  # inference, which misclassified low-depth count vectors.
  direct_binary_fallback <- FALSE
  # [FIX v7] pe_cov must exist even when the PE block below does not run
  # (no all_pairs / .pair_info() returned NULL): the k-mer block references
  # it when recomputing the combined percentage, which used to raise
  # "object 'pe_cov' not found" on that path.
  pe_cov <- NULL
  # Strict direct-alignment RLE; the combined display vector is direct
  # read-sequence alignment plus sequence-anchored k-mer evidence. PE
  # coordinate projection remains diagnostic and cannot fill allele gaps.
  itd_coverage_rle_direct <- NA_character_
  itd_uncovered_bases_direct <- NA_integer_
  itd_uncovered_bases_combined <- NA_integer_
  itd_uncovered_intervals_combined <- NA_character_
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
  length_pe_nbackground <- NA_integer_
  length_pe_nvalid <- NA_integer_
  length_pe_median_spanning <- NA_real_
  length_pe_median_background <- NA_real_
  length_pe_status <- if (is.null(all_pairs)) "no_pair_data" else "not_evaluated"
  length_pe_delta <- NA_real_
  length_pe_concordant <- NA
  length_pe_support <- NA_real_
  pe_est_support <- list()
  length_pe_support_nspanning <- NA_integer_
  length_pe_support_nspanning_all <- NA_integer_
  length_pe_support_nbackground <- NA_integer_
  length_pe_support_nvalid <- NA_integer_
  length_pe_support_nqname_matched <- NA_integer_
  length_pe_support_median_spanning <- NA_real_
  length_pe_support_median_background <- NA_real_
  length_pe_support_status <- if (is.null(all_pairs)) "no_pair_data" else "not_evaluated"
  length_pe_support_delta <- NA_real_
  length_pe_support_concordant <- NA
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
  # Keep a separate read-end-level view for direct ITD coverage. The clinical
  # support count, MAPQ, strand bias and soft-clip summaries remain one row per
  # fragment/QNAME, while base-level direct coverage includes both mates when
  # both provide independent aligned read sequences.
  coverage_support_rows <- .deduplicate_coverage_read_rows(support_rows)
  support_rows <- .deduplicate_support_rows(support_rows, len_specific_bp)
  if (is.null(support_rows) || nrow(support_rows) == 0L) return(NULL)
  support_qnames <- unique(support_rows$read_name)
  # Add the actual BAM sequences for the paired mates of support QNAMEs.
  # This makes mate contribution direct read alignment evidence; genomic
  # coordinate projection alone is not treated as ITD-allele coverage.
  mate_coverage_rows <- .extract_support_mate_coverage_rows(
    all_reads_for_coverage, support_qnames
  )
  itd_coverage_mate_read_ends_added <- as.integer(nrow(mate_coverage_rows))
  if (nrow(mate_coverage_rows)) {
    all_cols <- union(names(coverage_support_rows), names(mate_coverage_rows))
    for (nm in setdiff(all_cols, names(coverage_support_rows))) {
      coverage_support_rows[[nm]] <- NA
    }
    for (nm in setdiff(all_cols, names(mate_coverage_rows))) {
      mate_coverage_rows[[nm]] <- NA
    }
    coverage_support_rows <- .deduplicate_coverage_read_rows(rbind(
      coverage_support_rows[, all_cols, drop = FALSE],
      mate_coverage_rows[, all_cols, drop = FALSE]
    ))
  }
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
  # Terminal-ITD audit: summarize supporting alignments overlapping the final
  # 30 bp of the candidate interval. Since v0.4.20 the BAM loader retains the
  # QUAL column and candidate rows carry each read's mean Phred quality
  # (read_mean_bq), so base-quality metrics report real measurements; the
  # explicit "unavailable" statuses remain for frames built without QUAL.
  terminal_window_bp <- 30L
  terminal_start_genomic <- if (!is.na(best_len) && best_len > 0L) {
    len_specific_bp + best_len - terminal_window_bp
  } else NA_real_
  terminal_end_genomic <- if (!is.na(best_len) && best_len > 0L) {
    len_specific_bp + best_len - 1L
  } else NA_real_
  terminal_idx <- rep(FALSE, nrow(support_rows))
  if (is.finite(terminal_start_genomic) &&
      all(c("read_start", "read_end") %in% names(support_rows))) {
    rs <- suppressWarnings(as.numeric(support_rows$read_start))
    re <- suppressWarnings(as.numeric(support_rows$read_end))
    terminal_idx <- is.finite(rs) & is.finite(re) &
      rs <= terminal_end_genomic & re >= terminal_start_genomic
  }
  itd_terminal_supporting_alignments <- as.integer(sum(terminal_idx, na.rm = TRUE))
  itd_terminal_unique_qnames <- if (any(terminal_idx) && "read_name" %in% names(support_rows)) {
    as.integer(length(unique(stats::na.omit(.normalize_read_qname(support_rows$read_name[terminal_idx])))))
  } else 0L
  itd_terminal_reverse_alignments <- if (any(terminal_idx) && "is_reverse" %in% names(support_rows)) {
    as.integer(sum(as.logical(support_rows$is_reverse[terminal_idx]), na.rm = TRUE))
  } else NA_integer_
  itd_terminal_cigar_values <- if (any(terminal_idx)) as.character(support_cigars[terminal_idx]) else character(0)
  itd_terminal_cigar_softclip_alignments <- as.integer(sum(vapply(
    itd_terminal_cigar_values, function(z) {
      if (is.na(z)) return(FALSE)
      sc <- .get_softclip_lengths(z)
      isTRUE(sc$lead > 0L) || isTRUE(sc$trail > 0L)
    }, logical(1L)
  )))
  # Base-quality and adapter screening summaries (v0.4.20). Candidate rows
  # carry each read's mean Phred quality in `read_mean_bq`; the adapter flag
  # set comes from the region-wide .screen_adapter_reads() pass. Both are
  # computed, not hard-coded: when the QUAL column is absent the status
  # strings keep reporting unavailability explicitly.
  bq_vals <- if ("read_mean_bq" %in% names(support_rows)) {
    suppressWarnings(as.numeric(support_rows$read_mean_bq))
  } else numeric(0)
  bq_vals <- bq_vals[is.finite(bq_vals)]
  itd_mean_base_quality <- if (length(bq_vals)) round(mean(bq_vals), 1) else NA_real_
  itd_min_base_quality <- if (length(bq_vals)) round(min(bq_vals), 1) else NA_real_
  terminal_bq_vals <- if (any(terminal_idx) && "read_mean_bq" %in% names(support_rows)) {
    tv <- suppressWarnings(as.numeric(support_rows$read_mean_bq[terminal_idx]))
    tv[is.finite(tv)]
  } else numeric(0)
  itd_terminal_mean_base_quality <- if (length(terminal_bq_vals)) {
    round(mean(terminal_bq_vals), 1)
  } else NA_real_
  itd_terminal_base_quality_status <- if (!is.na(itd_terminal_mean_base_quality)) {
    "computed_over_terminal_window_support_rows"
  } else if (!isTRUE(base_quality_loaded)) {
    "unavailable_candidate_rows_do_not_retain_BAM_QUAL"
  } else {
    "unavailable_no_terminal_window_rows_with_quality"
  }
  itd_base_quality_status <- if (isTRUE(base_quality_loaded)) {
    if (isTRUE(min_base_quality > 0)) {
      paste0("applied_min_mean_base_quality_", min_base_quality, "_per_read")
    } else {
      "measured_not_filtered_min_base_quality_0"
    }
  } else {
    "unavailable_BAM_QUAL_not_loaded_into_candidate_rows"
  }
  adapter_support_hits <- if (length(adapter_positive_qnames)) {
    as.integer(sum(support_qnames %in% adapter_positive_qnames))
  } else 0L
  itd_adapter_screening_status <- if (adapter_n_configured > 0L) {
    paste0(
      "performed_", adapter_n_configured, "_adapters_min_overlap_",
      adapter_min_overlap, "_exact_prefix_match"
    )
  } else {
    "disabled_no_adapter_reference_configured"
  }
  terminal_evidence_rows <- which(terminal_idx)
  terminal_evidence_truncated <- length(terminal_evidence_rows) > 200L
  terminal_evidence_rows <- head(terminal_evidence_rows, 200L)
  terminal_evidence_field <- if (length(terminal_evidence_rows)) {
    paste(vapply(terminal_evidence_rows, function(i) {
      get_value <- function(nm, default = "NA") {
        if (!nm %in% names(support_rows)) return(default)
        value <- as.character(support_rows[[nm]][i])
        if (!length(value) || is.na(value) || !nzchar(value)) default else gsub("[|;\\t\\r\\n]", "_", value)
      }
      paste(c(
        get_value("read_name"), get_value("read_start"), get_value("read_end"),
        get_value("cigar"), get_value("is_reverse"), get_value("mapq"),
        get_value("read_seq")
      ), collapse = "|")
    }, character(1L)), collapse = ";")
  } else NA_character_

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
  sequence_unresolved <- FALSE
  observed_seq <- NA_character_
  observed_len <- 0L
  itd_observed_fragment_pid <- NA_real_
  itd_observed_fragment_query_coverage_pct <- NA_real_
  itd_observed_fragment_aligned_bases <- NA_integer_
  itd_observed_fragment_orientation <- NA_character_
  itd_observed_fragment_reference_start <- NA_integer_
  itd_observed_fragment_reference_end <- NA_integer_
  itd_observed_fragment_support_reads <- NA_integer_
  itd_observed_fragment_unique_sequences <- NA_integer_
  itd_observed_fragment_status <- "not_evaluated"

  valid_seqs <- support_rows$itd_seq[
    !is.na(support_rows$itd_seq) &
      nchar(support_rows$itd_seq) > 0L
  ]
  if (length(valid_seqs) > 0L) {
    seq_tab <- sort(table(valid_seqs), decreasing = TRUE)
    itd_seq <- names(seq_tab)[1L]
    itd_observed_fragment_support_reads <- as.integer(seq_tab[1L])
    itd_observed_fragment_unique_sequences <- as.integer(length(seq_tab))
    imputed <- FALSE
    observed_seq <- itd_seq
    observed_len <- nchar(itd_seq)
  }

  if (is.na(itd_seq)) {
    extracted_seqs <- character(0)
    extracted_qnames <- character(0)
    row_qnames <- .normalize_read_qname(as.character(support_rows$read_name))
    for (i in seq_len(nrow(support_rows))) {
      sc <- .get_softclips(support_rows$cigar[i], support_rows$read_seq[i])
      fragments <- c(sc$lead, sc$trail)
      fragments <- fragments[!is.na(fragments) & nzchar(fragments)]
      if (length(fragments)) {
        extracted_seqs <- c(extracted_seqs, fragments)
        extracted_qnames <- c(extracted_qnames, rep(row_qnames[i], length(fragments)))
      }
    }
    if (length(extracted_seqs) > 0L) {
      seq_tab <- sort(table(extracted_seqs), decreasing = TRUE)
      itd_seq <- names(seq_tab)[1L]
      selected_qnames <- unique(extracted_qnames[extracted_seqs == itd_seq])
      itd_observed_fragment_support_reads <- as.integer(length(selected_qnames))
      itd_observed_fragment_unique_sequences <- as.integer(length(seq_tab))
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
      !is_ptd_scale
  ) {
    # The observed soft-clip may originate inside the duplicated segment,
    # not necessarily at its first base. Appending the reference after the
    # clip therefore creates a coordinate-shifted chimera for internal clips.
    # Reconstruct the full duplication from the genomic reference only after
    # verifying that the observed fragment aligns well to SOMEWHERE within it.
    sequence_partial <- TRUE
    imputed <- TRUE
    partial_obs_len <- nchar(itd_seq)
    if (seg_in_ref) {
      ref_local_end <- local_seg_start + best_len - 1L
      ref_segment <- substr(ref_dna, local_seg_start, ref_local_end)
      observed_match <- .local_sequence_match(itd_seq, ref_segment)
      if (nchar(ref_segment) == best_len &&
          partial_obs_len >= 20L &&
          !is.na(observed_match$pid) &&
          observed_match$pid >= 85 &&
          observed_match$query_coverage >= 0.80 &&
          observed_match$aligned_width >= 20L) {
        itd_seq <- ref_segment
        observed_len <- partial_obs_len
      } else {
        # Keep the event/length call, but do not calculate a misleading
        # coverage percentage against a short fragment as if it were the
        # complete duplicated sequence.
        itd_seq <- NA_character_
        sequence_unresolved <- TRUE
        imputed <- FALSE
        observed_len <- partial_obs_len
      }
    } else {
      itd_seq <- NA_character_
      sequence_unresolved <- TRUE
      imputed <- FALSE
      observed_len <- partial_obs_len
    }
  }

  if (is.na(itd_seq) && observed_len == 0L && !is_ptd_scale && seg_in_ref) {
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

  # Observed-fragment identity/orientation is independent of whether a
  # full-length clinical ITD sequence could be validated. Compute it from the
  # actual observed fragment against the expected reference segment, and retain
  # its aligned coordinates and query coverage so a short local match cannot be
  # mistaken for a full-sequence validation.
  if (!is.na(observed_seq) && nzchar(observed_seq) && observed_len > 0L &&
      !is_ptd_scale && seg_in_ref) {
    observed_ref_segment <- substr(
      ref_dna, local_seg_start, local_seg_start + best_len - 1L
    )
    observed_match <- .local_sequence_match(observed_seq, observed_ref_segment)
    if (is.finite(observed_match$pid)) {
      ref_match_observed <- round(observed_match$pid, 1)
      itd_observed_fragment_pid <- as.numeric(observed_match$pid)
      itd_observed_fragment_query_coverage_pct <- 100 * as.numeric(observed_match$query_coverage)
      itd_observed_fragment_aligned_bases <- as.integer(observed_match$aligned_width)
      itd_observed_fragment_orientation <- as.character(observed_match$orientation)
      itd_observed_fragment_reference_start <- as.integer(observed_match$subject_start)
      itd_observed_fragment_reference_end <- as.integer(observed_match$subject_end)
      itd_observed_fragment_status <- "computed_local_fragment_alignment"
    } else {
      itd_observed_fragment_status <- "no_qualifying_local_fragment_alignment"
    }
  } else if (is_ptd_scale) {
    itd_observed_fragment_status <- "not_applicable_ptd_scale"
  } else if (!seg_in_ref) {
    itd_observed_fragment_status <- "unavailable_reference_segment_out_of_bounds"
  } else {
    itd_observed_fragment_status <- "unavailable_no_observed_fragment"
  }

  coverage_template_info <- .resolve_itd_coverage_template(
    itd_seq = itd_seq,
    best_len = best_len,
    ref_dna = ref_dna,
    local_seg_start = local_seg_start,
    is_ptd_scale = is_ptd_scale,
    seg_in_ref = seg_in_ref,
    sequence_imputed = imputed
  )
  coverage_template <- coverage_template_info$sequence
  coverage_template_source <- coverage_template_info$source

  if (!is.na(coverage_template) && nzchar(coverage_template)) {
    coverage_template_length <- as.integer(nchar(coverage_template))

    coverage_rows <- coverage_support_rows
    if (is.null(coverage_rows) || nrow(coverage_rows) == 0L) {
      coverage_rows <- support_rows
    }
    valid_mask <- !is.na(coverage_rows$read_seq) &
      nzchar(as.character(coverage_rows$read_seq))
    v_seqs <- as.character(coverage_rows$read_seq[valid_mask])
    v_qnames <- .normalize_read_qname(as.character(coverage_rows$read_name[valid_mask]))

    if (length(v_seqs) > 0L && (do_consistency || do_itd_coverage)) {
      # max_pairwise_alignments is a memory/chunk-size control, NOT a
      # total-read sampling cap. Every distinct supporting read-end alignment
      # is included; duplicate candidate rows are collapsed by read signature.
      batch_size <- suppressWarnings(as.integer(max_pairwise_alignments)[1L])
      if (is.na(batch_size) || batch_size < 1L) batch_size <- 100L
      batch_size <- min(500L, max(1L, batch_size))

      direct_result <- .align_reads_to_itd_template(
        read_seqs = v_seqs,
        template = coverage_template,
        chunk_size = batch_size,
        min_pid = 90,
        min_aligned = 20L
      )
      # Retain the legacy best-placement direct metric while calculating
      # repeat-aware lower/upper bounds across co-optimal placements.
      direct_ambiguity <- if (do_itd_coverage) {
        .align_reads_to_tandem_expanded_template(
          read_seqs = v_seqs, template = coverage_template,
          min_pid = 90, min_aligned = 20L,
          max_alternatives = 8L, score_delta = 1,
          best_result = direct_result, progress_every = 25L
        )
      } else NULL

      if (do_itd_coverage) {
        cov_direct <- direct_result$cov_vec
        cov_direct_certain <- direct_ambiguity$certain_cov
        cov_direct_possible <- direct_ambiguity$possible_cov
        itd_coverage_direct_certain_pct <- direct_ambiguity$certain_percent
        itd_coverage_direct_possible_pct <- direct_ambiguity$possible_percent
        itd_coverage_direct_ambiguous_pct <- direct_ambiguity$ambiguous_percent
        itd_coverage_direct_ambiguous_read_ends <- direct_ambiguity$ambiguous_read_ends
        itd_coverage_direct_multi_placement_read_ends <- direct_ambiguity$multi_placement_read_ends
        itd_coverage_direct_search_limit_read_ends <- direct_ambiguity$alternative_search_limit_read_ends
        itd_coverage_direct_mean_placements <- direct_ambiguity$mean_placements
        itd_coverage_direct_max_placements <- direct_ambiguity$max_placements
        itd_coverage_rle_direct_certain <- .rle_coverage_string(cov_direct_certain)
        itd_coverage_rle_direct_possible <- .rle_coverage_string(cov_direct_possible)
        itd_coverage_rle_direct_ambiguous <- .rle_coverage_string(direct_ambiguity$ambiguous_cov)
        itd_coverage_percent <- direct_result$percent
        itd_coverage_rle <- direct_result$rle_str
        itd_coverage_reads_processed <- direct_result$reads_processed
        itd_coverage_reads_matched <- direct_result$reads_matched
        itd_coverage_fragments_processed <- length(unique(stats::na.omit(v_qnames)))
        matched_read_idx <- is.finite(direct_result$pid_by_read) &
          direct_result$pid_by_read >= 90 &
          direct_result$aligned_query_width >= 20L
        itd_coverage_fragments_matched <- length(unique(stats::na.omit(v_qnames[matched_read_idx])))
        itd_coverage_direct_alignment_rate_pct <- if (direct_result$reads_processed > 0L) {
          100 * direct_result$reads_matched / direct_result$reads_processed
        } else NA_real_
        valid_pid <- direct_result$pid_by_read[
          is.finite(direct_result$pid_by_read) &
            direct_result$pid_by_read >= 90 &
            direct_result$aligned_query_width >= 20L
        ]
        valid_width <- direct_result$aligned_query_width[
          is.finite(direct_result$pid_by_read) &
            direct_result$pid_by_read >= 90 &
            direct_result$aligned_query_width >= 20L
        ]
        itd_coverage_direct_median_pid <- if (length(valid_pid)) stats::median(valid_pid) else NA_real_
        itd_coverage_direct_median_aligned_bases <- if (length(valid_width)) stats::median(valid_width) else NA_real_
        any_pid <- direct_result$pid_by_read[is.finite(direct_result$pid_by_read)]
        itd_coverage_direct_mean_pid <- if (length(any_pid)) mean(any_pid) else NA_real_
        valid_qcov <- direct_result$query_coverage_pct[
          is.finite(direct_result$pid_by_read) & direct_result$pid_by_read >= 90 &
            direct_result$aligned_query_width >= 20L
        ]
        itd_coverage_direct_median_query_coverage_pct <- if (length(valid_qcov)) stats::median(valid_qcov) else NA_real_
        itd_coverage_direct_n_any_alignment <- as.integer(sum(is.finite(direct_result$pid_by_read)))
        itd_coverage_direct_n_qualifying_alignment <- as.integer(direct_result$reads_matched)
        itd_coverage_direct_status <- if (direct_result$reads_matched > 0L) {
          "computed"
        } else {
          "computed_no_qualifying_direct_alignments"
        }
      }

      # Consistency is informative rather than a clinical filter. It is only
      # calculated against a fully read-derived sequence. A reference-derived
      # positional template can support coverage accounting, but cannot
      # independently validate the sequence from which it was copied.
      if (do_consistency) {
        if (identical(coverage_template_source, "read_derived_sequence")) {
          valid_pid <- is.finite(direct_result$pid_by_read)
          consistency_score <- if (any(valid_pid)) {
            sum(direct_result$pid_by_read[valid_pid] >= 90) /
              length(v_seqs) * 100
          } else {
            0
          }
          support_consistency_status <- if (direct_result$reads_matched > 0L) {
            "computed_read_derived_template_informational"
          } else {
            "computed_no_qualifying_alignments"
          }
        } else if (identical(coverage_template_source, "reference_imputed_sequence")) {
          # A reference-imputed sequence cannot independently validate read
          # consistency; keep the metric unavailable rather than scoring reads
          # against a template that was itself copied from the reference.
          support_consistency_status <- "not_comparable_reference_imputed_template"
        } else {
          # For unresolved sequences the reference template exists only to
          # coordinate coverage. ITDCoverageDirectAlignmentRatePct reports the
          # actual read-to-reference alignment rate separately.
          support_consistency_status <- "not_comparable_reference_only_coverage_template"
        }
      }
    } else if (do_itd_coverage || do_consistency) {
      if (do_consistency) {
        support_consistency_status <- "unavailable_no_support_read_sequences"
      }
      if (do_itd_coverage) {
        # No support sequence means direct read coverage cannot be measured.
        # Retain NA for the percentage/RLE rather than misreporting observed 0.
        cov_direct <- integer(0)
        itd_coverage_percent <- NA_real_
        itd_coverage_rle <- NA_character_
        itd_coverage_reads_processed <- 0L
        itd_coverage_reads_matched <- 0L
        itd_coverage_direct_status <- "unavailable_no_support_read_sequences"
      }
    }
  } else if (do_itd_coverage) {
    itd_coverage_status <- "unavailable_no_coordinate_template"
    itd_coverage_direct_status <- "unavailable_no_coordinate_template"
    itd_coverage_pe_status <- "unavailable_no_coordinate_template"
    itd_coverage_kmer_status <- "unavailable_no_coordinate_template"
    if (verbose) message(
      "[TALOS] ITD coverage unavailable: no full-length positional template ",
      "could be derived for breakpoint ", len_specific_bp,
      " and length ", best_len, " (reference interval ",
      genomic_start, "-", genomic_start + nchar(ref_dna) - 1L, ")."
    )
  }

  # [FIX v7] Remember the strict direct-alignment RLE here (the last point
  # where itd_coverage_rle still holds the direct-only vector). The display
  # fold below may overwrite itd_coverage_rle with the combined evidence.
  itd_coverage_rle_direct <- itd_coverage_rle
  if (do_itd_coverage && !is.na(coverage_template) && nzchar(coverage_template) &&
      !startsWith(itd_coverage_direct_status, "computed") &&
      !startsWith(itd_coverage_direct_status, "unavailable")) {
    itd_coverage_direct_status <- "unavailable_no_support_read_sequences"
  }

  # -- PE-spanning coverage of the ITD sequence ------------------------------
  # The BAM pair collection is requested whenever either PE length/discordant
  # metrics OR ITD coverage is enabled. QNAME linkage is required before a
  # pair can be attributed to the ITD-supporting read set; absent linkage is
  # unavailable (NA), not a measured 0% coverage result.
  if (do_itd_coverage && !is.na(coverage_template) &&
      nzchar(coverage_template) && !is_ptd_scale &&
      !is.na(best_len) && best_len > 0L &&
      nchar(coverage_template) >= best_len) {
    if (is.null(all_pairs) || length(all_pairs) == 0L) {
      itd_coverage_pe_status <- "no_pair_data"
      itd_pe_qname_status <- "no_pair_data"
    } else {
      pi_cov <- .pair_info(all_pairs, min_mapq = min_mapq)
      if (is.null(pi_cov)) {
        itd_coverage_pe_status <- "no_pair_data"
        itd_pe_qname_status <- "no_pair_data"
      } else {
        itd_pe_observation_pairs <- as.integer(pi_cov$n)
        itd_pe_valid_pairs <- as.integer(sum(pi_cov$valid, na.rm = TRUE))
        qname_match <- .match_support_qnames(pi_cov$qname, support_qnames)
        itd_pe_qname_usable_pairs <- qname_match$usable_pairs
        itd_pe_qname_matched_pairs <- qname_match$matched_pairs
        itd_pe_qname_status <- qname_match$status

        if (!identical(qname_match$status, "computed")) {
          itd_coverage_pe_status <- qname_match$status
          pe_cov <- list(cov_vec = integer(0), rle_str = NA_character_,
                         percent = NA_real_, percent_combined = NA_real_,
                         n_spanning_pairs = NA_integer_, n_support_pairs = NA_integer_)
        } else {
          valid_support <- pi_cov$valid & qname_match$matched &
            pi_cov$orientation %in% c("FR", "RF")
          valid_support[is.na(valid_support)] <- FALSE
          itd_pe_support_pairs <- as.integer(length(unique(qname_match$pair_qname[valid_support])))
          pe_cov <- .pe_spanning_coverage(
            qname = qname_match$pair_qname,
            s1 = pi_cov$s1, e1 = pi_cov$e1,
            s2 = pi_cov$s2, e2 = pi_cov$e2,
            orientation = pi_cov$orientation,
            valid = pi_cov$valid,
            support_qnames = .normalize_read_qname(support_qnames),
            breakpoint = len_specific_bp,
            span_len = as.integer(best_len),
            seq_len = nchar(coverage_template),
            direct_cov = cov_direct
          )
          itd_coverage_rle_pe <- pe_cov$rle_str
          cov_pe <- pe_cov$cov_vec
          itd_coverage_percent_pe <- pe_cov$percent
          itd_pe_spanning_pairs <- pe_cov$n_spanning_pairs
          itd_coverage_pe_status <- if (itd_pe_valid_pairs == 0L) {
            "computed_zero_valid_pairs_reference_projection_only"
          } else if (itd_pe_support_pairs == 0L) {
            "computed_no_valid_FR_RF_support_pairs_reference_projection_only"
          } else if (itd_pe_spanning_pairs == 0L) {
            "computed_no_spanning_support_pairs_reference_projection_only"
          } else {
            "computed_reference_coordinate_projection_not_allele_specific"
          }
          pe_summary <- .coverage_uncovered_summary(pe_cov$cov_vec)
          itd_uncovered_bases_pe <- pe_summary$bases
          itd_uncovered_intervals_pe <- pe_summary$intervals
          if (verbose && itd_pe_spanning_pairs > 0L) {
            message("[TALOS] PE-spanning ITD coverage: ", itd_pe_spanning_pairs,
                    " supporting pair(s) (", sprintf("%.1f%%", itd_coverage_percent_pe),
                    " of ITD bases)")
          }
        }
      }
    }
  } else if (do_itd_coverage && is_ptd_scale) {
    itd_coverage_pe_status <- "not_applicable_ptd_scale"
  } else if (do_itd_coverage && is.null(coverage_template)) {
    itd_coverage_pe_status <- "unavailable_no_coordinate_template"
  }

  # -- K-mer coverage of the ITD sequence ------------------------------------
  if (do_itd_coverage && !is.na(coverage_template) &&
      nzchar(coverage_template) && !is_ptd_scale &&
      !is.na(best_len) && best_len > 0L) {
    kmer_rows <- coverage_support_rows
    if (is.null(kmer_rows) || nrow(kmer_rows) == 0L) kmer_rows <- support_rows
    valid_mask <- !is.na(kmer_rows$read_seq) &
      nzchar(as.character(kmer_rows$read_seq))
    v_seqs_full <- as.character(kmer_rows$read_seq[valid_mask])
    v_qnames_full <- .normalize_read_qname(as.character(kmer_rows$read_name[valid_mask]))
    itd_kmer_tested_reads <- as.integer(length(v_seqs_full))
    itd_kmer_tested_fragments <- as.integer(length(unique(stats::na.omit(v_qnames_full))))
    if (length(v_seqs_full) > 0L) {
      km_cov <- .multi_kmer_itd_coverage(
        coverage_template, v_seqs_full,
        k_values = c(7L, 9L, 11L, 13L, 15L),
        min_offset_hits = c(3L, 2L, 2L, 2L, 2L)
      )
      cov_kmer <- km_cov$cov_vec
      cov_kmer_by_k <- lapply(km_cov$per_k, `[[`, "cov_vec")
      itd_coverage_rle_kmer <- km_cov$rle_str
      itd_coverage_rle_kmer_by_k <- lapply(km_cov$per_k, `[[`, "rle_str")
      itd_kmer_percent_by_k <- lapply(km_cov$per_k, `[[`, "percent")
      itd_kmer_reads_by_k <- lapply(km_cov$per_k, `[[`, "n_reads")
      itd_kmer_fragments_by_k <- lapply(km_cov$per_k, function(z) {
        matched <- z$matched_reads
        if (length(matched) != length(v_qnames_full)) return(NA_integer_)
        as.integer(length(unique(stats::na.omit(v_qnames_full[matched]))))
      })
      itd_kmer_reads <- km_cov$n_reads
      itd_kmer_matched_rate_pct <- if (itd_kmer_tested_reads > 0L) {
        100 * itd_kmer_reads / itd_kmer_tested_reads
      } else NA_real_
      matched_kmer_idx <- km_cov$matched_reads
      if (length(matched_kmer_idx) == length(v_qnames_full)) {
        itd_kmer_matched_fragments <- as.integer(length(unique(stats::na.omit(v_qnames_full[matched_kmer_idx]))))
      }
      itd_kmer_matched_fragment_rate_pct <- if (!is.na(itd_kmer_tested_fragments) && itd_kmer_tested_fragments > 0L) {
        100 * itd_kmer_matched_fragments / itd_kmer_tested_fragments
      } else NA_real_
      itd_coverage_percent_kmer <- km_cov$percent
      itd_coverage_kmer_status <- if (km_cov$n_reads > 0L) {
        "computed"
      } else {
        "computed_no_position_unique_kmer_matches"
      }
      km_summary <- .coverage_uncovered_summary(cov_kmer)
      itd_uncovered_bases_kmer <- km_summary$bases
      itd_uncovered_intervals_kmer <- km_summary$intervals
    } else {
      itd_coverage_kmer_status <- "unavailable_no_support_read_sequences"
      itd_coverage_rle_kmer <- NA_character_
      itd_coverage_percent_kmer <- NA_real_
    }
  } else if (do_itd_coverage && is_ptd_scale) {
    itd_coverage_kmer_status <- "not_applicable_ptd_scale"
  } else if (do_itd_coverage && is.null(coverage_template)) {
    itd_coverage_kmer_status <- "unavailable_no_coordinate_template"
  }

  # -- Display fold: direct read-sequence alignment + sequence-anchored k-mer
  if (do_itd_coverage && !is.na(coverage_template) && nzchar(coverage_template)) {
    n_template <- nchar(coverage_template)
    # PE coordinate projection is retained as a diagnostic only: it can
    # count a mate mapped to the original WT copy as if it covered the ITD
    # copy. Allele-specific combined coverage therefore uses actual read
    # sequence alignments (both mates) and sequence-anchored k-mers only.
    source_available <-
      startsWith(itd_coverage_direct_status, "computed") ||
      startsWith(itd_coverage_kmer_status, "computed")
    if (!source_available) {
      itd_coverage_status <- "unavailable_no_allele_specific_read_sequence_evidence"
      itd_coverage_percent_combined <- NA_real_
      itd_coverage_rle <- NA_character_
    } else {
      if (length(cov_direct) != n_template) cov_direct <- integer(0)
      if (length(cov_pe) != n_template) cov_pe <- integer(0)
      if (length(cov_kmer) != n_template) cov_kmer <- integer(0)
      disp <- .fold_display_coverage(
        cov_direct, integer(0), cov_kmer,
        direct_is_binary = direct_binary_fallback
      )
      if (length(disp) == n_template) {
        cov_display <- disp
        direct_summary <- .coverage_uncovered_summary(cov_direct)
        combined_summary <- .coverage_uncovered_summary(disp)
        itd_uncovered_bases_direct <- direct_summary$bases
        itd_uncovered_bases_combined <- combined_summary$bases
        itd_uncovered_intervals_combined <- combined_summary$intervals
        itd_coverage_percent_combined <- sum(disp > 0L) / n_template * 100
        itd_coverage_rle <- .rle_coverage_string(disp)
        itd_coverage_status <- "computed"
        if (length(cov_direct_certain) == n_template &&
            length(cov_direct_possible) == n_template) {
          certain_fold <- .fold_display_coverage(
            cov_direct_certain, integer(0), cov_kmer, direct_is_binary = FALSE
          )
          possible_fold <- .fold_display_coverage(
            cov_direct_possible, integer(0), cov_kmer, direct_is_binary = FALSE
          )
          if (length(certain_fold) == n_template && length(possible_fold) == n_template) {
            cov_display_certain <- certain_fold
            cov_display_possible <- possible_fold
            itd_coverage_combined_certain_pct <- sum(certain_fold > 0L) / n_template * 100
            itd_coverage_combined_possible_pct <- sum(possible_fold > 0L) / n_template * 100
            itd_coverage_rle_combined_certain <- .rle_coverage_string(certain_fold)
            itd_coverage_rle_combined_possible <- .rle_coverage_string(possible_fold)
          }
        }
      } else {
        itd_coverage_status <- "incomplete_source_vectors"
        itd_coverage_percent_combined <- NA_real_
        itd_coverage_rle <- NA_character_
      }
    }
  } else if (do_itd_coverage) {
    itd_coverage_status <- if (is_ptd_scale) "not_applicable_ptd_scale" else
      "unavailable_no_coordinate_template"
  }
  if (do_consistency && identical(coverage_template_source, "unavailable")) {
    support_consistency_status <- "unavailable_no_coordinate_template"
  }

  if (do_alignment_score && !is.na(itd_seq) && !imputed &&
      !is_ptd_scale && seg_in_ref) {
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
      }
    }
  }

  sequence_dependent_metrics_status <- if (is_ptd_scale) {
    "not_applicable_ptd_scale"
  } else if (is.na(itd_seq) || !nzchar(itd_seq)) {
    "unavailable_unresolved_full_length_sequence"
  } else if (isTRUE(imputed)) {
    "not_comparable_reference_imputed_full_length_sequence"
  } else {
    "full_length_read_derived_sequence_available"
  }
  ref_match_total_status <- if (!do_alignment_score) {
    "disabled"
  } else if (!is.na(ref_match_total)) {
    "computed"
  } else if (is.na(itd_seq) || !nzchar(itd_seq)) {
    "unavailable_unresolved_full_length_sequence"
  } else if (isTRUE(imputed)) {
    "not_comparable_reference_imputed_full_length_sequence"
  } else {
    "unavailable_no_valid_reference_comparison"
  }
  alignment_score_status <- ref_match_total_status

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
  raw_support_count <- vaf_metrics$raw_support
  depth_at_bp <- vaf_metrics$depth_at_bp
  raw_af <- vaf_metrics$raw_af
  corrected_af <- vaf_metrics$corrected_af
  size_correction_factor <- vaf_metrics$size_correction_factor
  wildtype_reads <- round(vaf_metrics$avg_wt)

  # Wild-type QNAME set carried on the result row (v0.4.20). Same strict
  # window (breakpoint +/- 10 bp, primary alignments) and same raw-string
  # setdiff semantics as .compute_left_right_coverage(), so its length
  # matches wildtype_reads; .merge_ptd_intervals() takes the exact union of
  # these sets across merged component rows instead of summing counts.
  wildtype_qnames <- character(0)
  if (!is.null(wildtype_info$gr) && length(wildtype_info$gr) > 0L &&
      length(len_specific_bp) == 1L && !is.na(len_specific_bp)) {
    widx <- wildtype_info$is_primary &
      GenomicRanges::start(wildtype_info$gr) <= (len_specific_bp - 10L) &
      GenomicRanges::end(wildtype_info$gr) >= (len_specific_bp + 10L)
    widx[is.na(widx)] <- FALSE
    wildtype_qnames <- setdiff(
      unique(stats::na.omit(as.character(wildtype_info$qnames[widx]))),
      support_qnames
    )
  }

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
      discordant_ratio <- compute_discordant_ratio(all_pairs, len_specific_bp, min_mapq = min_mapq)
    }

    # The all-pair and support-QNAME PE estimates share the same background
    # insert-size median; only their spanning-pair populations differ.
    pe_est <- compute_pe_itd_length(
      all_pairs, len_specific_bp, nominal_read_len = nominal_read_len,
      min_mapq = min_mapq
    )
    pe_est_support <- compute_pe_itd_length(
      all_pairs, len_specific_bp, nominal_read_len = nominal_read_len,
      min_mapq = min_mapq, support_qnames = support_qnames
    )
    length_pe <- pe_est$length_pe
    length_pe_nspanning <- as.integer(pe_est$n_spanning)
    length_pe_nbackground <- as.integer(pe_est$n_background)
    length_pe_nvalid <- as.integer(pe_est$n_valid %||% NA_integer_)
    length_pe_median_spanning <- as.numeric(pe_est$median_spanning %||% NA_real_)
    length_pe_median_background <- as.numeric(pe_est$median_background %||% NA_real_)
    length_pe_status <- pe_est$status %||% if (!is.na(length_pe)) "computed" else "insufficient_pairs"
    length_pe_support <- pe_est_support$length_pe
    length_pe_support_nspanning <- as.integer(pe_est_support$n_spanning)
    length_pe_support_nspanning_all <- as.integer(pe_est_support$n_spanning_all %||% NA_integer_)
    length_pe_support_nbackground <- as.integer(pe_est_support$n_background)
    length_pe_support_nvalid <- as.integer(pe_est_support$n_valid)
    length_pe_support_nqname_matched <- as.integer(pe_est_support$n_support_qname_matched %||% NA_integer_)
    length_pe_support_median_spanning <- as.numeric(pe_est_support$median_spanning %||% NA_real_)
    length_pe_support_median_background <- as.numeric(pe_est_support$median_background %||% NA_real_)
    length_pe_support_status <- pe_est_support$status %||% if (!is.na(length_pe_support)) "computed" else "insufficient_pairs"
    if (!is.na(length_pe_support) && !is.na(best_len) && best_len > 0L) {
      length_pe_support_delta <- as.numeric(length_pe_support - best_len)
      support_tolerance <- max(20, 0.15 * best_len)
      length_pe_support_concordant <- abs(length_pe_support_delta) <= support_tolerance
    }
    if (!is.na(length_pe) && !is.na(best_len) && best_len > 0L) {
      length_pe_delta <- as.numeric(length_pe - best_len)
      concordance_tolerance <- max(20, 0.15 * best_len)
      length_pe_concordant <- abs(length_pe_delta) <= concordance_tolerance
      if (verbose && !isTRUE(length_pe_concordant)) {
        message(
          "[TALOS] WARNING: anchor/k-mer length (", best_len,
          " bp) and PE insert-size estimate (", length_pe,
          " bp) disagree by ", abs(length_pe_delta),
          " bp; review pair geometry and breakpoint before reporting."
        )
      }
    }
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
      best_len = best_len,
      min_mapq = min_mapq
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

  if (do_detect_orientation && !is.na(itd_seq) && !imputed &&
      !is_ptd_scale && seg_in_ref) {
    local_start <- len_specific_bp - genomic_start + 1L
    local_end <- local_start + best_len - 1L
    if (!is.na(local_end) && local_end <= nchar(ref_dna)) {
      orientation <- detect_orientation(
        itd_seq,
        substr(ref_dna, local_start, local_end)
      )
    }
  }
  orientation_status <- if (!do_detect_orientation) {
    "disabled"
  } else if (!is.na(orientation)) {
    "computed"
  } else if (is.na(itd_seq) || !nzchar(itd_seq)) {
    "unavailable_unresolved_full_length_sequence"
  } else if (isTRUE(imputed)) {
    "not_comparable_reference_imputed_full_length_sequence"
  } else {
    "unavailable_no_orientation_match"
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
      raw_support_count,
      ", raw AF=",
      round(raw_af, 4),
      ", corrected AF=",
      round(corrected_af, 4),
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
    CoverageTemplateSource = coverage_template_source,
    CoverageTemplateLength = coverage_template_length,
    Length = best_len,
    InFrame = if (!is.na(best_len) && best_len > 0L) (as.integer(best_len) %% 3L == 0L) else NA,
    DupExons = dup_exons,
    CandidateTypes = candidate_types,
    SupportingReads = as.integer(raw_support_count),
    AdjustedSupportingReads = as.integer(corrected_support),
    SizeCorrectionFactor = size_correction_factor,
    RawAlleleFrequency = raw_af,
    CorrectedAlleleFrequency = corrected_af,
    IndependentSupportQNAMEs = as.integer(raw_support),
    AnchorConfirmedQNAMEs = as.integer(anchor_confirmed_qnames),
    AnchorConfirmedQNAMEsStatus = if (is.na(anchor_confirmed_qnames)) "not_available" else "computed_from_extension_anchor_source",
    WildtypeReads = wildtype_reads,
    DepthAtBreakpoint = depth_at_bp,
    CorrectedDepthAtBreakpoint = as.integer(vaf_metrics$corrected_depth),
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
    SequenceSource = if (sequence_unresolved) {
      "partial_read_unresolved"
    } else if (!imputed) {
      "observed"
    } else if (sequence_partial && !is.na(itd_seq)) {
      "partial_read+ref"
    } else if (sequence_partial) {
      "partial_read_unresolved"
    } else {
      "ref_imputed"
    },
    RefMatch_Observed = ref_match_observed,
    RefMatch_Total = ref_match_total,
    RefMatchObservedStatus = if (!is.na(ref_match_observed)) "computed_from_observed_fragment" else
      if (identical(itd_observed_fragment_status, "not_applicable_ptd_scale")) "not_applicable_ptd_scale" else
        itd_observed_fragment_status,
    RefMatchTotalStatus = ref_match_total_status,
    AlignmentScoreStatus = alignment_score_status,
    OrientationStatus = orientation_status,
    SequenceDependentMetricsStatus = sequence_dependent_metrics_status,
    SupportConsistencyStatus = support_consistency_status,
    ITDReadCoverage = itd_coverage_percent,
    ITDReadCoverage_Direct = itd_coverage_percent,
    ITDReadCoverage_PE = itd_coverage_percent_pe,
    ITDCoverageStatus = itd_coverage_status,
    ITDCoverageDirectStatus = itd_coverage_direct_status,
    ITDCoveragePEStatus = itd_coverage_pe_status,
    ITDCoverageKmerStatus = itd_coverage_kmer_status,
    ITDCoverageRLE = itd_coverage_rle,
    ITDCoverageRLE_Direct = itd_coverage_rle_direct,
    ITDCoverageRLE_PE = itd_coverage_rle_pe,
    ITDCoverageRLE_Kmer = itd_coverage_rle_kmer,
    ITDKmerAmbiguousCoveragePercent = if (exists("km_cov", inherits = FALSE) && !is.null(km_cov$ambiguous_percent)) km_cov$ambiguous_percent else NA_real_,
    ITDKmerAmbiguousCoverageRLE = if (exists("km_cov", inherits = FALSE) && length(km_cov$ambiguous_cov_vec)) .rle_coverage_string(km_cov$ambiguous_cov_vec) else NA_character_,
    ITDReadCoverage_Combined = itd_coverage_percent_combined,
    ITDReadCoverage_Kmer = itd_coverage_percent_kmer,
    ITDReadCoverage_DirectCertain = itd_coverage_direct_certain_pct,
    ITDReadCoverage_DirectPossible = itd_coverage_direct_possible_pct,
    ITDReadCoverage_DirectAmbiguous = itd_coverage_direct_ambiguous_pct,
    ITDReadCoverage_CombinedCertain = itd_coverage_combined_certain_pct,
    ITDReadCoverage_CombinedPossible = itd_coverage_combined_possible_pct,
    ITDDirectAmbiguousReadEnds = itd_coverage_direct_ambiguous_read_ends,
    ITDDirectMultiPlacementReadEnds = itd_coverage_direct_multi_placement_read_ends,
    ITDDirectAlternativeSearchLimitReadEnds = itd_coverage_direct_search_limit_read_ends,
    ITDDirectMeanPlacements = itd_coverage_direct_mean_placements,
    ITDDirectMaxPlacements = itd_coverage_direct_max_placements,
    ITDKmerKValues = itd_kmer_k_values,
    ITDReadCoverage_Kmer7 = itd_kmer_percent_by_k[["7"]] %||% NA_real_,
    ITDReadCoverage_Kmer9 = itd_kmer_percent_by_k[["9"]] %||% NA_real_,
    ITDReadCoverage_Kmer11 = itd_kmer_percent_by_k[["11"]] %||% NA_real_,
    ITDReadCoverage_Kmer13 = itd_kmer_percent_by_k[["13"]] %||% NA_real_,
    ITDReadCoverage_Kmer15 = itd_kmer_percent_by_k[["15"]] %||% NA_real_,
    ITDKmerMatchedReadEnds7 = itd_kmer_reads_by_k[["7"]] %||% NA_integer_,
    ITDKmerMatchedReadEnds9 = itd_kmer_reads_by_k[["9"]] %||% NA_integer_,
    ITDKmerMatchedReadEnds11 = itd_kmer_reads_by_k[["11"]] %||% NA_integer_,
    ITDKmerMatchedReadEnds13 = itd_kmer_reads_by_k[["13"]] %||% NA_integer_,
    ITDKmerMatchedReadEnds15 = itd_kmer_reads_by_k[["15"]] %||% NA_integer_,
    ITDKmerMatchedFragments7 = itd_kmer_fragments_by_k[["7"]] %||% NA_integer_,
    ITDKmerMatchedFragments9 = itd_kmer_fragments_by_k[["9"]] %||% NA_integer_,
    ITDKmerMatchedFragments11 = itd_kmer_fragments_by_k[["11"]] %||% NA_integer_,
    ITDKmerMatchedFragments13 = itd_kmer_fragments_by_k[["13"]] %||% NA_integer_,
    ITDKmerMatchedFragments15 = itd_kmer_fragments_by_k[["15"]] %||% NA_integer_,
    ITDCoverageRLE_DirectCertain = itd_coverage_rle_direct_certain,
    ITDCoverageRLE_DirectPossible = itd_coverage_rle_direct_possible,
    ITDCoverageRLE_DirectAmbiguous = itd_coverage_rle_direct_ambiguous,
    ITDCoverageRLE_CombinedCertain = itd_coverage_rle_combined_certain,
    ITDCoverageRLE_CombinedPossible = itd_coverage_rle_combined_possible,
    ITDCoverageRLE_Kmer7 = itd_coverage_rle_kmer_by_k[["7"]] %||% NA_character_,
    ITDCoverageRLE_Kmer9 = itd_coverage_rle_kmer_by_k[["9"]] %||% NA_character_,
    ITDCoverageRLE_Kmer11 = itd_coverage_rle_kmer_by_k[["11"]] %||% NA_character_,
    ITDCoverageRLE_Kmer13 = itd_coverage_rle_kmer_by_k[["13"]] %||% NA_character_,
    ITDCoverageRLE_Kmer15 = itd_coverage_rle_kmer_by_k[["15"]] %||% NA_character_,
    ITDCoverageReadsProcessed = itd_coverage_reads_processed,
    ITDCoverageMateReadEndsAdded = itd_coverage_mate_read_ends_added,
    ITDCoverageReadsMatched = itd_coverage_reads_matched,
    ITDCoverageReadAlignmentsProcessed = itd_coverage_reads_processed,
    ITDCoverageFragmentsProcessed = itd_coverage_fragments_processed,
    ITDCoverageFragmentsWithDirectAlignment = itd_coverage_fragments_matched,
    ITDCoverageDirectAlignmentRatePct = itd_coverage_direct_alignment_rate_pct,
    ITDCoverageDirectMedianPID = itd_coverage_direct_median_pid,
    ITDCoverageDirectMedianAlignedBases = itd_coverage_direct_median_aligned_bases,
    ITDCoverageDirectMeanPID = itd_coverage_direct_mean_pid,
    ITDCoverageDirectMedianQueryCoveragePct = itd_coverage_direct_median_query_coverage_pct,
    ITDCoverageDirectReadsAnyAlignment = itd_coverage_direct_n_any_alignment,
    ITDCoverageDirectReadsQualifying = itd_coverage_direct_n_qualifying_alignment,
    ITDObservedFragmentPID = itd_observed_fragment_pid,
    ITDObservedFragmentQueryCoveragePct = itd_observed_fragment_query_coverage_pct,
    ITDObservedFragmentAlignedBases = itd_observed_fragment_aligned_bases,
    ITDObservedFragmentOrientation = itd_observed_fragment_orientation,
    ITDObservedFragmentReferenceStart = itd_observed_fragment_reference_start,
    ITDObservedFragmentReferenceEnd = itd_observed_fragment_reference_end,
    ITDObservedFragmentTemplateStart = itd_observed_fragment_reference_start,
    ITDObservedFragmentTemplateEnd = itd_observed_fragment_reference_end,
    ITDObservedFragmentSupportReads = itd_observed_fragment_support_reads,
    ITDObservedFragmentUniqueSequences = itd_observed_fragment_unique_sequences,
    ITDObservedFragmentStatus = itd_observed_fragment_status,
    ITDPE_ObservationPairs = itd_pe_observation_pairs,
    ITDPE_ValidPairs = itd_pe_valid_pairs,
    ITDPE_SupportPairs = itd_pe_support_pairs,
    ITDPE_QnameUsablePairs = itd_pe_qname_usable_pairs,
    ITDPE_QnameMatchedPairs = itd_pe_qname_matched_pairs,
    ITDPE_QnameMatchStatus = itd_pe_qname_status,
    ITDCoverageDirectCoveredBases = .coverage_depth_summary(cov_direct)$covered_bases,
    ITDCoverageDirectTotalDepth = .coverage_depth_summary(cov_direct)$total_depth,
    ITDCoverageDirectMeanDepth = .coverage_depth_summary(cov_direct)$mean_depth,
    ITDCoverageDirectMaxDepth = .coverage_depth_summary(cov_direct)$max_depth,
    ITDCoveragePECoveredBases = .coverage_depth_summary(cov_pe)$covered_bases,
    ITDCoveragePETotalDepth = .coverage_depth_summary(cov_pe)$total_depth,
    ITDCoveragePEMeanDepth = .coverage_depth_summary(cov_pe)$mean_depth,
    ITDCoveragePEMaxDepth = .coverage_depth_summary(cov_pe)$max_depth,
    ITDCoverageKmerCoveredBases = .coverage_depth_summary(cov_kmer)$covered_bases,
    ITDCoverageKmerTotalDepth = .coverage_depth_summary(cov_kmer)$total_depth,
    ITDCoverageKmerMeanDepth = .coverage_depth_summary(cov_kmer)$mean_depth,
    ITDCoverageKmerMaxDepth = .coverage_depth_summary(cov_kmer)$max_depth,
    ITDCoverageCombinedCoveredBases = .coverage_depth_summary(cov_display)$covered_bases,
    ITDCoverageCombinedTotalDepth = .coverage_depth_summary(cov_display)$total_depth,
    ITDCoverageCombinedMeanDepth = .coverage_depth_summary(cov_display)$mean_depth,
    ITDCoverageCombinedMaxDepth = .coverage_depth_summary(cov_display)$max_depth,
    ITDCoverageCombinedDepthMethod = "per_base_max_non_additive",
    ITDKmerTestedReads = itd_kmer_tested_reads,
    ITDKmerMatchedReadPct = itd_kmer_matched_rate_pct,
    ITDKmerTestedFragments = itd_kmer_tested_fragments,
    ITDKmerMatchedFragments = itd_kmer_matched_fragments,
    ITDKmerMatchedFragmentPct = itd_kmer_matched_fragment_rate_pct,
    ITDUncoveredBases_Direct = itd_uncovered_bases_direct,
    ITDUncoveredIntervals_Direct = .coverage_uncovered_summary(cov_direct)$intervals,
    ITDUncoveredBases_PE = itd_uncovered_bases_pe,
    ITDUncoveredIntervals_PE = itd_uncovered_intervals_pe,
    ITDUncoveredBases_Kmer = itd_uncovered_bases_kmer,
    ITDUncoveredIntervals_Kmer = itd_uncovered_intervals_kmer,
    ITDUncoveredBases_Combined = itd_uncovered_bases_combined,
    ITDUncoveredIntervals_Combined = itd_uncovered_intervals_combined,
    ITDKmerReads = itd_kmer_reads,
    ITDPE_SpanningPairs = itd_pe_spanning_pairs,
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
    LengthPE_NBackground = length_pe_nbackground,
    LengthPE_NValid = length_pe_nvalid,
    LengthPE_MedianSpanningInsert = length_pe_median_spanning,
    LengthPE_MedianBackgroundInsert = length_pe_median_background,
    LengthPEStatus = length_pe_status,
    LengthPE_Delta = length_pe_delta,
    LengthPE_Concordant = length_pe_concordant,
    LengthPE_AllPairs = length_pe,
    LengthPE_AllPairs_NSpanning = length_pe_nspanning,
    LengthPE_AllPairs_MedianSpanningInsert = length_pe_median_spanning,
    LengthPE_AllPairs_MedianBackgroundInsert = length_pe_median_background,
    LengthPE_AllPairs_Status = length_pe_status,
    LengthPE_AllPairs_Eligibility = "Valid primary pairs; MAPQ >= 20; outer insert <= 5000 bp; >=5 strict spanning pairs (pair starts <= breakpoint - read length and ends >= breakpoint + read length); >=10 background pairs",
    LengthPE_SupportQNAME = length_pe_support,
    LengthPE_SupportQNAME_NSpanning = length_pe_support_nspanning,
    LengthPE_SupportQNAME_NSpanningAll = length_pe_support_nspanning_all,
    LengthPE_SupportQNAME_NBackground = length_pe_support_nbackground,
    LengthPE_SupportQNAME_NValid = length_pe_support_nvalid,
    LengthPE_SupportQNAME_NMatched = length_pe_support_nqname_matched,
    LengthPE_SupportQNAME_MedianSpanningInsert = length_pe_support_median_spanning,
    LengthPE_SupportQNAME_MedianBackgroundInsert = length_pe_support_median_background,
    LengthPE_SupportQNAME_Status = length_pe_support_status,
    LengthPE_SupportQNAME_Eligibility = "Same strict pair filters as all-pair estimate, restricted to candidate-support QNAMEs; requires usable pair QNAMEs and >=5 strict spanning support pairs",
    LengthPE_SupportQNAME_Delta = length_pe_support_delta,
    LengthPE_SupportQNAME_Concordant = length_pe_support_concordant,
    LengthPE_SupportInsert = as.numeric(pe_est_support$length_pe_support_insert %||% NA_real_),
    LengthPE_SupportInsert_NPairs = as.integer(pe_est_support$n_support_insert_pairs %||% NA_integer_),
    LengthPE_SupportInsert_MedianInsert = as.numeric(pe_est_support$median_support_insert %||% NA_real_),
    LengthPE_SupportInsert_Status = as.character(pe_est_support$support_insert_status %||% "not_computed"),
    LengthPE_SupportInsert_Eligibility = "Candidate-support QNAME pairs whose outer fragment spans the breakpoint; >=5 pairs and valid background median; relaxed geometry, diagnostic only",
    ITDTerminalWindowStart = as.numeric(terminal_start_genomic),
    ITDTerminalWindowEnd = as.numeric(terminal_end_genomic),
    ITDTerminalSupportAlignments = itd_terminal_supporting_alignments,
    ITDTerminalUniqueQNAMEs = itd_terminal_unique_qnames,
    ITDTerminalReverseAlignments = itd_terminal_reverse_alignments,
    ITDTerminalSoftclipAlignments = itd_terminal_cigar_softclip_alignments,
    ITDTerminalBaseQualityStatus = itd_terminal_base_quality_status,
    ITDBaseQualityStatus = itd_base_quality_status,
    ITDMeanBaseQuality = itd_mean_base_quality,
    ITDMinBaseQuality = itd_min_base_quality,
    ITDTerminalMeanBaseQuality = itd_terminal_mean_base_quality,
    ITDAdapterScreeningStatus = itd_adapter_screening_status,
    ITDAdapterScreenedReads = as.integer(adapter_screened_reads),
    ITDAdapterPositiveReads = as.integer(adapter_positive_reads),
    ITDAdapterPositiveSupportReads = adapter_support_hits,
    ITDTerminalReadEvidence = terminal_evidence_field,
    ITDTerminalReadEvidenceTruncated = terminal_evidence_truncated,
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
    SupportQNAMEs = .qname_set_to_string(support_qnames),
    WildtypeQNAMEs = .qname_set_to_string(wildtype_qnames),
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
  th <- function(name) {
    v <- thresholds[[name]]
    if (is.null(v) || length(v) != 1L || is.na(v)) {
      return(if (startsWith(name, "max")) Inf else 0)
    }
    as.numeric(v)
  }

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
    "InFrame",
    "LengthExt",
    "SupportingReads",
    "WildtypeReads",
    "DepthAtBreakpoint",
    "CorrectedDepthAtBreakpoint",
    "AlleleFrequency",
    "RawAlleleFrequency",
    "CorrectedAlleleFrequency",
    "SizeCorrectionFactor",
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
    "ITDReadCoverage_Combined",
    "ITDReadCoverage_Kmer",
    "ITDUncoveredBases_Direct",
    "ITDUncoveredBases_Combined",
    "ITDUncoveredIntervals_Combined",
    "ITDKmerReads",
    "ITDKmerTestedReads",
    "ITDKmerMatchedReadPct",
    "ITDKmerTestedFragments",
    "ITDKmerMatchedFragments",
    "ITDKmerMatchedFragmentPct",
    "ITDPE_SpanningPairs",
    "ITDMeanBaseQuality",
    "ITDAdapterPositiveSupportReads",
    "LengthPE_Delta",
    "LengthPE_Concordant",
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

#' Load ALU consensus sequences from a FASTA file
#'
#' Expects a FASTA with entries named after ALU subfamilies, e.g.:
#'   >AluSx
#'   GGCCGGGCGCGGTGGCTCAC...
#'   >AluY
#'   ...
#'
#' Ships with the package at inst/extdata/alu_consensus.fa.
#' Users can supply a custom FASTA via `consensus_fa`.
#'
#' @param consensus_fa Path to FASTA (default: bundled alu_consensus.fa).
#' @return Named DNAStringSet.
#' @keywords internal
.load_alu_consensus <- function(consensus_fa = NULL) {
  if (is.null(consensus_fa)) {
    consensus_fa <- system.file(
      "extdata",
      "alu_consensus.fa",
      package = "TALOS"
    )
  }
  if (!file.exists(consensus_fa)) {
    stop(
      "ALU consensus FASTA not found at: ",
      consensus_fa,
      "\n",
      "  Supply a custom path via consensus_fa= or install the full TALOS package."
    )
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings required for ALU sequence matching.")
  }
  Biostrings::readDNAStringSet(consensus_fa)
}


#' Map the soft-clip side to the reference breakpoint coordinate
#'
#' Coordinates returned by GAlignments are 1-based. A leading soft clip
#' precedes the first aligned reference base, so its insertion breakpoint
#' is immediately before read_start (read_start - 1). A trailing soft clip
#' follows the last aligned reference base, so its breakpoint is read_end.
#'
#' [FIX Finding 2] Chromosome-boundary guard: if side == "lead" and
#' read_start == 1, read_start - 1 == 0 is not a valid genomic coordinate.
#' Return NA_integer_ in that case so the candidate can be excluded before
#' clustering and annotation, rather than carrying a 0 into the pipeline.
#'
#' @param side Character scalar: "lead" or "trail".
#' @param read_start 1-based start of the aligned read.
#' @param read_end 1-based end of the aligned read.
#' @return Integer genomic breakpoint coordinate, or NA_integer_ when the
#'   leading-clip breakpoint would fall before position 1.
#' @keywords internal
.alu_clip_breakpoint <- function(side, read_start, read_end) {
  if (
    length(side) != 1L ||
      is.na(side) ||
      !side %in% c("lead", "trail")
  ) {
    stop("side must be 'lead' or 'trail'.", call. = FALSE)
  }

  if (
    length(read_start) != 1L ||
      is.na(read_start) ||
      length(read_end) != 1L ||
      is.na(read_end)
  ) {
    stop("read_start and read_end must be non-missing scalars.", call. = FALSE)
  }

  if (side == "lead") {
    bp <- as.integer(read_start - 1L)
    # [FIX Finding 2] Chromosome-boundary guard: a leading clip at
    # read_start == 1 produces bp = 0, which is not a valid 1-based genomic
    # coordinate. Return NA so .extract_alu_candidates() can exclude it.
    if (bp < 1L) {
      return(NA_integer_)
    }
    bp
  } else {
    as.integer(read_end)
  }
}


#' Match a single soft-clip sequence against all ALU consensus entries
#'
#' Uses local pairwise alignment (Smith-Waterman).  Returns a list with:
#'   subtype       - best-matching ALU family name
#'   score         - normalised alignment score (0-1)
#'   strand        - "+" (sense) or "-" (antisense)
#'   aln_start     - start position in the ALU consensus (for 5' truncation calc)
#'   aln_end       - end position in the ALU consensus
#'   pattern_start - start position of the local alignment in the *clip*
#'                   sequence (1-based). Bases before this position in the clip
#'                   are not aligned to ALU and may be target-site sequence.
#'   pattern_end   - end position of the local alignment in the *clip*
#'                   sequence. Bases after this position in the clip are not
#'                   aligned to ALU and may be target-site sequence or poly-A.
#'
#' [ITEM 2] Pattern alignment positions are now stored alongside subject
#' (consensus) positions. The non-ALU bases at clip boundaries
#' (clip[1 .. pattern_start-1] for trail clips, clip[pattern_end+1 .. end]
#' for lead clips) frequently correspond to the target-site duplication
#' (TSD) created at insertion time. Storing these positions enables more
#' precise TSD detection without re-running the alignment. See
#' .extract_alu_boundary_bases() for the helper that recovers them.
#'
#' @param clip_seq  Character string (soft-clip bases).
#' @param alu_seqs  Named DNAStringSet from .load_alu_consensus().
#' @param min_clip_len Minimum clip length to attempt alignment (default 25).
#' @param min_score    Minimum normalised score to accept a match (default 0.6).
#' @param min_aligned_bases Minimum aligned clip width required (default 20).
#' @keywords internal
.match_clip_to_alu <- function(
  clip_seq, alu_seqs, min_clip_len = 25L, min_score = 0.60,
  min_aligned_bases = 20L
) {
  if (length(clip_seq) != 1L || is.na(clip_seq) || nchar(clip_seq) < min_clip_len ||
      !requireNamespace("Biostrings", quietly = TRUE) ||
      !requireNamespace("pwalign", quietly = TRUE)) return(NULL)
  clip_seq <- toupper(as.character(clip_seq))
  n <- nchar(clip_seq)
  clip_dna <- Biostrings::DNAString(clip_seq)
  clip_rc <- Biostrings::reverseComplement(clip_dna)
  best <- NULL; best_score <- -Inf
  orientations <- list(`+` = clip_dna, `-` = clip_rc)
  for (nm in names(alu_seqs)) {
    ref <- alu_seqs[[nm]]
    for (strand in names(orientations)) {
      aln <- tryCatch(pwalign::pairwiseAlignment(
        orientations[[strand]], ref, type = "local", gapOpening = -10, gapExtension = -0.5
      ), error = function(e) NULL)
      if (is.null(aln)) next
      pat <- pwalign::pattern(aln); subj <- pwalign::subject(aln)
      ps <- as.integer(pwalign::start(pat)); pe <- as.integer(pwalign::end(pat))
      ss <- as.integer(pwalign::start(subj)); se <- as.integer(pwalign::end(subj))
      if (anyNA(c(ps, pe, ss, se)) || pe < ps) next
      aligned_bases <- pe - ps + 1L
      if (aligned_bases < min_aligned_bases) next
      norm_score <- as.numeric(pwalign::score(aln)) / aligned_bases
      if (!is.finite(norm_score) || norm_score <= best_score) next
      if (strand == "-") {
        ps_fwd <- n - pe + 1L; pe_fwd <- n - ps + 1L
      } else {
        ps_fwd <- ps; pe_fwd <- pe
      }
      best_score <- norm_score
      best <- list(subtype = nm, score = norm_score, strand = strand,
                   aln_start = ss, aln_end = se, pattern_start = ps_fwd,
                   pattern_end = pe_fwd, aligned_bases = aligned_bases,
                   clip_fraction = aligned_bases / n)
    }
  }
  if (is.null(best) || best$score < min_score) return(NULL)
  best
}


#' Extract the non-ALU boundary bases from a soft-clip sequence
#'
#' [ITEM 2] Helper that uses pattern_start / pattern_end stored by
#' .match_clip_to_alu() to recover the bases at the clip boundary that did
#' NOT align to the ALU consensus. For a freshly inserted ALU with a TSD
#' (inserted material layout: upstream_TSD + ALU + downstream_TSD, in the
#' donor genome between two reference positions):
#'
#'   - Trail clip (read crosses the UPSTREAM junction; aligned portion is
#'     the upstream reference, 3' end of the read is the inserted material):
#'       clip = upstream_TSD + ALU_start
#'       The ALU pattern alignment starts at pattern_start, so the bases
#'       clip[1 .. pattern_start-1] are candidate upstream TSD.
#'
#'   - Lead clip (read crosses the DOWNSTREAM junction; aligned portion is
#'     the downstream reference, 5' end of the read is the inserted material):
#'       clip = ALU_end + downstream_TSD
#'       The ALU pattern alignment ends at pattern_end, so the bases
#'       clip[pattern_end+1 .. nchar(clip)] are candidate downstream TSD.
#'
#' poly-A tails appended to ALU 3' ends also land in the lead-clip suffix;
#' callers can strip a terminal poly-A first via .detect_poly_a() if the
#' read orientation makes this ambiguous.
#'
#' .match_clip_to_alu() maps reverse-orientation alignment coordinates back
#' to the original forward clip, so these boundary coordinates are always
#' interpreted in the original clip sequence.
#'
#' @param clip_seq     Character - the soft-clip sequence.
#' @param pattern_start Integer - start of ALU alignment in clip (1-based).
#' @param pattern_end   Integer - end of ALU alignment in clip (1-based).
#' @param clip_side    Character - "lead" or "trail".
#' @return Character scalar with the boundary bases ("" if none), or
#'   NA_character_ if inputs are missing / invalid.
#' @keywords internal
.extract_alu_boundary_bases <- function(
  clip_seq,
  pattern_start,
  pattern_end,
  clip_side
) {
  if (
    is.na(clip_seq) || !nzchar(clip_seq) ||
      is.na(pattern_start) || is.na(pattern_end) ||
      !clip_side %in% c("lead", "trail")
  ) {
    return(NA_character_)
  }
  n <- nchar(clip_seq)
  if (pattern_start < 1L || pattern_end > n || pattern_start > pattern_end) {
    return(NA_character_)
  }

  if (clip_side == "trail") {
    # Trail clip = upstream_TSD + ALU_start. TSD is the PREFIX.
    if (pattern_start <= 1L) return("")
    substr(clip_seq, 1L, pattern_start - 1L)
  } else {
    # Lead clip = ALU_end + downstream_TSD. TSD is the SUFFIX.
    if (pattern_end >= n) return("")
    substr(clip_seq, pattern_end + 1L, n)
  }
}


#' Detect poly-A tail in a soft-clip sequence
#'
#' An ALU 3'-end clip shows a stretch of A (sense) or T (antisense).
#' Returns the length of the longest run found (0 if absent).
#'
#' @param clip_seq Character string.
#' @param min_run  Minimum consecutive A/T to report (default 6).
#' @keywords internal
.detect_poly_a <- function(clip_seq, min_run = 6L) {
  min_run <- suppressWarnings(as.integer(min_run))
  if (length(min_run) != 1L || is.na(min_run) || min_run < 1L ||
      length(clip_seq) != 1L || is.na(clip_seq) || !nzchar(clip_seq)) {
    return(0L)
  }

  clip_seq <- toupper(as.character(clip_seq))

  # A poly-A/T tail is terminal, not merely an A/T-rich run somewhere
  # inside the soft clip. Searching the whole clip used to annotate internal
  # ALU bases or unrelated low-complexity sequence as a tail.
  terminal_run <- function(base, seq) {
    pattern <- paste0(base, "{", min_run, ",}")
    left <- regexpr(paste0("^", pattern), seq, perl = TRUE)
    right <- regexpr(paste0(pattern, "$"), seq, perl = TRUE)
    left_n <- if (left[1L] < 0L) 0L else attr(left, "match.length")[[1L]]
    right_n <- if (right[1L] < 0L) 0L else attr(right, "match.length")[[1L]]
    as.integer(max(left_n, right_n))
  }

  best <- max(terminal_run("A", clip_seq), terminal_run("T", clip_seq))
  if (best < min_run) 0L else as.integer(best)
}


#' Estimate 5' truncation of an ALU insertion
#'
#' New ALU copies are frequently 5'-truncated. If the clip matches the
#' ALU consensus starting at position > 1, the unmatched 5'-end is lost.
#'
#' @param aln_start Integer - alignment start in consensus (1-based).
#' @param consensus_len Integer - total consensus length.
#' @keywords internal
.alu_5p_truncation <- function(aln_start, consensus_len) {
  if (is.na(aln_start) || is.na(consensus_len)) {
    return(NA_integer_)
  }
  as.integer(max(0L, aln_start - 1L))
}


#' Detect target-site duplication (TSD) from read evidence
#'
#' [FIX Finding 1] The previous implementation compared the upstream suffix
#' and downstream prefix of the REFERENCE sequence and reported their longest
#' match as the TSD. This was biologically unreliable: a genuine TSD is
#' created during insertion, and the reference genome contains the original
#' target sequence only once. Reference-flank comparison can only detect
#' pre-existing local repeats, not a TSD created by insertion.
#'
#' New approach: use read evidence. A genuine TSD causes the trail-clip
#' breakpoints (reads from upstream, aligning through the upstream TSD copy)
#' and lead-clip breakpoints (reads from downstream, aligning from the
#' downstream TSD copy) to differ by exactly the TSD length:
#'
#'   trail_bp  = end of TSD in reference (last aligned base of upstream reads)
#'   lead_bp   = start of TSD - 1 in reference (base before first aligned base
#'               of downstream reads)
#'   TSD_length = trail_bp - lead_bp
#'   TSD_sequence = reference[lead_bp + 1 .. trail_bp]
#'
#' When both sides are present and the gap is in the expected TSD range
#' (5-20 bp), the TSD is reported as "read_supported". When only one side
#' is present, or the gap is outside the expected range, the TSD is reported
#' as "unknown" rather than falling back to the unreliable reference-flank
#' comparison.
#'
#' @param ref_dna       Character - full reference sequence for the region.
#' @param ins_pos       Integer - insertion site (1-based, local coordinates).
#' @param cluster_df    data.frame of supporting reads (from
#'                      .extract_alu_candidates). Must contain columns
#'                      genomic_pos and clip_side. May be NULL.
#' @param genomic_start Integer - genomic start of ref_dna window, for
#'                      converting genomic_pos to local coordinates.
#' @param max_tsd       Integer - maximum TSD length to accept (default 20).
#' @param min_tsd       Integer - minimum TSD length to accept (default 5).
#' @return Named list: tsd_len (integer), tsd_seq (character),
#'   tsd_status (character: "read_supported" or "unknown").
#' @keywords internal
.detect_alu_tsd <- function(
  ref_dna,
  ins_pos,
  cluster_df = NULL,
  genomic_start = NA_integer_,
  max_tsd = 20L,
  min_tsd = 5L
) {
  ref_len <- nchar(ref_dna)
  if (is.na(ins_pos) || ins_pos < 1L || ins_pos > ref_len) {
    return(list(tsd_len = NA_integer_, tsd_seq = NA_character_,
                tsd_status = "unknown"))
  }

  # Without read evidence we cannot establish a TSD. Report as unknown.
  if (is.null(cluster_df) || nrow(cluster_df) == 0L) {
    return(list(tsd_len = NA_integer_, tsd_seq = NA_character_,
                tsd_status = "unknown"))
  }

  trail_bps <- cluster_df$genomic_pos[
    !is.na(cluster_df$clip_side) & cluster_df$clip_side == "trail"
  ]
  lead_bps <- cluster_df$genomic_pos[
    !is.na(cluster_df$clip_side) & cluster_df$clip_side == "lead"
  ]

  # Exclude NA / out-of-range breakpoints (e.g. chromosome-boundary reads
  # flagged by .alu_clip_breakpoint()).
  trail_bps <- trail_bps[!is.na(trail_bps) & trail_bps > 0L]
  lead_bps  <- lead_bps[!is.na(lead_bps) & lead_bps > 0L]

  if (length(trail_bps) == 0L || length(lead_bps) == 0L) {
    # Only one side present -- cannot determine TSD from read evidence.
    return(list(tsd_len = NA_integer_, tsd_seq = NA_character_,
                tsd_status = "unknown"))
  }

  trail_bp <- as.integer(round(stats::median(trail_bps)))
  lead_bp  <- as.integer(round(stats::median(lead_bps)))

  tsd_len <- trail_bp - lead_bp

  if (is.na(tsd_len) || tsd_len < min_tsd || tsd_len > max_tsd) {
    # Gap outside the expected TSD range. Report as unknown.
    return(list(tsd_len = NA_integer_, tsd_seq = NA_character_,
                tsd_status = "unknown"))
  }

  # Extract the TSD sequence from the reference.
  # Convert genomic breakpoints to local coordinates.
  local_start <- lead_bp + 1L - genomic_start + 1L
  local_end   <- trail_bp - genomic_start + 1L

  if (local_start < 1L || local_end > ref_len || local_start > local_end) {
    return(list(tsd_len = NA_integer_, tsd_seq = NA_character_,
                tsd_status = "unknown"))
  }

  tsd_seq <- substr(ref_dna, local_start, local_end)

  list(
    tsd_len = as.integer(tsd_len),
    tsd_seq = tsd_seq,
    tsd_status = "read_supported"
  )
}


#' Extract ALU-positive soft-clip candidates from a set of reads
#'
#' For each read, both the leading and trailing soft-clips are checked
#' against every ALU consensus.  Reads with at least one matching clip
#' are retained with the best-hit metadata.
#'
#' @param reads       GAlignments object (from .load_bam_data_streaming).
#' @param alu_seqs    Named DNAStringSet.
#' @param genomic_start Integer - regional offset for local-coord conversion.
#' @param min_clip_len  Minimum clip length (default 25 bp).
#' @param min_alu_score Minimum normalised alignment score (default 0.60).
#' @param verbose       Print progress (default FALSE).
#' @return data.frame with columns:
#'   read_name, genomic_pos, local_pos, clip_side, clip_seq,
#'   alu_subtype, alu_score, alu_strand, aln_start, aln_end,
#'   pattern_start, pattern_end, boundary_bases, poly_a_len, mapq, is_reverse
#' @keywords internal
.extract_alu_candidates <- function(
  reads,
  alu_seqs,
  genomic_start,
  min_clip_len = 25L,
  min_alu_score = 0.60,
  verbose = FALSE
) {
  n_reads <- length(reads)
  if (n_reads == 0L) {
    return(data.frame())
  }
  if (verbose) {
    message("[ALU] Checking ", n_reads, " reads for ALU soft-clips")
  }

  all_qnames <- .safe_qnames(reads)
  all_mapqs <- S4Vectors::mcols(reads)$mapq
  all_mapqs[is.na(all_mapqs)] <- 0L
  all_flags <- S4Vectors::mcols(reads)$flag
  is_reverse <- bitwAnd(all_flags, 0x10L) != 0L
  read_starts <- BiocGenerics::start(reads)
  read_ends <- BiocGenerics::end(reads)
  read_cigars <- GenomicAlignments::cigar(reads)
  read_seqs <- as.character(S4Vectors::mcols(reads)$seq)

  has_clip <- .has_softclip(read_cigars)
  idxs <- which(has_clip)

  if (verbose) {
    message("[ALU] ", length(idxs), " soft-clipped reads to scan")
  }

  rows <- vector("list", 2L * length(idxs))
  k <- 0L

  for (i in idxs) {
    clips <- .get_softclips(read_cigars[i], read_seqs[i])

    for (side in c("lead", "trail")) {
      clip_seq <- clips[[if (side == "lead") "lead" else "trail"]]
      if (is.na(clip_seq) || nchar(clip_seq) < min_clip_len) {
        next
      }

      hit <- .match_clip_to_alu(
        clip_seq,
        alu_seqs,
        min_clip_len = min_clip_len,
        min_score = min_alu_score
      )
      if (is.null(hit)) {
        next
      }

      poly_a <- .detect_poly_a(clip_seq)

      g_pos <- .alu_clip_breakpoint(
        side = side,
        read_start = read_starts[i],
        read_end = read_ends[i]
      )

      # [FIX Finding 2] Exclude out-of-range breakpoints (e.g. leading clip
      # at chromosome position 1 producing NA from .alu_clip_breakpoint()).
      # This prevents coordinate 0 or NA from entering clustering and
      # annotation.
      if (is.na(g_pos) || g_pos < 1L) {
        next
      }

      l_pos <- g_pos - genomic_start + 1L

      # [FIX] Also guard the local coordinate: if it falls outside the
      # reference window, the candidate cannot be annotated.
      if (l_pos < 1L) {
        next
      }

      # [ITEM 2] Extract the non-ALU boundary bases from the clip using the
      # stored pattern alignment positions. For a lead clip these are the
      # bases before pattern_start; for a trail clip, the bases after
      # pattern_end. They often contain the TSD copy carried on the read.
      boundary_bases <- .extract_alu_boundary_bases(
        clip_seq = clip_seq,
        pattern_start = hit$pattern_start,
        pattern_end = hit$pattern_end,
        clip_side = side
      )

      k <- k + 1L
      rows[[k]] <- data.frame(
        read_name = all_qnames[i],
        genomic_pos = g_pos,
        local_pos = l_pos,
        clip_side = side,
        clip_seq = clip_seq,
        alu_subtype = hit$subtype,
        alu_score = hit$score,
        alu_strand = hit$strand,
        aln_start = hit$aln_start,
        aln_end = hit$aln_end,
        pattern_start = hit$pattern_start,
        pattern_end = hit$pattern_end,
        aln_bases = hit$aligned_bases,
        clip_fraction = hit$clip_fraction,
        boundary_bases = boundary_bases,
        poly_a_len = poly_a,
        mapq = all_mapqs[i],
        is_reverse = is_reverse[i],
        stringsAsFactors = FALSE
      )
    }
  }

  if (k == 0L) {
    return(data.frame())
  }
  do.call(rbind, rows[seq_len(k)])
}


#' Summarise one ALU insertion cluster into a result row
#'
#' @param cluster_df  Subset of the candidate data.frame for one cluster.
#' @param ref_dna     Reference sequence for the region (character).
#' @param gene_config Resolved gene config list.
#' @param wt_info     Wildtype read info list (from .prepare_wildtype_info).
#' @param alu_seqs    Named DNAStringSet (for ALU consensus length lookup).
#' @param min_support Minimum supporting reads to report.
#' @keywords internal
.summarise_alu_cluster <- function(
  cluster_df,
  ref_dna,
  gene_config,
  wt_info,
  alu_seqs,
  min_support = 3L,
  min_tsd = 0L
) {
  cluster_df <- cluster_df[
    order(cluster_df$alu_score, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  # Preserve both clip sides for TSD reconstruction before collapsing QNAMEs
  # for independent support counts and per-read summary metrics.
  cluster_all <- cluster_df
  cluster_df <- cluster_df[!duplicated(cluster_df$read_name), , drop = FALSE]
  if (nrow(cluster_df) < min_support) {
    return(NULL)
  }

  best_idx <- which.max(cluster_df$alu_score)
  best_hit <- cluster_df[best_idx, ]

  # [FIX] Use round() before as.integer() to avoid the systematic downward
  # truncation bias of as.integer() on even-length vectors (e.g. median of
  # 100 and 101 is 100.5, which as.integer() truncates to 100).
  median_pos <- as.integer(round(stats::median(cluster_df$genomic_pos)))
  local_pos <- median_pos - gene_config$genomic_start + 1L

  n_fwd <- sum(!cluster_df$is_reverse)
  n_rev <- sum(cluster_df$is_reverse)

  # [FIX Finding 1] TSD is now detected from read evidence (trail vs lead
  # breakpoint gap) rather than reference-flank comparison. The cluster_df
  # and genomic_start are passed so the function can use both sides' reads.
  tsd <- .detect_alu_tsd(
    ref_dna,
    local_pos,
    cluster_df = cluster_all,
    genomic_start = gene_config$genomic_start,
    min_tsd = if (isTRUE(min_tsd > 0L)) as.integer(min_tsd) else 5L
  )

  # [ITEM 9] max(..., na.rm = TRUE) on an all-NA vector is -Inf with a warning.
  poly_a_vals <- cluster_df$poly_a_len
  poly_a_len <- if (length(poly_a_vals) == 0L || all(is.na(poly_a_vals))) {
    NA_integer_
  } else {
    as.integer(max(poly_a_vals, na.rm = TRUE))
  }

  consensus_len <- tryCatch(
    Biostrings::nchar(alu_seqs[[best_hit$alu_subtype]]),
    error = function(e) NA_integer_
  )
  trunc_5p <- .alu_5p_truncation(best_hit$aln_start, consensus_len)

  wt_at_site <- .count_wildtype_at_position(
    wt_info,
    median_pos,
    support_qnames = cluster_df$read_name
  )

  n_support <- nrow(cluster_df)
  depth <- n_support + wt_at_site
  vaf <- if (depth > 0L) n_support / depth else NA_real_

  # Sequence, subtype, strand and truncation must all refer to the same hit.
  best_clip <- as.character(best_hit$clip_seq)

  data.frame(
    InsertionSite = median_pos,
    ALUSubtype = best_hit$alu_subtype,
    ALUOrientation = ifelse(best_hit$alu_strand == "+", "sense", "antisense"),
    TSD_Length = tsd$tsd_len,
    TSD_Sequence = tsd$tsd_seq,
    TSD_Status = tsd$tsd_status,
    PolyALength = poly_a_len,
    ALU5pTruncation = trunc_5p,
    ALU_Sequence = best_clip,
    ALU_RawScore = best_hit$alu_score,
    SupportingReads = as.integer(n_support),
    WildtypeReads = wt_at_site,
    DepthAtBreakpoint = as.integer(depth),
    AlleleFrequency = round(vaf, 4L),
    MeanSupportMAPQ = round(mean(cluster_df$mapq, na.rm = TRUE), 1L),
    StrandBias = round(n_rev / max(1L, n_support), 3L),
    stringsAsFactors = FALSE
  )
}


#' Apply ALU-specific output filters to a result row
#'
#' @param row           Single-row data.frame from .summarise_alu_cluster.
#' @param min_support   Minimum supporting reads (default 3).
#' @param min_alu_score Minimum alignment score (default 0.60).
#' @param vaf_threshold Minimum allele frequency (default 0.01).
#' @param min_tsd       Minimum TSD length; 0 disables (default 0).
#' @keywords internal
.filter_alu_call <- function(
  row,
  min_support = 3L,
  min_alu_score = 0.60,
  vaf_threshold = 0.01,
  min_tsd = 0L
) {
  if (is.null(row)) {
    return(FALSE)
  }

  # Explicit min_tsd requires a read-supported TSD; unknown evidence fails.
  tsd_ok <- if (min_tsd == 0L) {
    TRUE
  } else {
    # An explicit minimum is a real filter: unknown/unresolved TSD evidence
    # cannot satisfy a requested minimum and must not silently pass.
    isTRUE(!is.null(row$TSD_Status) && length(row$TSD_Status) == 1L &&
           !is.na(row$TSD_Status) && row$TSD_Status == "read_supported" &&
           !is.na(row$TSD_Length) && row$TSD_Length >= min_tsd)
  }

  isTRUE(row$SupportingReads >= min_support) &&
    isTRUE(row$ALU_RawScore >= min_alu_score) &&
    isTRUE(
      !is.na(row$AlleleFrequency) && row$AlleleFrequency >= vaf_threshold
    ) &&
    tsd_ok
}


#' Write ALU result artefacts (TSV + HTML report)
#'
#' Shared by the empty and non-empty exit paths of detect_alu() so negative
#' runs still leave an audit artefact (the HTML report has an explicit "No
#' variants detected" state).
#' @keywords internal
.write_alu_output <- function(final_df, output_prefix, output_folder, sample_name,
                              gene_name, gene_config, html_report, verbose) {
  if (is.null(output_prefix) || !nzchar(output_prefix)) return(invisible(NULL))
  sample_folder <- .ensure_dir(file.path(output_folder, sample_name))
  base_name <- paste(output_prefix, sample_name, gene_name, .ts_stamp(), sep = "_")
  if (ncol(final_df) > 0L) {
    tsv_path <- file.path(sample_folder, paste0(base_name, ".tsv"))
    .write_results_tsv(final_df, tsv_path)
    if (verbose) message("[ALU] Results written to: ", tsv_path)
  }
  if (html_report) {
    report_path <- file.path(sample_folder, paste0(base_name, "_report.html"))
    tryCatch({
      talos_html_report(
        result_df = final_df,
        gene_configs = stats::setNames(list(gene_config), gene_name),
        output_file = report_path,
        title = sprintf("TALOS ALU Report - %s | %s", gene_name, sample_name),
        mode = "alu"
      )
      if (verbose) message("[ALU] HTML report written to: ", report_path)
    }, error = function(e) warning("[ALU] HTML report skipped: ", conditionMessage(e)))
  }
  invisible(NULL)
}

#' Detect ALU/SINE mobile-element insertions from a BAM file
#'
#' Analyses soft-clipped reads within the targeted genomic window and
#' identifies candidate ALU insertions by aligning clip sequences to
#' bundled ALU consensus sequences (AluSx, AluY, AluJb families).
#'
#' @param bam_path       Path to an indexed BAM file.
#' @param gene_config    Resolved gene config list (from \code{\link{get_gene_config}}).
#' @param consensus_fa   Path to ALU consensus FASTA (default: bundled file).
#' @param min_support    Minimum supporting reads per cluster (default 3).
#' @param min_alu_score  Minimum normalised alignment score (default 0.60).
#' @param min_clip_len   Minimum soft-clip length to attempt alignment (default 25).
#' @param min_mapq       Minimum MAPQ (default 20).
#' @param cluster_tolerance  Clustering radius in bp (default 15).
#' @param vaf_threshold  Minimum allele frequency to report (default 0.01).
#' @param min_tsd        Minimum TSD length; 0 disables (default 0).
#' @param max_reads_in_region Safety cap on reads loaded (default 200000).
#' @param do_annotate_hotspots Run hotspot annotation (default TRUE).
#' @param hotspot_db_path Optional path to custom hotspot CSV.
#' @param output_prefix  Base prefix for output files (default "TALOS_ALU").
#' @param output_folder  Output directory.
#' @param sample_name    Sample identifier; derived from BAM if NULL.
#' @param html_report    Write a self-contained HTML report (default TRUE).
#' @param verbose        Print progress messages (default TRUE).
#' @return data.frame with one row per ALU insertion detected.
#' @export
detect_alu <- function(
  bam_path,
  gene_config,
  consensus_fa = NULL,
  min_support = 3L,
  min_alu_score = 0.60,
  min_clip_len = 25L,
  min_mapq = 20L,
  cluster_tolerance = 15L,
  vaf_threshold = 0.01,
  min_tsd = 0L,
  max_reads_in_region = 200000L,
  do_annotate_hotspots = TRUE,
  hotspot_db_path = NULL,
  output_prefix = "TALOS_ALU",
  output_folder = "./results",
  sample_name = NULL,
  html_report = TRUE,
  verbose = TRUE
) {
  start_time <- Sys.time()
  if (length(min_tsd) != 1L || is.na(min_tsd) || min_tsd < 0L) {
    stop("min_tsd must be one non-negative integer (0 disables the filter).", call. = FALSE)
  }
  if (length(min_mapq) != 1L || is.na(min_mapq) || min_mapq < 0L) {
    stop("min_mapq must be one non-negative integer.", call. = FALSE)
  }
  gene_name <- gene_config$gene %||% "UNKNOWN"
  gene_config$min_mapq <- as.integer(min_mapq)
  if (is.null(sample_name)) {
    sample_name <- .sample_from_bam(bam_path)
  }

  if (verbose) {
    message(sprintf(
      "[ALU] Starting detection | Gene: %s | Sample: %s",
      gene_name,
      sample_name
    ))
  }

  if (min_tsd > 0L && verbose) {
    message("[ALU] min_tsd is active: calls without a read-supported TSD >= ",
            min_tsd, " bp (including unknown TSD status) will be rejected.")
  }

  alu_seqs <- tryCatch(
    .load_alu_consensus(consensus_fa),
    error = function(e) {
      stop("[ALU] Could not load ALU consensus: ", conditionMessage(e))
    }
  )
  if (verbose) {
    message(sprintf(
      "[ALU] Loaded %d ALU consensus sequences: %s",
      length(alu_seqs),
      paste(names(alu_seqs), collapse = ", ")
    ))
  }

  bam_data <- .load_bam_data_streaming(
    bam_path,
    gene_config,
    compute_pairs = FALSE,
    max_reads = max_reads_in_region,
    verbose = verbose
  )
  all_reads <- bam_data$reads

  if (length(all_reads) == 0L) {
    if (verbose) {
      message("[ALU] No reads in target region.")
    }
    .write_alu_output(data.frame(), output_prefix, output_folder, sample_name,
                      gene_name, gene_config, html_report, verbose)
    return(invisible(data.frame()))
  }

  mapqs <- S4Vectors::mcols(all_reads)$mapq
  mapqs[is.na(mapqs)] <- 0L
  all_reads_for_wt <- all_reads[mapqs >= min_mapq]

  candidates <- .extract_alu_candidates(
    reads = all_reads_for_wt,
    alu_seqs = alu_seqs,
    genomic_start = gene_config$genomic_start,
    min_clip_len = min_clip_len,
    min_alu_score = min_alu_score,
    verbose = verbose
  )

  if (nrow(candidates) == 0L) {
    if (verbose) {
      message("[ALU] No ALU-positive soft-clips found.")
    }
    .write_alu_output(data.frame(), output_prefix, output_folder, sample_name,
                      gene_name, gene_config, html_report, verbose)
    return(invisible(data.frame()))
  }
  if (verbose) {
    message(sprintf(
      "[ALU] %d ALU-positive clip(s) from %d unique reads.",
      nrow(candidates),
      length(unique(candidates$read_name))
    ))
  }

  wt_info <- .prepare_wildtype_info(
    all_reads_for_wt,
    gene_config$genomic_start,
    gene_config$genomic_end
  )
  wt_info$sample_name <- sample_name

  clusters <- .cluster_breakpoints(candidates$genomic_pos, cluster_tolerance)
  ref_dna <- gene_config$genomic_ref_seq
  if (is.null(ref_dna) || is.na(ref_dna)) {
    stop("[ALU] No genomic reference sequence in gene_config.")
  }

  results <- list()
  for (cl in clusters) {
    cl_df <- candidates[candidates$genomic_pos %in% cl, , drop = FALSE]
    row <- .summarise_alu_cluster(
      cluster_df = cl_df,
      ref_dna = ref_dna,
      gene_config = gene_config,
      wt_info = wt_info,
      alu_seqs = alu_seqs,
      min_support = min_support,
      min_tsd = min_tsd
    )
    if (
      .filter_alu_call(row, min_support, min_alu_score, vaf_threshold, min_tsd)
    ) {
      results[[length(results) + 1L]] <- row
    }
  }

  if (length(results) == 0L) {
    if (verbose) {
      message("[ALU] No insertions passed all filters.")
    }
    .write_alu_output(data.frame(), output_prefix, output_folder, sample_name,
                      gene_name, gene_config, html_report, verbose)
    return(invisible(data.frame()))
  }

  final_df <- do.call(rbind, lapply(results, as.data.frame))
  final_df$Sample <- sample_name
  final_df$Gene <- gene_name
  final_df$Genome <- gene_config$build %||% "unknown"

  col_order <- c(
    "Sample",
    "Gene",
    "Genome",
    setdiff(names(final_df), c("Sample", "Gene", "Genome"))
  )
  final_df <- final_df[, col_order, drop = FALSE]

  if (do_annotate_hotspots) {
    final_df <- annotate_hotspots(
      final_df,
      db_path = hotspot_db_path,
      genome_build = gene_config$build,
      pos_col = "InsertionSite"
    )
  } else {
    final_df$Hotspot <- FALSE
    final_df$HotspotName <- NA_character_
  }

  if (!is.null(gene_config$target_exons)) {
    final_df <- .annotate_exonic_region(
      final_df,
      gene_config$target_exons,
      pos_col = "InsertionSite"
    )
  }

  .write_alu_output(final_df, output_prefix, output_folder, sample_name,
                    gene_name, gene_config, html_report, verbose)

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (verbose) {
    message(sprintf(
      "[ALU] Done | %d insertion(s) | %.1f s",
      nrow(final_df),
      elapsed
    ))
  }
  invisible(final_df)
}


#' Simplified entry point for ALU insertion detection
#'
#' Calls \code{\link{get_gene_config}} then \code{\link{detect_alu}}.
#' Gene-specific ALU settings can be provided in \code{gene_config.yaml}
#' under an \code{alu_settings:} block with keys:
#'   \code{min_support}, \code{min_alu_score}, \code{min_clip_len},
#'   \code{min_mapq}, \code{cluster_tolerance}, \code{vaf_threshold}, \code{min_tsd}.
#'
#' @inheritParams detect_alu
#' @param gene       Gene symbol (must be present in \code{yaml_path}).
#' @param build      Genome build: \code{"hg19"} or \code{"hg38"}.
#' @param padding    Padding around targeted exons (default 500).
#' @param yaml_path  Path to gene config YAML (defaults to bundled file).
#' @param bsgenome   Optional BSgenome object or package string.
#' @param exon_padding Number of flanking exons to add to the window.
#' @param ...        Additional arguments forwarded to \code{detect_alu}.
#' @return data.frame of ALU insertions (invisibly).
#' @export
talos_alu <- function(
  bam_path,
  gene,
  build = "hg19",
  padding = 500L,
  consensus_fa = NULL,
  min_support = NULL,
  min_alu_score = NULL,
  min_clip_len = NULL,
  min_mapq = NULL,
  cluster_tolerance = NULL,
  vaf_threshold = NULL,
  min_tsd = NULL,
  sample_name = NULL,
  output_prefix = "TALOS_ALU",
  output_folder = "./results",
  html_report = TRUE,
  do_annotate_hotspots = TRUE,
  hotspot_db_path = NULL,
  yaml_path = system.file("extdata", "gene_config.yaml", package = "TALOS"),
  bsgenome = NULL,
  exon_padding = 0L,
  verbose = TRUE,
  ...
) {
  # Match talos(): explicit arguments override gene_settings, which override defaults.
  yaml_gs <- tryCatch(yaml::read_yaml(yaml_path)[[gene]]$gene_settings,
                      error = function(e) NULL)
  if (missing(padding) && !is.null(yaml_gs$padding)) padding <- as.integer(yaml_gs$padding)
  if (missing(exon_padding) && !is.null(yaml_gs$exon_padding)) exon_padding <- as.integer(yaml_gs$exon_padding)

  config <- get_gene_config(
    gene = gene,
    build = build,
    config_path = yaml_path,
    bsgenome = bsgenome,
    padding = padding,
    exon_padding = exon_padding
  )
  config$build <- build

  defaults <- formals(detect_alu)[c("min_support", "min_alu_score", "min_clip_len",
                                     "min_mapq", "cluster_tolerance", "vaf_threshold", "min_tsd")]
  p <- .resolve_settings(environment(), defaults, config$alu_settings)

  detect_alu(
    bam_path = bam_path,
    gene_config = config,
    consensus_fa = consensus_fa,
    min_support = p$min_support,
    min_alu_score = p$min_alu_score,
    min_clip_len = p$min_clip_len,
    min_mapq = p$min_mapq,
    cluster_tolerance = p$cluster_tolerance,
    vaf_threshold = p$vaf_threshold,
    min_tsd = p$min_tsd,
    do_annotate_hotspots = do_annotate_hotspots,
    hotspot_db_path = hotspot_db_path,
    output_prefix = output_prefix,
    output_folder = output_folder,
    sample_name = sample_name,
    html_report = html_report,
    verbose = verbose,
    ...
  )
}

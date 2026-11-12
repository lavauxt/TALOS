#' Resolve the BAM sequence name matching a configured chromosome
#'
#' Gene configs use UCSC names ("chr11"); BAMs may use Ensembl names ("11",
#' "MT"). Querying a BAM with a name absent from its header is an error, so
#' the name is resolved against the header first.
#'
#' Replaces the post-hoc mapSeqlevels() conversion (ITEM 10): that branch
#' could only run after a successful query, i.e. when the name already
#' matched. The collision guard of ITEM 10 is kept implicitly: only the one
#' target contig is ever renamed.
#'
#' @param bam_path Path to an indexed BAM.
#' @param chrom Configured chromosome name.
#' @return The BAM's name for \code{chrom}; \code{chrom} itself (with a
#'   warning) when no candidate is found.
#' @keywords internal
.resolve_bam_seqname <- function(bam_path, chrom) {
  bam_names <- tryCatch(
    GenomeInfoDb::seqnames(Rsamtools::seqinfo(Rsamtools::BamFile(bam_path))),
    error = function(e) character(0)
  )
  if (!length(bam_names) || chrom %in% bam_names) return(chrom)
  bare <- sub("^chr", "", chrom)
  candidates <- unique(c(
    bare, paste0("chr", bare),
    if (bare %in% c("M", "MT")) c("MT", "M", "chrM", "chrMT")
  ))
  hit <- candidates[candidates %in% bam_names]
  if (length(hit)) return(hit[[1L]])
  warning("Chromosome '", chrom, "' not found in BAM header of ", basename(bam_path),
          ". Coverage will be NA.")
  chrom
}

# Rename the BAM contig back to the configured name so every downstream
# GRanges (exons, coverage lookups, overlaps) uses one naming style.
.rename_to_config <- function(x, bam_chrom, chrom) {
  if (is.null(x) || identical(bam_chrom, chrom)) return(x)
  lv <- GenomeInfoDb::seqlevels(x)
  if (!bam_chrom %in% lv || chrom %in% lv) return(x)
  lv[lv == bam_chrom] <- chrom
  GenomeInfoDb::renameSeqlevels(x, stats::setNames(lv, GenomeInfoDb::seqlevels(x)))
}

# Chunked reader shared by single reads and read pairs (DRY).
.read_bam_chunks <- function(bam_path, reader, max_n, chunk_size, label, verbose) {
  bam_file <- Rsamtools::BamFile(bam_path, yieldSize = chunk_size)
  open(bam_file)
  on.exit(close(bam_file), add = TRUE)
  chunks <- list()
  total <- 0L
  repeat {
    chunk <- reader(bam_file)
    if (length(chunk) == 0L) break
    total <- total + length(chunk)
    if (!is.null(max_n) && total >= max_n) {
      chunk <- chunk[seq_len(max_n - (total - length(chunk)))]
      chunks[[length(chunks) + 1L]] <- chunk
      warning("[TALOS] max_reads_in_region limit reached (", max_n,
              "); this analysis may be truncated. Increase the limit or set it to NULL for an uncapped run.",
              call. = FALSE)
      if (verbose) message("[TALOS] Stopping early at max_reads = ", max_n)
      break
    }
    chunks[[length(chunks) + 1L]] <- chunk
    if (verbose) {
      message("[TALOS] Loaded ", label, " chunk ", length(chunks), " (",
              length(chunk), ", total ", total, ")")
    }
  }
  if (length(chunks) == 0L) NULL else do.call(c, chunks)
}

.load_bam_data_streaming <- function(
    bam_path,
    gene_config,
    compute_pairs = FALSE,
    max_reads = NULL,
    chunk_size = 50000L,
    verbose = FALSE
) {
  target_chrom <- gene_config$chrom
  bam_chrom <- .resolve_bam_seqname(bam_path, target_chrom)
  if (verbose) {
    message("[TALOS] Loading BAM region: ", bam_chrom, ":",
            gene_config$genomic_start, "-", gene_config$genomic_end)
  }
  which_range <- GenomicRanges::GRanges(
    seqnames = bam_chrom,
    ranges = IRanges::IRanges(start = gene_config$genomic_start, end = gene_config$genomic_end)
  )
  # "qual" is loaded since v0.4.20 so per-read base-quality statistics can be
  # computed and the base-quality filter actually applies to candidate rows
  # (previously ITDBaseQualityStatus was hard-coded to a "not loaded" state).
  param <- Rsamtools::ScanBamParam(
    which = which_range,
    what = c("seq", "qname", "cigar", "flag", "mapq", "isize", "qual"),
    flag = Rsamtools::scanBamFlag(
      isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE,
      isDuplicate = FALSE, isNotPassingQualityControls = FALSE
    )
  )

  all_reads <- .read_bam_chunks(
    bam_path, function(bf) GenomicAlignments::readGAlignments(bf, param = param, use.names = TRUE),
    max_reads, chunk_size, "read", verbose
  )
  if (verbose) message("[TALOS] Total reads loaded: ", length(all_reads %||% integer(0L)))
  if (is.null(all_reads) || length(all_reads) == 0L) {
    if (verbose) message("[TALOS] No reads found in target region.")
    return(list(reads = GenomicAlignments::GAlignments(), cov = NULL, pairs = NULL, cov_list = NULL))
  }
  all_reads <- .rename_to_config(all_reads, bam_chrom, target_chrom)

  if (verbose) message("[TALOS] Computing coverage...")
  cov_rle_list <- GenomicAlignments::coverage(all_reads)
  cov <- if (target_chrom %in% names(cov_rle_list)) cov_rle_list[[target_chrom]] else NULL

  pairs <- NULL
  if (compute_pairs) {
    if (verbose) message("[TALOS] Loading read pairs...")
    pairs <- .read_bam_chunks(
      bam_path,
      function(bf) suppressWarnings(GenomicAlignments::readGAlignmentPairs(bf, param = param, use.names = TRUE)),
      max_reads, chunk_size, "pair", verbose
    )
    pairs <- .rename_to_config(pairs, bam_chrom, target_chrom)
    if (verbose) message("[TALOS] Total pairs loaded: ", length(pairs %||% integer(0L)))
  }
  list(reads = all_reads, cov = cov, pairs = pairs, cov_list = cov_rle_list)
}

.filter_reads_by_cigar <- function(
  reads,
  min_mapq,
  min_ins_filter,
  verbose = FALSE
) {
  if (verbose) {
    message(
      "[TALOS] Pre-filtering reads: min_mapq=",
      min_mapq,
      ", min_ins_filter=",
      min_ins_filter
    )
  }

  flags <- S4Vectors::mcols(reads)$flag
  mapqs <- S4Vectors::mcols(reads)$mapq
  mapqs[is.na(mapqs)] <- 0L
  cigars <- GenomicAlignments::cigar(reads)
  unmapped <- bitwAnd(flags, 0x4) != 0L

  # [FIX] Unmapped reads must be excluded: they have no valid genomic
  # coordinates and would produce NA starts/ends in downstream candidate
  # extraction. The previous code OR'd `unmapped` into the keep condition,
  # allowing unmapped reads to pass when min_mapq == 0. This is a latent
  # bug: unmapped reads contribute NA breakpoints, corrupt coverage, and
  # can crash pairwiseAlignment when their sequence is NA.
  op_table <- if (.has_cigarillo) {
    cigarillo::tabulate_cigar_ops(cigars)
  } else {
    suppressWarnings(GenomicAlignments::cigarOpTable(cigars))
  }

  has_softclip <- .has_softclip(cigars)

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

  # [FIX] Removed `unmapped |` from the OR: unmapped reads must NOT pass
  # the filter. Only mapped reads with soft-clips or sufficient net insertion
  # (and adequate MAPQ) are kept.
  keep <- (!unmapped & (has_softclip | (net_ins >= min_ins_filter))) &
    (mapqs >= min_mapq)
  keep[is.na(keep)] <- FALSE

  filtered <- reads[keep]
  if (verbose) {
    message(
      "[TALOS] ",
      sum(keep),
      " reads passed pre-filter (",
      length(reads),
      " total)"
    )
  }
  filtered
}

#' Build Optimized ScanBamParam for Evidence Extraction
#' @param gene_config Resolved gene config list from get_gene_config()
#' @param buffer_bp Flanking window size in base pairs
#' @param bam_chrom Chromosome name as spelled in the BAM header.
#' @return Rsamtools::ScanBamParam object
.build_optimized_bam_param <- function(gene_config, buffer_bp = 15000L,
                                       bam_chrom = gene_config$chrom) {
  target_region <- suppressWarnings(GenomicRanges::trim(GenomicRanges::GRanges(
    seqnames = bam_chrom,
    ranges = IRanges::IRanges(
      start = max(1L, gene_config$genomic_start - buffer_bp),
      end = gene_config$genomic_end + buffer_bp
    ),
    strand = "*"
  )))
  bam_flags <- Rsamtools::scanBamFlag(
    isUnmappedQuery = FALSE,
    isSecondaryAlignment = FALSE,
    isSupplementaryAlignment = NA,
    isDuplicate = FALSE,
    isNotPassingQualityControls = FALSE
  )
  Rsamtools::ScanBamParam(
    flag = bam_flags,
    what = c("qname", "flag", "rname", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize"),
    tag = c("SA", "NM"),
    which = target_region
  )
}

#' Mean Phred base quality per read.
#'
#' @param quals Character vector of Phred+33 ASCII quality strings (as
#'   obtained from \code{as.character()} on the BAM \code{qual} column).
#' @return Numeric vector of mean Phred scores; \code{NA_real_} for reads
#'   with a missing or empty quality string.
#' @keywords internal
.mean_base_quality <- function(quals) {
  qchar <- if (is.null(quals)) character(0) else as.character(quals)
  out <- rep(NA_real_, length(qchar))
  ok <- !is.na(qchar) & nzchar(qchar)
  if (any(ok)) {
    # Phred+33: each ASCII byte encodes Q = code - 33. utf8ToInt() is safe
    # here because the quality alphabet is plain ASCII.
    out[ok] <- vapply(
      qchar[ok],
      function(s) mean(as.integer(utf8ToInt(s) - 33L)),
      numeric(1L)
    )
  }
  out
}

#' Does the loaded read set carry usable BAM QUAL strings?
#' @keywords internal
.bam_has_base_quality <- function(reads) {
  mc <- S4Vectors::mcols(reads)
  if (is.null(mc) || !"qual" %in% names(mc) || length(mc$qual) == 0L) return(FALSE)
  q1 <- as.character(mc$qual[1L])
  length(q1) == 1L && !is.na(q1) && nzchar(q1)
}

#' Drop reads whose mean base quality is below a threshold
#'
#' Reads without a quality string are retained (they cannot be evaluated);
#' the run-level status fields report that base quality was unavailable for
#' them rather than silently assuming they passed.
#'
#' @param reads GAlignments loaded with a \code{qual} mcols column.
#' @param min_base_quality Minimum mean Phred score per read.
#' @return Filtered \code{reads}.
#' @keywords internal
.filter_reads_by_base_quality <- function(reads, min_base_quality, verbose = FALSE) {
  if (length(reads) == 0L) return(reads)
  quals <- S4Vectors::mcols(reads)$qual
  if (is.null(quals) || length(quals) != length(reads)) {
    if (verbose) message("[TALOS] Base-quality filter skipped: no QUAL column loaded.")
    return(reads)
  }
  mean_bq <- .mean_base_quality(quals)
  keep <- is.na(mean_bq) | mean_bq >= min_base_quality
  keep[is.na(keep)] <- TRUE
  if (verbose) {
    message(
      "[TALOS] Base-quality filter (min mean Phred = ", min_base_quality, "): ",
      sum(keep), " / ", length(reads), " reads retained"
    )
  }
  reads[keep]
}

#' Flag sequences carrying adapter read-through
#'
#' A read is flagged when any configured adapter prefix of at least
#' \code{min_overlap} bases occurs contiguously in the read (exact match).
#' Adapter read-through reaches the 3' end of the stored sequence, so this
#' covers both soft-clipped adapter segments and adapter bases that partially
#' aligned to the reference.
#'
#' @param seqs Character vector of read sequences.
#' @param adapters Character vector of adapter sequences.
#' @param min_overlap Minimum contiguous adapter bases to flag a read.
#' @return Logical vector, one flag per input sequence.
#' @keywords internal
.screen_adapter_seqs <- function(seqs, adapters, min_overlap = 10L) {
  seqs <- if (is.null(seqs)) character(0) else as.character(seqs)
  seqs[is.na(seqs)] <- ""
  hit <- rep(FALSE, length(seqs))
  adapters <- adapters[!is.na(adapters) & nzchar(adapters)]
  for (ad in toupper(adapters)) {
    # Prefixes longer than min_overlap are only tested where needed: a hit at
    # length p implies a hit at p - 1, so the OR accumulates the widest match.
    for (p in seq.int(min_overlap, nchar(ad))) {
      hit <- hit | grepl(substr(ad, 1L, p), seqs, fixed = TRUE)
    }
  }
  hit
}

#' Screen loaded reads for adapter contamination (report-only)
#'
#' Screening does not remove reads; it produces the region-wide counts and
#' the affected QNAME set reported per variant via
#' \code{ITDAdapterPositiveSupportReads} and
#' \code{ITDAdapterScreeningStatus}.
#'
#' @param reads GAlignments with a \code{seq} mcols column.
#' @param adapter_sequences Character vector of adapter sequences.
#' @param min_overlap Minimum contiguous adapter bases to flag a read.
#' @return list(screened, positive, qnames).
#' @keywords internal
.screen_adapter_reads <- function(reads, adapter_sequences, min_overlap = 10L,
                                  verbose = FALSE) {
  adapters <- adapter_sequences %||% character(0)
  adapters <- adapters[!is.na(adapters) & nzchar(adapters)]
  if (length(reads) == 0L || length(adapters) == 0L) {
    return(list(screened = length(reads), positive = 0L, qnames = character(0)))
  }
  seqs <- toupper(as.character(S4Vectors::mcols(reads)$seq))
  hit <- .screen_adapter_seqs(seqs, adapters, min_overlap)
  qnames <- unique(.safe_qnames(reads)[hit])
  if (verbose) {
    message(
      "[TALOS] Adapter screening: ", length(hit), " / ", length(reads),
      " reads carry >= ", min_overlap, " bases of ", length(adapters),
      " configured adapter(s) (", length(qnames), " unique QNAMEs)"
    )
  }
  list(screened = length(reads), positive = as.integer(sum(hit)), qnames = qnames)
}

#' Extract BAM Evidence into a Plain data.frame
#' @param bam_file Path to BAM file
#' @param gene_config Resolved gene config list
#' @param buffer_bp Flanking window size in base pairs (default 15000)
#' @param min_mapq Minimum MAPQ of a record (default 20)
#' @return data.frame of parsed BAM records (one row per alignment record).
#'   \code{rname}/\code{mrnm} keep the BAM spelling; only coordinates are
#'   used downstream.
.get_evidence_dt <- function(bam_file, gene_config, buffer_bp = 15000L, min_mapq = 20L) {
  bam_chrom <- .resolve_bam_seqname(bam_file, gene_config$chrom)
  param <- .build_optimized_bam_param(gene_config, buffer_bp = buffer_bp, bam_chrom = bam_chrom)
  res <- Rsamtools::scanBam(bam_file, param = param)[[1L]]
  if (length(res$qname) == 0L) return(data.frame())
  n <- length(res$qname)
  sa_col <- if (!is.null(res$tag) && !is.null(res$tag$SA)) res$tag$SA else rep(NA_character_, n)
  out <- data.frame(
    qname = res$qname, flag = res$flag, rname = as.character(res$rname),
    pos = res$pos, qwidth = res$qwidth, mapq = res$mapq, cigar = res$cigar,
    mrnm = as.character(res$mrnm), mpos = res$mpos, isize = res$isize,
    sa = sa_col, stringsAsFactors = FALSE
  )
  out[!is.na(out$mapq) & out$mapq >= min_mapq, , drop = FALSE]
}

# Empty evidence table (one schema for every evidence source).
.empty_dup_evidence <- function() {
  data.frame(qname = character(), start = integer(), end = integer(), prob = numeric(),
             evidence_type = character(), precise = logical(), stringsAsFactors = FALSE)
}

#' Tandem-duplication evidence from split reads (SA tags)
#'
#' A read crossing the novel junction of a tandem duplication [S, E] has one
#' alignment ending at E with its trailing part clipped, and one starting at S
#' with its leading part clipped, with S < E (a deletion gives the mirror
#' image: the left-clipped piece lies downstream). The previous version kept
#' only \code{pmax/pmin} of the two positions, which discarded this direction:
#' on a + strand gene every split read (deletions, inversions, fusions to
#' another chromosome) became a PTD back-edge, and on a - strand gene none
#' did. Here:
#' \itemize{
#'   \item SA entries on another chromosome or strand are ignored
#'     (translocation / inversion, not a tandem duplication);
#'   \item the clip side of each piece decides which end is S and which is E;
#'   \item the primary and supplementary records of a read describe the same
#'     junction once (deduplicated by read name).
#' }
#' @param dt data.frame of BAM records (from .get_evidence_dt)
#' @return data.frame(qname, start = S, end = E, prob, evidence_type, precise)
.parse_sa_tags <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0L) return(.empty_dup_evidence())
  split_dt <- dt[!is.na(dt$sa) & dt$sa != "", , drop = FALSE]
  if (nrow(split_dt) == 0L) return(.empty_dup_evidence())

  sa_list <- strsplit(split_dt$sa, ";", fixed = TRUE)
  rep_idx <- rep(seq_len(nrow(split_dt)), lengths(sa_list))
  sa_flat <- unlist(sa_list, use.names = FALSE)
  keep <- nzchar(sa_flat)
  rep_idx <- rep_idx[keep]
  tokens <- strsplit(sa_flat[keep], ",", fixed = TRUE)
  valid <- lengths(tokens) >= 5L
  rep_idx <- rep_idx[valid]
  tokens <- tokens[valid]
  if (length(tokens) == 0L) return(.empty_dup_evidence())

  tok <- function(i) vapply(tokens, `[`, character(1L), i)
  sa_rname <- tok(1L)
  sa_pos <- suppressWarnings(as.integer(tok(2L)))
  sa_strand <- tok(3L)
  sa_cigar <- tok(4L)
  sa_mapq <- suppressWarnings(as.numeric(tok(5L)))

  rec <- split_dt[rep_idx, , drop = FALSE]
  rec_strand <- ifelse(bitwAnd(rec$flag, 0x10L) != 0L, "-", "+")
  ok <- !is.na(sa_pos) & !is.na(sa_mapq) & sa_rname == rec$rname & sa_strand == rec_strand
  if (!any(ok)) return(.empty_dup_evidence())
  rec <- rec[ok, , drop = FALSE]
  sa_pos <- sa_pos[ok]; sa_cigar <- sa_cigar[ok]; sa_mapq <- sa_mapq[ok]

  a <- .cigar_layout(rec$cigar)
  b <- .cigar_layout(sa_cigar)
  a_start <- rec$pos; a_end <- rec$pos + a$width - 1L
  b_start <- sa_pos;  b_end <- sa_pos + b$width - 1L
  a_right <- a$trail > a$lead; a_left <- a$lead > a$trail
  b_right <- b$trail > b$lead; b_left <- b$lead > b$trail

  S <- ifelse(a_right & b_left, b_start, ifelse(a_left & b_right, a_start, NA_integer_))
  E <- ifelse(a_right & b_left, a_end, ifelse(a_left & b_right, b_end, NA_integer_))
  is_dup <- !is.na(S) & !is.na(E) & S < E
  if (!any(is_dup)) return(.empty_dup_evidence())

  mean_mapq <- (rec$mapq + sa_mapq) / 2
  ev <- data.frame(
    qname = rec$qname[is_dup],
    start = as.integer(S[is_dup]),
    end = as.integer(E[is_dup]),
    prob = pmin(0.99, pmax(0.10, (mean_mapq[is_dup] / 60) * 0.95)),
    evidence_type = "split_read",
    precise = TRUE,
    stringsAsFactors = FALSE
  )
  ev <- ev[order(ev$qname, -ev$prob), , drop = FALSE]
  ev[!duplicated(ev[, c("qname", "start", "end")]), , drop = FALSE]
}

#' Tandem-duplication evidence from discordant read pairs
#'
#' Pairs spanning a tandem-duplication junction are everted (RF): the
#' leftmost mate maps on the reverse strand, the rightmost on the forward
#' strand, far apart. The previous test \code{abs(isize) > max | isize < 0}
#' selected the reverse mate of every normal pair (isize is negative for the
#' rightmost mate by definition) and ignored orientation and chromosome.
#'
#' Shadow reads (mate unmapped) carry no partner coordinate, so they no
#' longer create graph edges (the approximated mate position of BUG 5 created
#' spurious short edges instead).
#'
#' [FIX] The end coordinate previously used \code{d$qwidth} (the current
#' read's width) regardless of which mate is rightmost. When the current
#' read is the leftmost mate, the end should use the mate's width. Since the
#' mate's width is not in the data, we now use \code{pmax(pos, mpos) +
#' max(qwidth) - 1L} as an upper bound (slight over-estimate is safer than
#' under-estimate for exon-overlap mapping). For same-length pairs (the
#' common case) this is exact.
#'
#' [ITEM 5] For mixed-length read pairs the upper bound above is off by the
#' width difference. We now look up the mate's qwidth from the same data
#' frame: each pair typically contributes two rows (one per mate) to \code{dt},
#' so a hash from \code{paste(qname, pos)} to \code{qwidth} lets us recover
#' the mate's actual width by keying on \code{paste(qname, mpos)}. When the
#' mate's record is missing from the scanned region (only one mate in
#' \code{dt}), we fall back to the current read's qwidth, which is exact for
#' same-length pairs and the original approximation otherwise.
#'
#' @param dt data.frame of BAM records
#' @param expected_isize_max Minimum mate separation (bp) to call a pair
#'   discordant (default 1000)
#' @return data.frame(qname, start, end, prob, evidence_type, precise)
.extract_discordant_and_shadow <- function(dt, expected_isize_max = 1000L) {
  if (is.null(dt) || nrow(dt) == 0L) return(.empty_dup_evidence())
  f <- dt$flag
  mate_chr <- ifelse(dt$mrnm %in% c("=", NA), dt$rname, dt$mrnm)
  usable <- bitwAnd(f, 0x1L) != 0L & bitwAnd(f, 0x4L) == 0L & bitwAnd(f, 0x8L) == 0L &
    bitwAnd(f, 0x800L) == 0L & !is.na(dt$mpos) & mate_chr == dt$rname
  if (!any(usable)) return(.empty_dup_evidence())
  d <- dt[usable, , drop = FALSE]
  read_rev <- bitwAnd(d$flag, 0x10L) != 0L
  mate_rev <- bitwAnd(d$flag, 0x20L) != 0L
  read_left <- d$pos <= d$mpos
  left_rev <- ifelse(read_left, read_rev, mate_rev)
  right_rev <- ifelse(read_left, mate_rev, read_rev)
  everted <- left_rev & !right_rev & abs(d$mpos - d$pos) > expected_isize_max
  if (!any(everted)) return(.empty_dup_evidence())
  d <- d[everted, , drop = FALSE]

  # Build a (qname, pos) -> CIGAR reference-span lookup for both mates.
  # Query width is not a genomic span: soft clips inflate it and skipped
  # introns (N operations) are absent from it. If a mate record is missing,
  # use the current read's reference span as the documented approximation.
  # Use CIGAR reference span for each mate. qwidth counts query bases and
  # therefore includes soft clips while excluding introns (N); it is not a
  # genomic span for clipped or spliced alignments.
  ref_width <- .cigar_layout(dt$cigar)$width
  bad_width <- is.na(ref_width) | ref_width <= 0L
  ref_width[bad_width] <- as.integer(dt$qwidth[bad_width])
  pos_key <- paste(dt$qname, dt$pos, sep = "\1")
  width_by_pos <- stats::setNames(as.integer(ref_width), pos_key)
  own_key <- paste(d$qname, d$pos, sep = "\1")
  mate_key <- paste(d$qname, d$mpos, sep = "\1")
  own_width <- as.integer(unname(width_by_pos[own_key]))
  mate_width <- as.integer(unname(width_by_pos[mate_key]))
  own_width[is.na(own_width)] <- as.integer(d$qwidth[is.na(own_width)])
  mate_width[is.na(mate_width)] <- own_width[is.na(mate_width)]
  rightmost_width <- ifelse(d$pos >= d$mpos, own_width, mate_width)
  rightmost_width[is.na(rightmost_width)] <- own_width[is.na(rightmost_width)]

  ev <- data.frame(
    qname = d$qname,
    start = as.integer(pmin(d$pos, d$mpos)),
    end = as.integer(pmax(d$pos, d$mpos) + rightmost_width - 1L),
    prob = pmin(0.85, pmax(0.15, (d$mapq / 60) * 0.80)),
    evidence_type = "discordant_pair",
    precise = FALSE,
    stringsAsFactors = FALSE
  )
  ev <- ev[order(ev$qname, -ev$prob), , drop = FALSE]
  ev[!duplicated(ev$qname), , drop = FALSE]
}

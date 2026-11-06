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
  param <- Rsamtools::ScanBamParam(
    which = which_range,
    what = c("seq", "qname", "cigar", "flag", "mapq", "isize"),
    flag = Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE)
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

  keep <- (unmapped | has_softclip | (net_ins >= min_ins_filter)) &
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

#' Extract BAM Evidence into a Plain data.frame
#' @param bam_file Path to BAM file
#' @param gene_config Resolved gene config list
#' @param buffer_bp Flanking window size in base pairs (default 15000)
#' @param min_mapq Minimum MAPQ of a record (default 10)
#' @return data.frame of parsed BAM records (one row per alignment record).
#'   \code{rname}/\code{mrnm} keep the BAM spelling; only coordinates are
#'   used downstream.
.get_evidence_dt <- function(bam_file, gene_config, buffer_bp = 15000L, min_mapq = 10L) {
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
  ev <- data.frame(
    qname = d$qname,
    start = as.integer(pmin(d$pos, d$mpos)),
    end = as.integer(pmax(d$pos, d$mpos) + d$qwidth - 1L),
    prob = pmin(0.85, pmax(0.15, (d$mapq / 60) * 0.80)),
    evidence_type = "discordant_pair",
    precise = FALSE,
    stringsAsFactors = FALSE
  )
  ev <- ev[order(ev$qname, -ev$prob), , drop = FALSE]
  ev[!duplicated(ev$qname), , drop = FALSE]
}

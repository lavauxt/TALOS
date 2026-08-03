#' Build Optimized ScanBamParam for Evidence Extraction
#' @param gene_config Resolved gene config list from get_gene_config()
#' @param buffer_bp Flanking window size in base pairs
#' @return Rsamtools::ScanBamParam object
.build_optimized_bam_param <- function(gene_config, buffer_bp = 15000L) {
  target_region <- suppressWarnings(GenomicRanges::trim(GenomicRanges::GRanges(
    seqnames = gene_config$chrom,
    ranges   = IRanges::IRanges(
      start = max(1L, gene_config$genomic_start - buffer_bp),
      end   = gene_config$genomic_end + buffer_bp
    ),
    strand = "*"
  )))

  bam_flags <- Rsamtools::scanBamFlag(
    isUnmappedQuery              = FALSE,
    isSecondaryAlignment         = FALSE,
    isSupplementaryAlignment     = NA,  
    isDuplicate                  = FALSE,
    isNotPassingQualityControls  = FALSE
  )

  Rsamtools::ScanBamParam(
    flag  = bam_flags,
    what  = c("qname", "flag", "rname", "pos", "qwidth", "mapq", "cigar", "mrnm", "mpos", "isize"),
    tag   = c("SA", "NM"),
    which = target_region
  )
}

#' Extract BAM Evidence into a Plain data.frame
#' @param bam_file Path to BAM file
#' @param gene_config Resolved gene config list
#' @param buffer_bp Flanking window size in base pairs (default 15000)
#' @return data.frame of parsed BAM records (one row per alignment record)
.get_evidence_dt <- function(bam_file, gene_config, buffer_bp = 15000L) {
  param <- .build_optimized_bam_param(gene_config, buffer_bp = buffer_bp)
  res <- Rsamtools::scanBam(bam_file, param = param)[[1]]

  if (length(res$qname) == 0L) return(data.frame())

  sa_col <- if (!is.null(res$tag) && !is.null(res$tag$SA)) res$tag$SA else NA_character_

  out <- data.frame(
    qname  = res$qname,
    flag   = res$flag,
    rname  = as.character(res$rname),
    pos    = res$pos,
    qwidth = res$qwidth,
    mapq   = res$mapq,
    cigar  = res$cigar,
    mpos   = res$mpos,
    isize  = res$isize,
    sa     = sa_col,
    stringsAsFactors = FALSE
  )


  out[!is.na(out$mapq) & out$mapq >= 10, , drop = FALSE]
}

#' Parse Split-Read Supplementary Alignments (SA Tags)
#'
#' Vectorized (no per-row R loop): expands each read's ';'-delimited SA
#' entries into rows, then parses the fixed 6-field SA format
#' ("rname,pos,strand,CIGAR,mapq,NM") with vapply instead of a nested loop.
#'
#' @param dt data.frame of BAM records (from .get_evidence_dt)
#' @return data.frame of split-read coordinate pairs with raw probabilities
.parse_sa_tags <- function(dt) {
  if (nrow(dt) == 0L) return(data.frame())
  split_dt <- dt[!is.na(dt$sa) & dt$sa != "", , drop = FALSE]
  if (nrow(split_dt) == 0L) return(data.frame())

  sa_list <- strsplit(split_dt$sa, ";", fixed = TRUE)
  n_per   <- lengths(sa_list)
  rep_idx <- rep(seq_len(nrow(split_dt)), n_per)
  sa_flat <- unlist(sa_list, use.names = FALSE)
  keep    <- nzchar(sa_flat)
  rep_idx <- rep_idx[keep]; sa_flat <- sa_flat[keep]
  if (length(sa_flat) == 0L) return(data.frame())

  tokens  <- strsplit(sa_flat, ",", fixed = TRUE)
  n_tok   <- lengths(tokens)
  valid   <- n_tok >= 5L
  rep_idx <- rep_idx[valid]; tokens <- tokens[valid]
  if (length(tokens) == 0L) return(data.frame())

  sa_pos  <- suppressWarnings(as.numeric(vapply(tokens, `[`, character(1L), 2L)))
  sa_mapq <- suppressWarnings(as.numeric(vapply(tokens, `[`, character(1L), 5L)))
  ok      <- !is.na(sa_pos) & !is.na(sa_mapq)
  rep_idx <- rep_idx[ok]; sa_pos <- sa_pos[ok]; sa_mapq <- sa_mapq[ok]
  if (length(rep_idx) == 0L) return(data.frame())

  read_pos  <- split_dt$pos[rep_idx]
  read_mapq <- split_dt$mapq[rep_idx]
  mean_mapq <- (read_mapq + sa_mapq) / 2

  data.frame(
    qname         = split_dt$qname[rep_idx],
    from_coord    = pmax(read_pos, sa_pos),
    to_coord      = pmin(read_pos, sa_pos),
    prob          = pmin(0.99, pmax(0.10, (mean_mapq / 60) * 0.95)),
    evidence_type = "split_read",
    stringsAsFactors = FALSE
  )
}

#' Extract Discordant Pairs and Shadow (Mate-Unmapped) Reads
#' @param dt data.frame of BAM records
#' @param expected_isize_max Insert size threshold for anomaly detection
#' @return data.frame of paired-end evidence edges
.extract_discordant_and_shadow <- function(dt, expected_isize_max = 1000L) {
  if (nrow(dt) == 0L) return(data.frame())
  out <- list()

  disc_ok <- !is.na(dt$mpos) & (abs(dt$isize) > expected_isize_max | dt$isize < 0)
  if (any(disc_ok)) {
    d <- dt[disc_ok, , drop = FALSE]
    out[[length(out) + 1L]] <- data.frame(
      qname = d$qname,
      from_coord = pmax(d$pos, d$mpos), to_coord = pmin(d$pos, d$mpos),
      prob = pmin(0.85, pmax(0.15, (d$mapq / 60) * 0.80)),
      evidence_type = "discordant_pair",
      stringsAsFactors = FALSE
    )
  }

  shadow_ok <- bitwAnd(dt$flag, 8L) > 0L
  if (any(shadow_ok)) {
    s <- dt[shadow_ok, , drop = FALSE]
    out[[length(out) + 1L]] <- data.frame(
      qname = s$qname, from_coord = s$pos, to_coord = s$pos,
      prob = 0.40, evidence_type = "shadow_read",
      stringsAsFactors = FALSE
    )
  }

  if (length(out) == 0L) return(data.frame())
  do.call(rbind, out)
}

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

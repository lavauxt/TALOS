# ============================================================================
# TALOS – BAM I/O and read pre-filtering
# ============================================================================

# ---------------------------------------------------------------------------
# Streaming BAM loader (with optional read-pair pass and coverage)
# ---------------------------------------------------------------------------
.load_bam_data_streaming <- function(bam_path, gene_config,
                                      compute_pairs = FALSE,
                                      max_reads = NULL, chunk_size = 50000L,
                                      verbose = FALSE) {
  if (verbose) message(
    "[TALOS] Loading BAM region: ",
    gene_config$chrom, ":", gene_config$genomic_start, "-", gene_config$genomic_end
  )

  target_chrom <- gene_config$chrom
  which_range  <- GenomicRanges::GRanges(
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

  .read_chunks <- function(bam_file, max_n) {
    chunk_list  <- list()
    total       <- 0L
    chunk_count <- 0L
    open(bam_file)
    on.exit(close(bam_file), add = TRUE)
    repeat {
      chunk <- GenomicAlignments::readGAlignments(bam_file, param = param,
                                                  use.names = TRUE)
      if (length(chunk) == 0L) break
      chunk_count <- chunk_count + 1L
      total       <- total + length(chunk)
      if (verbose) message("[TALOS] Loaded chunk ", chunk_count,
                           " (", length(chunk), " reads, total ", total, ")")
      chunk_list[[chunk_count]] <- chunk
      if (!is.null(max_n) && total >= max_n) {
        keep <- max_n - (total - length(chunk))
        chunk_list[[chunk_count]] <- chunk[seq_len(keep)]
        if (verbose) message("[TALOS] Stopping early at max_reads = ", max_n)
        break
      }
    }
    if (length(chunk_list) == 0L) NULL else do.call(c, chunk_list)
  }

  bam_file  <- Rsamtools::BamFile(bam_path, yieldSize = chunk_size)
  all_reads <- .read_chunks(bam_file, max_reads)

  if (verbose) message("[TALOS] Total reads loaded: ", length(all_reads %||% integer(0L)))

  if (is.null(all_reads) || length(all_reads) == 0L) {
    if (verbose) message("[TALOS] No reads found in target region.")
    return(list(reads = GenomicAlignments::GAlignments(), cov = NULL,
                pairs = NULL, cov_list = NULL))
  }

  # ---- Robust chromosome name normalisation ----
  current_levels <- GenomeInfoDb::seqlevels(all_reads)
  if (!target_chrom %in% current_levels) {
    converted <- GenomeInfoDb::mapSeqlevels(current_levels, "UCSC")
    if (any(!is.na(converted))) {
      new_levels <- ifelse(is.na(converted), current_levels, converted)
      all_reads  <- GenomeInfoDb::renameSeqlevels(all_reads, new_levels)
      if (!target_chrom %in% GenomeInfoDb::seqlevels(all_reads))
        warning("Chromosome '", target_chrom, "' not found after seqlevel conversion. Coverage will be NA.")
    } else {
      warning("Chromosome '", target_chrom, "' not found in BAM seqlevels. Coverage will be NA.")
    }
  }

  if (verbose) message("[TALOS] Computing coverage...")
  cov_rle_list <- GenomicAlignments::coverage(all_reads)
  cov          <- cov_rle_list[[target_chrom]]

  # ---- Optional paired-end pass (for discordant ratio) ----
  pairs <- NULL
  if (compute_pairs) {
    if (verbose) message("[TALOS] Loading read pairs...")
    bam_pairs <- Rsamtools::BamFile(bam_path, yieldSize = chunk_size)
    pair_list   <- list()
    total_pairs <- 0L
    chunk_count <- 0L
    open(bam_pairs)
    on.exit(close(bam_pairs), add = TRUE)   
    repeat {
      p_chunk <- suppressWarnings(
        GenomicAlignments::readGAlignmentPairs(bam_pairs, param = param)
      )
      if (length(p_chunk) == 0L) break
      chunk_count <- chunk_count + 1L
      total_pairs <- total_pairs + length(p_chunk)
      if (verbose) message("[TALOS] Loaded pair chunk ", chunk_count,
                           " (", length(p_chunk), " pairs, total ", total_pairs, ")")
      pair_list[[chunk_count]] <- p_chunk
      if (!is.null(max_reads) && total_pairs >= max_reads) {
        keep <- max_reads - (total_pairs - length(p_chunk))
        pair_list[[chunk_count]] <- p_chunk[seq_len(keep)]
        break
      }
    }
    if (length(pair_list) > 0L) pairs <- do.call(c, pair_list)
    if (verbose) message("[TALOS] Total pairs loaded: ", length(pairs %||% integer(0L)))
  }

  list(reads = all_reads, cov = cov, pairs = pairs, cov_list = cov_rle_list)
}


# ---------------------------------------------------------------------------
# Pre-filter reads: keep only those with soft-clips or net insertions ≥ threshold
# ---------------------------------------------------------------------------
.filter_reads_by_cigar <- function(reads, min_mapq, min_ins_filter,
                                    verbose = FALSE) {
  if (verbose) message("[TALOS] Pre-filtering reads: min_mapq=", min_mapq,
                       ", min_ins_filter=", min_ins_filter)

  flags  <- S4Vectors::mcols(reads)$flag
  mapqs  <- S4Vectors::mcols(reads)$mapq
  mapqs[is.na(mapqs)] <- 0L
  cigars   <- GenomicAlignments::cigar(reads)
  unmapped <- bitwAnd(flags, 0x4) != 0L

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
  if (verbose) message("[TALOS] ", sum(keep), " reads passed pre-filter (",
                       length(reads), " total)")
  filtered
}
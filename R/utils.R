# ── Null-coalescing operator (used throughout the package) ───────────────────
`%||%` <- function(x, y) if (is.null(x)) y else x

# ── Shared helper for softclip extraction ───────────────────────────────────
.get_softclips <- function(cig, r_seq) {
  lead <- NA_character_; trail <- NA_character_
  if (is.na(cig) || is.na(r_seq) || !grepl("S", cig, fixed = TRUE)) return(list(lead=lead, trail=trail))
  ops <- strsplit(cig, "(?<=[MIDNSHP=X])", perl = TRUE)[[1L]]
  ops <- ops[nzchar(ops)]
  if (length(ops) >= 1L && grepl("S$", ops[1L])) {
    op_len <- as.integer(sub("S$", "", ops[1L]))
    if (!is.na(op_len) && op_len > 0 && op_len <= nchar(r_seq)) lead <- substr(r_seq, 1L, op_len)
  }
  if (length(ops) >= 2L && grepl("S$", ops[length(ops)])) {
    op_len <- as.integer(sub("S$", "", ops[length(ops)]))
    if (!is.na(op_len) && op_len > 0 && op_len <= nchar(r_seq)) trail <- substr(r_seq, nchar(r_seq) - op_len + 1L, nchar(r_seq))
  }
  list(lead = lead, trail = trail)
}

# ── 3-letter amino acid lookup for HGVS protein notation ────────────────────
.aa1_to_aa3 <- c(
  A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys",
  Q = "Gln", E = "Glu", G = "Gly", H = "His", I = "Ile",
  L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro",
  S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val",
  `*` = "Ter"
)

.to_three_letter <- function(aa1_str) {
  chars <- strsplit(aa1_str, "")[[1L]]
  paste(vapply(chars, function(ch) .aa1_to_aa3[[ch]] %||% ch, character(1L)),
        collapse = "")
}

#' Annotate ITDs with known hotspot information
#' @param itd_df Data frame with columns \code{Gene}, \code{GenomicPosition}, \code{Length}.
#' @param db_path Optional path to CSV file (forces CSV mode).
#' @param genome_build \code{"hg19"} or \code{"hg38"}.
#' @param db_conn Optional DBI connection (overrides CSV).
#' @param fallback_csv_path Optional fallback CSV path.
#' @return Data frame with added columns \code{Hotspot} and \code{HotspotName}.
#' @export
annotate_hotspots <- function(itd_df, db_path = NULL, genome_build = NULL,
                               db_conn = NULL, fallback_csv_path = NULL) {
  itd_df$Hotspot     <- FALSE
  itd_df$HotspotName <- NA_character_

  if (!is.null(db_conn)) {
    if (!requireNamespace("DBI", quietly = TRUE))
      stop("Package 'DBI' required for database mode.")

    genes_in_df <- unique(itd_df$Gene)
    hs_list <- lapply(genes_in_df, function(g) {
      build_clause <- if (!is.null(genome_build)) " AND Build = ?" else ""
      params <- if (!is.null(genome_build)) list(g, genome_build) else list(g)
      DBI::dbGetQuery(db_conn,
        paste0("SELECT Gene, Start, End, Name FROM hotspots WHERE Gene = ?", build_clause),
        params = params
      )
    })
    hotspots_all <- do.call(rbind, hs_list)

    if (!is.null(hotspots_all) && nrow(hotspots_all) > 0) {
      hs_gr <- GenomicRanges::GRanges(
        seqnames = hotspots_all$Gene,
        ranges = IRanges::IRanges(start = hotspots_all$Start, end = hotspots_all$End),
        Name = hotspots_all$Name
      )
      itd_end <- itd_df$GenomicPosition +
                 ifelse(is.na(itd_df$Length), 0L, itd_df$Length) - 1L
      itd_end <- pmax(itd_end, itd_df$GenomicPosition)
      itd_gr <- GenomicRanges::GRanges(
        seqnames = itd_df$Gene,
        ranges = IRanges::IRanges(start = itd_df$GenomicPosition, end = itd_end)
      )
      hits <- GenomicRanges::findOverlaps(itd_gr, hs_gr, type = "any")
      if (length(hits) > 0) {
        qh <- S4Vectors::queryHits(hits)
        sh <- S4Vectors::subjectHits(hits)
        for (qi in unique(qh)) {
          names_i <- paste(unique(hotspots_all$Name[sh[qh == qi]]), collapse = ";")
          itd_df$Hotspot[qi]     <- TRUE
          itd_df$HotspotName[qi] <- names_i
        }
      }
    }
    return(itd_df)
  }

  csv_path <- db_path
  if (is.null(csv_path))
    csv_path <- system.file("extdata", "hotspots.csv", package = "TALOS")

  if (!is.null(csv_path) && file.exists(csv_path)) {
    hotspots <- read.csv(csv_path, stringsAsFactors = FALSE)
    if (!is.null(genome_build) && "Build" %in% colnames(hotspots))
      hotspots <- hotspots[hotspots$Build == genome_build, ]

    if (nrow(hotspots) == 0L) return(itd_df)

    hs_gr <- GenomicRanges::GRanges(
      seqnames = hotspots$Gene,
      ranges = IRanges::IRanges(start = hotspots$Start, end = hotspots$End),
      Name = hotspots$Name, Gene = hotspots$Gene
    )

    itd_end <- itd_df$GenomicPosition +
               ifelse(is.na(itd_df$Length), 0L, itd_df$Length) - 1L
    itd_end <- pmax(itd_end, itd_df$GenomicPosition)

    itd_gr <- GenomicRanges::GRanges(
      seqnames = itd_df$Gene,
      ranges = IRanges::IRanges(
        start = itd_df$GenomicPosition,
        end   = itd_end
      )
    )

    hits <- GenomicRanges::findOverlaps(itd_gr, hs_gr, type = "any")
    if (length(hits) > 0) {
      qh <- S4Vectors::queryHits(hits)
      sh <- S4Vectors::subjectHits(hits)
      for (qi in unique(qh)) {
        names_i <- paste(unique(hotspots$Name[sh[qh == qi]]), collapse = ";")
        itd_df$Hotspot[qi]     <- TRUE
        itd_df$HotspotName[qi] <- names_i
      }
    }
    return(itd_df)
  }

  if (!is.null(fallback_csv_path) && file.exists(fallback_csv_path))
    return(annotate_hotspots(itd_df, db_path = fallback_csv_path,
                             genome_build = genome_build))

  message("No hotspot data available. Hotspot annotation skipped.")
  return(itd_df)
}

.annotate_exonic_region <- function(final_df, exons_gr) {
  if (is.null(exons_gr) || nrow(final_df) == 0) {
    final_df$Region     <- NA_character_
    final_df$ExonNumber <- NA_integer_
    return(final_df)
  }

  exon_num_vec <- NULL
  if (!is.null(exons_gr$exon_number) && length(exons_gr$exon_number) > 0)
    exon_num_vec <- as.integer(unlist(exons_gr$exon_number))
  else if (!is.null(exons_gr$exon_idx) && length(exons_gr$exon_idx) > 0)
    exon_num_vec <- as.integer(unlist(exons_gr$exon_idx))
  else if (!is.null(names(exons_gr)))
    exon_num_vec <- suppressWarnings(as.integer(names(exons_gr)))

  bp_gr <- GenomicRanges::GRanges(
    seqnames = GenomeInfoDb::seqnames(exons_gr)[1L],
    ranges   = IRanges::IRanges(start = final_df$GenomicPosition,
                                end   = final_df$GenomicPosition)
  )
  hits <- GenomicRanges::findOverlaps(bp_gr, exons_gr, type = "within")
  final_df$Region     <- "intronic"
  final_df$ExonNumber <- NA_integer_

  if (length(hits) > 0) {
    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)
    keep  <- !duplicated(q_idx)
    q_idx <- q_idx[keep];  s_idx <- s_idx[keep]
    final_df$Region[q_idx] <- "exonic"
    if (!is.null(exon_num_vec)) final_df$ExonNumber[q_idx] <- exon_num_vec[s_idx]
  }
  final_df
}


compute_hgvs_annotations <- function(gene_config, genomic_pos,
                                      dup_seq  = NA_character_,
                                      dup_len  = NA_integer_,
                                      debug    = FALSE) {

  genomic_pos <- suppressWarnings(as.numeric(genomic_pos))
  if (is.na(genomic_pos))
    return(list(c_notation = NA_character_, p_notation = NA_character_))

  # ---- Check if transcript reference is available ----
  if (is.null(gene_config$transcript) || is.na(gene_config$transcript) ||
      length(gene_config$transcript) == 0 || gene_config$transcript == "") {
    return(list(c_notation = NA_character_, p_notation = NA_character_))
  }

  full_cdna <- gene_config$full_cdna
  if (is.null(full_cdna) || is.na(full_cdna) || nchar(full_cdna) == 0) {
    if (!is.null(gene_config$all_exons) && !is.null(gene_config$bsgenome_obj)) {
      if (debug) message("[HGVS] Rebuilding full_cdna from exons and BSgenome")
      exons_gr <- gene_config$all_exons
      strand   <- gene_config$strand
      bs_obj   <- gene_config$bsgenome_obj
      exon_seq_list <- lapply(seq_along(exons_gr), function(i) {
        seq <- as.character(BSgenome::getSeq(bs_obj, exons_gr[i]))
        if (strand == "-")
          seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(seq)))
        seq
      })
      full_cdna <- paste(unlist(exon_seq_list), collapse = "")
      if (nchar(full_cdna) == 0) {
        warning("Could not rebuild full_cdna. cDNA annotation may be incomplete.")
        full_cdna <- NA_character_
      }
    } else {
      warning("full_cdna missing and cannot be rebuilt. cDNA annotation may be incomplete.")
      full_cdna <- NA_character_
    }
  }

  seq_len <- if (!is.na(dup_seq) && nchar(dup_seq) > 0) nchar(dup_seq) else NA_integer_
  dup_len <- if (!is.na(dup_len) && dup_len > 0) as.integer(dup_len)
             else if (!is.na(seq_len)) seq_len
             else NA_integer_
  if (is.na(dup_len) || dup_len == 0L)
    return(list(c_notation = NA_character_, p_notation = NA_character_))

  exons_gr <- gene_config$all_exons
  if (is.null(exons_gr) || length(exons_gr) == 0)
    return(list(c_notation = NA_character_, p_notation = NA_character_))

  strand <- gene_config$strand

  if (!is.null(exons_gr$exon_idx)) {
    exons_gr <- exons_gr[order(exons_gr$exon_idx)]
  } else {
    if (strand == "+")
      exons_gr <- exons_gr[order(GenomicRanges::start(exons_gr))]
    else
      exons_gr <- exons_gr[order(GenomicRanges::start(exons_gr), decreasing = TRUE)]
  }

  exons_df <- data.frame(
    start    = GenomicRanges::start(exons_gr),
    end      = GenomicRanges::end(exons_gr),
    width    = GenomicRanges::width(exons_gr),
    exon_idx = if (!is.null(exons_gr$exon_idx)) exons_gr$exon_idx
               else seq_along(exons_gr),
    stringsAsFactors = FALSE
  )
  exons_df$cum_len_before <- c(0L, cumsum(exons_df$width)[-nrow(exons_df)])

  cds_offset <- if (is.null(gene_config$cds_offset)) 0L else gene_config$cds_offset

  genomic_end <- if (strand == "+") genomic_pos + dup_len - 1L
                 else               genomic_pos - dup_len + 1L

  .genomic_to_cdna_hgvs <- function(pos, exons_df, strand, cds_offset) {

    exon_hit <- which(exons_df$start <= pos & exons_df$end >= pos)
    if (length(exon_hit) > 0) {
      idx <- exon_hit[1L]
      pos_in_exon <- if (strand == "+")
        pos - exons_df$start[idx] + 1L
      else
        exons_df$end[idx] - pos + 1L
      c_pos <- exons_df$cum_len_before[idx] + pos_in_exon - cds_offset
      if (c_pos <= 0L) return(NA_character_)
      return(as.character(c_pos))
    }

    if (strand == "+") {
      up_mask <- exons_df$end   < pos
      dn_mask <- exons_df$start > pos
    } else {
      up_mask <- exons_df$start > pos
      dn_mask <- exons_df$end   < pos
    }

    up_note <- NULL; up_off <- Inf
    dn_note <- NULL; dn_off <- Inf

    if (any(up_mask)) {
      up_idx  <- tail(which(up_mask), 1L)
      up_c    <- exons_df$cum_len_before[up_idx] + exons_df$width[up_idx] - cds_offset
      up_off  <- if (strand == "+") pos - exons_df$end[up_idx]
                 else exons_df$start[up_idx] - pos
      up_note <- paste0(up_c, "+", up_off)
    }
    if (any(dn_mask)) {
      dn_idx  <- head(which(dn_mask), 1L)
      dn_c    <- exons_df$cum_len_before[dn_idx] + 1L - cds_offset
      dn_off  <- if (strand == "+") exons_df$start[dn_idx] - pos
                 else pos - exons_df$end[dn_idx]
      dn_note <- paste0(dn_c, "-", dn_off)
    }

    if (is.null(up_note) && is.null(dn_note)) return(paste0("?+", pos))
    if (is.null(up_note)) return(dn_note)
    if (is.null(dn_note)) return(up_note)
    if (up_off <= dn_off) up_note else dn_note
  }

  start_notation <- .genomic_to_cdna_hgvs(genomic_pos, exons_df, strand, cds_offset)
  end_notation   <- .genomic_to_cdna_hgvs(genomic_end,  exons_df, strand, cds_offset)

  if (is.na(start_notation) || is.na(end_notation)) {

    return(list(c_notation = "c.?+?", p_notation = "intronic"))
  }

  c_notation <- paste0("c.", start_notation, "_", end_notation, "dup")

  p_notation <- NA_character_
   if (!is.na(c_notation) &&
      !is.na(dup_len) && !is.na(full_cdna) && nchar(full_cdna) > cds_offset + 3L) {

    coding_seq <- substr(full_cdna, cds_offset + 1L, nchar(full_cdna))
    start_c    <- suppressWarnings(as.numeric(start_notation))
    end_c      <- suppressWarnings(as.numeric(end_notation))

    if (!is.na(start_c) && !is.na(end_c) && start_c >= 1L && end_c >= 1L) {
      c_start <- as.integer(start_c)
      c_end   <- as.integer(end_c)

      if (c_start >= 1L && c_end <= nchar(coding_seq)) {
        dup_cds <- substr(coding_seq, c_start, c_end)

        if (nchar(dup_cds) >= 3L) {
          trim_len <- nchar(dup_cds) - (nchar(dup_cds) %% 3L)

          aa_start_pos <- ceiling(c_start / 3L)
          aa_end_pos   <- ceiling(c_end   / 3L)

          cds_trim <- substr(coding_seq, 1L,
                             nchar(coding_seq) - nchar(coding_seq) %% 3L)
          if (nchar(cds_trim) >= 3L) {
            full_aa <- suppressWarnings(
              as.character(Biostrings::translate(Biostrings::DNAString(cds_trim)))
            )
            if (aa_end_pos <= nchar(full_aa)) {
              aa1_1 <- substr(full_aa, aa_start_pos, aa_start_pos)
              aa1_2 <- substr(full_aa, aa_end_pos,   aa_end_pos)
              aa3_1 <- .to_three_letter(aa1_1)
              aa3_2 <- .to_three_letter(aa1_2)

              if (nchar(dup_cds) %% 3L != 0L) {
                p_notation <- paste0("p.", aa3_1, aa_start_pos, "fsTer?")
              } else if (aa_start_pos == aa_end_pos) {
                p_notation <- paste0("p.", aa3_1, aa_start_pos, "dup")
              } else {
                p_notation <- paste0("p.", aa3_1, aa_start_pos,
                                     "_", aa3_2, aa_end_pos, "dup")
              }
            }
          }
        }
      }
    }
  }


  if (is.na(p_notation)) {
    if (grepl("[+-]", start_notation)) {
      p_notation <- "intronic"
    }
  }

  list(c_notation = c_notation, p_notation = p_notation)
}


compute_coverage_drop <- function(cov_all, chrom, breakpoint,
                                   flank = 200, min_depth = 5, 
                                   min_flank_bp = 5,
                                   verbose = FALSE) {
  if (is.null(cov_all)) {
    if (verbose) message("Coverage drop: cov_all is NULL")
    return(NA_real_)
  }

  if (inherits(cov_all, "Rle") || is.numeric(cov_all)) {
    cov_rle <- cov_all
  } else if (chrom %in% names(cov_all)) {
    cov_rle <- cov_all[[chrom]]
  } else {
    if (verbose) message("Coverage drop: chrom missing in cov_all list")
    return(NA_real_)
  }

  cov_len   <- length(cov_rle)
  start_pos <- max(1L, breakpoint - flank)
  end_pos   <- min(cov_len, breakpoint + flank)
  if (start_pos >= end_pos) {
    if (verbose) message("Coverage drop: start_pos >= end_pos")
    return(NA_real_)
  }
  
  cov_vec         <- as.numeric(cov_rle[start_pos:end_pos])
  bp_relative     <- breakpoint - start_pos + 1L
  left_end_rel    <- max(1L, bp_relative - 10L)
  right_start_rel <- min(length(cov_vec), bp_relative + 10L)
  
  if (left_end_rel < min_flank_bp || (length(cov_vec) - right_start_rel) < min_flank_bp) {
    if (verbose) message("Coverage drop: insufficient flanking bases (left=", left_end_rel,
                         ", right=", length(cov_vec) - right_start_rel, ")")
    return(NA_real_)
  }
  
  left_cov  <- mean(cov_vec[1L:left_end_rel],                 na.rm = TRUE)
  right_cov <- mean(cov_vec[right_start_rel:length(cov_vec)], na.rm = TRUE)
  
  if (left_cov < min_depth || right_cov < min_depth) {
    if (verbose) message("Coverage drop: left_cov=", left_cov, " right_cov=", right_cov, " below min_depth=", min_depth)
    return(NA_real_)
  }
  
  ratio <- max(left_cov, right_cov) / min(left_cov, right_cov)
  if (verbose) message("Coverage drop ratio = ", ratio)
  ratio
}

compute_local_coverage <- function(cov_all, chrom, breakpoint,
                                    buffer = 10L, verbose = FALSE) {
  if (is.null(cov_all)) {
    if (verbose) message("Local coverage: cov_all is NULL")
    return(NA_real_)
  }

  if (inherits(cov_all, "Rle") || is.numeric(cov_all)) {
    cov_rle <- cov_all
  } else if (chrom %in% names(cov_all)) {
    cov_rle <- cov_all[[chrom]]
  } else {
    if (verbose) message("Local coverage: chrom missing in cov_all list")
    return(NA_real_)
  }

  cov_len   <- length(cov_rle)
  start_pos <- breakpoint - buffer
  end_pos   <- breakpoint + buffer
  if (start_pos < 1L || end_pos > cov_len || start_pos > end_pos) {
    if (verbose) message("Local coverage: window out of range")
    return(NA_real_)
  }

  min_cov <- min(as.numeric(cov_rle[start_pos:end_pos]))
  if (verbose) message("Local coverage (min over bp+/-", buffer, ") = ", min_cov)
  min_cov
}

compute_span_coverage <- function(cov_all, chrom, start, span_len,
                                   verbose = FALSE) {
  empty_result <- list(min_cov = NA_real_, mean_cov = NA_real_)
  if (is.null(cov_all) || is.na(start) || is.na(span_len) || span_len < 1L) {
    if (verbose) message("Span coverage: missing inputs")
    return(empty_result)
  }

  if (inherits(cov_all, "Rle") || is.numeric(cov_all)) {
    cov_rle <- cov_all
  } else if (chrom %in% names(cov_all)) {
    cov_rle <- cov_all[[chrom]]
  } else {
    if (verbose) message("Span coverage: chrom missing in cov_all list")
    return(empty_result)
  }

  cov_len   <- length(cov_rle)
  start_pos <- start
  end_pos   <- start + span_len - 1L
  if (start_pos < 1L || end_pos > cov_len || start_pos > end_pos) {
    if (verbose) message("Span coverage: span out of range")
    return(empty_result)
  }

  cov_vec <- as.numeric(cov_rle[start_pos:end_pos])
  result  <- list(min_cov = min(cov_vec), mean_cov = mean(cov_vec))
  if (verbose) message("Span coverage over ", span_len, "bp: min=", result$min_cov,
                        " mean=", round(result$mean_cov, 1))
  result
}

compute_span_depth_fold_change <- function(cov_all, chrom, start, span_len,
                                            flank = 300L, verbose = FALSE) {
  if (is.null(cov_all) || is.na(start) || is.na(span_len) || span_len < 1L) {
    return(NA_real_)
  }

  if (inherits(cov_all, "Rle") || is.numeric(cov_all)) {
    cov_rle <- cov_all
  } else if (chrom %in% names(cov_all)) {
    cov_rle <- cov_all[[chrom]]
  } else {
    return(NA_real_)
  }

  cov_len    <- length(cov_rle)
  span_start <- start
  span_end   <- start + span_len - 1L
  if (span_start < 1L || span_end > cov_len || span_start > span_end) {
    return(NA_real_)
  }

  left_start  <- max(1L, span_start - flank)
  left_end    <- span_start - 1L
  right_start <- span_end + 1L
  right_end   <- min(cov_len, span_end + flank)

  flank_vals <- c(
    if (left_end >= left_start) as.numeric(cov_rle[left_start:left_end]) else numeric(0),
    if (right_end >= right_start) as.numeric(cov_rle[right_start:right_end]) else numeric(0)
  )
  if (length(flank_vals) == 0) return(NA_real_)

  background_cov <- mean(flank_vals)
  if (is.na(background_cov) || background_cov <= 0) return(NA_real_)

  span_cov <- mean(as.numeric(cov_rle[span_start:span_end]))
  span_cov / background_cov
}


compute_discordant_ratio <- function(all_pairs, breakpoint, flank = 5000,
                                     min_mapq = 20, insert_size_factor = 2.0,
                                     min_everted_separation = 100L) {
  if (is.null(all_pairs) || length(all_pairs) == 0) return(NA_real_)
  
  start_region <- max(1L, breakpoint - flank)
  end_region   <- breakpoint + flank
  
  pairs_in_window <- all_pairs[
    GenomicRanges::start(GenomicAlignments::first(all_pairs)) >= start_region &
    GenomicRanges::end(GenomicAlignments::first(all_pairs))   <= end_region
  ]
  if (length(pairs_in_window) == 0) return(NA_real_)
  
  r1 <- GenomicAlignments::first(pairs_in_window)
  r2 <- GenomicAlignments::last(pairs_in_window)
  
  flag1 <- S4Vectors::mcols(r1)$flag; if (is.null(flag1)) flag1 <- rep(0L, length(r1))
  flag2 <- S4Vectors::mcols(r2)$flag; if (is.null(flag2)) flag2 <- rep(0L, length(r2))
  mapq1 <- S4Vectors::mcols(r1)$mapq; if (is.null(mapq1)) mapq1 <- rep(0L, length(r1))
  mapq2 <- S4Vectors::mcols(r2)$mapq; if (is.null(mapq2)) mapq2 <- rep(0L, length(r2))
  
  mapq1[is.na(mapq1)] <- 0L; mapq2[is.na(mapq2)] <- 0L
  
  valid_mask <- mapq1 >= min_mapq & mapq2 >= min_mapq &
                bitwAnd(flag1, 0x100) == 0 & bitwAnd(flag1, 0x800) == 0 &
                bitwAnd(flag2, 0x100) == 0 & bitwAnd(flag2, 0x800) == 0
                
  if (sum(valid_mask) == 0) return(NA_real_)
  
  start1 <- GenomicRanges::start(r1)[valid_mask]
  start2 <- GenomicRanges::start(r2)[valid_mask]
  
  end1 <- GenomicRanges::end(r1)[valid_mask]
  end2 <- GenomicRanges::end(r2)[valid_mask]
  mapped_isize <- pmax(end1, end2) - pmin(start1, start2)
  
  med_insert <- median(mapped_isize, na.rm = TRUE)
  sd_insert  <- sd(mapped_isize, na.rm = TRUE)
  if (is.na(sd_insert)) sd_insert <- 0
  
  large_insert <- mapped_isize > (med_insert + (insert_size_factor * sd_insert))
  
  spanning <- (start1 <= breakpoint & start2 > breakpoint) |
              (start2 <= breakpoint & start1 > breakpoint)
  
  rev1_full <- bitwAnd(flag1, 0x10) != 0
  rev2_full <- bitwAnd(flag2, 0x10) != 0
  rev1 <- rev1_full[valid_mask]
  rev2 <- rev2_full[valid_mask]
  
  left_is_r1 <- start1 <= start2
  left_rev   <- ifelse(left_is_r1, rev1, rev2)
  right_rev  <- ifelse(left_is_r1, rev2, rev1)
  mate_sep   <- abs(start1 - start2)
  everted    <- left_rev & !right_rev & mate_sep >= min_everted_separation
  
  discordant_spanning <- sum(spanning & (large_insert | everted))
  total_valid_pairs   <- sum(valid_mask)
  
  discordant_spanning / total_valid_pairs
}


#' Estimate ITD length from paired-end insert sizes
#'
#' For ITDs larger than the read length, neither CIGAR insertions nor k-mer
#' backward jumps can span the full duplication in a single read.  Read pairs
#' that straddle the breakpoint carry an inflated outer insert size equal to
#' the normal insert size plus the ITD length.  This function estimates the
#' ITD length as:
#'   \code{median(spanning-pair outer insert) - median(background outer insert)}
#'
#' This estimate is complementary to the k-mer length; it is more reliable when
#' the ITD exceeds the read length and should be preferred in that regime.
#' It is stored in the \code{LengthPE} output column and does NOT replace the
#' primary \code{Length} value (which is always CIGAR-based when available,
#' k-mer-based otherwise).
#'
#' @param all_pairs A \code{GAlignmentPairs} object.
#' @param breakpoint Genomic breakpoint position (integer).
#' @param nominal_read_len Typical read length used to determine whether the
#'   ITD is likely larger than a single read (default 150).
#' @param flank Flanking window (bp) for collecting background pairs (default 5000).
#' @param min_mapq Minimum MAPQ for both reads of a pair (default 20).
#' @param min_spanning Minimum spanning pairs required for a valid estimate
#'   (default 5).
#' @param max_insert Maximum plausible outer insert size to include (default 5000).
#' @return Named list with elements \code{length_pe} (estimated length,
#'   \code{NA_real_} if insufficient data), \code{n_spanning} (number of
#'   spanning pairs used), and \code{n_background} (number of background
#'   pairs used).
#' @export
compute_pe_itd_length <- function(all_pairs, breakpoint,
                                   nominal_read_len = 150L,
                                   flank            = 5000L,
                                   min_mapq         = 20L,
                                   min_spanning     = 5L,
                                   max_insert       = 5000L) {

  empty <- list(length_pe = NA_real_, n_spanning = 0L, n_background = 0L)
  if (is.null(all_pairs) || length(all_pairs) == 0L) return(empty)

  r1 <- GenomicAlignments::first(all_pairs)
  r2 <- GenomicAlignments::last(all_pairs)

  mq1 <- S4Vectors::mcols(r1)$mapq; mq1[is.na(mq1)] <- 0L
  mq2 <- S4Vectors::mcols(r2)$mapq; mq2[is.na(mq2)] <- 0L
  f1  <- S4Vectors::mcols(r1)$flag; if (is.null(f1)) f1 <- rep(0L, length(r1))
  f2  <- S4Vectors::mcols(r2)$flag; if (is.null(f2)) f2 <- rep(0L, length(r2))

  primary_mask <- bitwAnd(f1, 0x100L) == 0L & bitwAnd(f1, 0x800L) == 0L &
                  bitwAnd(f2, 0x100L) == 0L & bitwAnd(f2, 0x800L) == 0L
  mapq_mask    <- mq1 >= min_mapq & mq2 >= min_mapq
  valid        <- primary_mask & mapq_mask

  s1 <- GenomicRanges::start(r1); e2 <- GenomicRanges::end(r2)
  
  pair_start <- pmin(s1, GenomicRanges::start(r2))
  pair_end   <- pmax(GenomicRanges::end(r1), e2)
  isize      <- pair_end - pair_start + 1L

  size_mask <- isize > 0L & isize <= max_insert

 
  bg_mask <- valid & size_mask &
             pair_start >= (breakpoint - flank) &
             pair_end   <= (breakpoint + flank) &
             !(pair_start < breakpoint & pair_end > breakpoint)

  span_mask <- valid & size_mask &
               pair_start <= (breakpoint - nominal_read_len) &
               pair_end   >= (breakpoint + nominal_read_len)

  n_span <- sum(span_mask, na.rm = TRUE)
  n_bg   <- sum(bg_mask,   na.rm = TRUE)

  if (n_span < min_spanning || n_bg < 10L) return(empty)

  med_span <- stats::median(isize[span_mask], na.rm = TRUE)
  med_bg   <- stats::median(isize[bg_mask],   na.rm = TRUE)

  est <- round(med_span - med_bg)
  if (is.na(est) || est <= 0L) return(empty)

  list(length_pe = est, n_spanning = n_span, n_background = n_bg)
}

compute_microhomology <- function(support_rows, ref_dna, breakpoint,
                                   genomic_start, debug = FALSE) {
  if (nrow(support_rows) == 0) return(NA_integer_)
  mh_lengths    <- integer()
  local_bp      <- breakpoint - genomic_start + 1L
  ref_len_total <- nchar(ref_dna)
  
  has_biostrings <- requireNamespace("Biostrings", quietly = TRUE)
  
  for (i in seq_len(nrow(support_rows))) {
    cig      <- support_rows$cigar[i]
    read_seq <- support_rows$read_seq[i]
    
    sc <- .get_softclips(cig, read_seq)
    
    if (!is.na(sc$lead)) {
      len <- nchar(sc$lead)
      ref_start <- local_bp + 1L
      ref_end   <- min(ref_len_total, local_bp + len)
      ref_down  <- substr(ref_dna, ref_start, ref_end)
      if (nchar(ref_down) >= len) {
        mm <- if (has_biostrings) {
          Biostrings::lcprefix(sc$lead, ref_down)
        } else {
          min_len <- min(nchar(sc$lead), nchar(ref_down))
          if (substr(sc$lead, 1, min_len) == substr(ref_down, 1, min_len)) min_len else 0
        }
        if (mm > 0) mh_lengths <- c(mh_lengths, mm)
      }
    }
    
    if (!is.na(sc$trail)) {
      len <- nchar(sc$trail)
      ref_start <- max(1L, local_bp - len)
      ref_end   <- local_bp - 1L
      ref_up    <- substr(ref_dna, ref_start, ref_end)
      if (nchar(ref_up) >= len) {
        mm <- if (has_biostrings) {
          Biostrings::lcprefix(sc$trail, ref_up)
        } else {
          min_len <- min(nchar(sc$trail), nchar(ref_up))
          if (substr(sc$trail, 1, min_len) == substr(ref_up, 1, min_len)) min_len else 0
        }
        if (mm > 0) mh_lengths <- c(mh_lengths, mm)
      }
    }
  }
  if (length(mh_lengths) == 0) return(NA_integer_)
  as.integer(median(mh_lengths))
}

compute_repeat_entropy <- function(ref_dna, breakpoint, genomic_start,
                                    window = 50) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    warning("Biostrings not available; repeat entropy set to 0")
    return(0)
  }
  local_pos  <- breakpoint - genomic_start + 1L
  start_pos  <- max(1L, local_pos - window)
  end_pos    <- min(nchar(ref_dna), local_pos + window)
  seq_window <- substr(ref_dna, start_pos, end_pos)
  if (nchar(seq_window) < 10L) return(0)
  freq <- Biostrings::oligonucleotideFrequency(Biostrings::DNAString(seq_window), 2L)
  freq <- freq / sum(freq)
  -sum(freq * log2(freq + 1e-12))
}

#' Detect orientation of a duplication (tandem vs inverted) using alignment
#' @param itd_seq Duplicated sequence (character)
#' @param ref_seq Reference segment (character)
#' @param min_pid Minimum percent identity to call (0-1)
#' @export
detect_orientation <- function(itd_seq, ref_seq, min_pid = 0.90) {
  if (is.na(itd_seq) || is.na(ref_seq) ||
      nchar(itd_seq) == 0 || nchar(ref_seq) == 0) return(NA_character_)
  itd_seq <- toupper(itd_seq); ref_seq <- toupper(ref_seq)

  has_pwalign   <- requireNamespace("pwalign",   quietly = TRUE)
  has_biostrings <- requireNamespace("Biostrings", quietly = TRUE)

  if (nchar(itd_seq) < 10 || nchar(ref_seq) < 10) {
    if (itd_seq == ref_seq) return("tandem")
    if (has_biostrings) {
      rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(ref_seq)))
      if (itd_seq == rc) return("inverted")
    }
    return(NA_character_)
  }

  if (has_pwalign && has_biostrings) {
    aln_fwd <- pwalign::pairwiseAlignment(itd_seq, ref_seq, type = "local")
    if (pwalign::pid(aln_fwd) >= min_pid * 100) return("tandem")
    rc <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(ref_seq))
    )
    aln_rev <- pwalign::pairwiseAlignment(itd_seq, rc, type = "local")
    if (pwalign::pid(aln_rev) >= min_pid * 100) return("inverted")
  } else if (has_biostrings) {
    if (itd_seq == ref_seq) return("tandem")
    rc <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(ref_seq))
    )
    if (itd_seq == rc) return("inverted")
  } else {
    if (itd_seq == ref_seq) return("tandem")
  }
  NA_character_
}


.resolve_ref_base <- function(pos, ref_dna, genomic_start) {
  local_pos <- pos - genomic_start + 1L
  ref_len   <- nchar(ref_dna)
  if (local_pos > 0L && local_pos <= ref_len) substr(ref_dna, local_pos, local_pos) else "N"
}

.format_info <- function(row) {
  cds_val <- if (is.na(row$HGVS_cDNA))    "." else row$HGVS_cDNA
  info    <- sprintf(
    "SVTYPE=DUP;SVLEN=%d;END=%d;GENE=%s;CDS=%s;AF=%.4f;DP=%d;SUPPORT=%d;WT=%d",
    row$Length, row$GenomicPosition + row$Length - 1L, row$Gene,
    cds_val, row$AlleleFrequency, row$DepthAtBreakpoint,
    row$SupportingReads, row$WildtypeReads
  )
  if (!is.na(row$HGVS_Protein) && row$HGVS_Protein != "")
    info <- paste0(info, sprintf(";AA=%s", row$HGVS_Protein))

  append_if <- function(info, cond, fmt, val) {
    if (cond) paste0(info, sprintf(fmt, val)) else info
  }
  info <- append_if(info, !is.null(row$StrandBias)           && !is.na(row$StrandBias),           ";SB=%.4f",       row$StrandBias)
  info <- append_if(info, !is.null(row$MeanSupportMAPQ)      && !is.na(row$MeanSupportMAPQ),      ";MMAPQ=%.1f",    row$MeanSupportMAPQ)
  info <- append_if(info, !is.null(row$BreakpointSpread)     && !is.na(row$BreakpointSpread),     ";BPSPREAD=%d",   as.integer(row$BreakpointSpread))
  info <- append_if(info, !is.null(row$SoftclipFraction)     && !is.na(row$SoftclipFraction),     ";SCFRAC=%.4f",   row$SoftclipFraction)
  info <- append_if(info, !is.null(row$UniqueBreakpoints)    && !is.na(row$UniqueBreakpoints),    ";UBP=%d",        as.integer(row$UniqueBreakpoints))
  info <- append_if(info, !is.null(row$CoverageDrop)         && !is.na(row$CoverageDrop),         ";CD=%.2f",       row$CoverageDrop)
  info <- append_if(info, !is.null(row$DiscordantRatio)      && !is.na(row$DiscordantRatio),      ";DRATIO=%.4f",   row$DiscordantRatio)
  info <- append_if(info, !is.null(row$MedianMicrohomology)  && !is.na(row$MedianMicrohomology),  ";MH=%d",         as.integer(row$MedianMicrohomology))
  info <- append_if(info, !is.null(row$RepeatEntropy)        && !is.na(row$RepeatEntropy),        ";ENT=%.2f",      row$RepeatEntropy)
  info <- append_if(info, !is.null(row$Orientation)          && !is.na(row$Orientation),          ";ORIENT=%s",     row$Orientation)
  info <- append_if(info, !is.null(row$SupportConsistency)   && !is.na(row$SupportConsistency),   ";CONSIST=%.4f",  row$SupportConsistency)
  info <- append_if(info, !is.null(row$RefMatch_Observed)    && !is.na(row$RefMatch_Observed),    ";REFMATCH_OBS=%.1f", row$RefMatch_Observed)
  info <- append_if(info, !is.null(row$RefMatch_Total)       && !is.na(row$RefMatch_Total),       ";REFMATCH_TOT=%.1f", row$RefMatch_Total)
  info <- append_if(info, !is.null(row$ITDReadCoverage)      && !is.na(row$ITDReadCoverage),      ";ITDCOV=%.1f",   row$ITDReadCoverage)
  if (!is.null(row$Hotspot) && isTRUE(row$Hotspot) && !is.na(row$HotspotName))
    info <- paste0(info, sprintf(";HOTSPOT=%s", row$HotspotName))
  info <- append_if(info, !is.null(row$AlignmentScore)       && !is.na(row$AlignmentScore),       ";ALNSCORE=%.4f", row$AlignmentScore)
  info <- append_if(info, !is.null(row$TotalSupportBases)    && !is.na(row$TotalSupportBases),    ";SUPBASES=%d",   as.integer(row$TotalSupportBases))
  info <- append_if(info, !is.null(row$LeftSoftclipCount)    && !is.na(row$LeftSoftclipCount),    ";LSC=%d",        row$LeftSoftclipCount)
  info <- append_if(info, !is.null(row$RightSoftclipCount)   && !is.na(row$RightSoftclipCount),   ";RSC=%d",        row$RightSoftclipCount)
  info <- append_if(info, !is.null(row$LeftSoftclipPctSupport) && !is.na(row$LeftSoftclipPctSupport), ";LSCPCT=%.1f", row$LeftSoftclipPctSupport)
  info <- append_if(info, !is.null(row$RightSoftclipPctSupport) && !is.na(row$RightSoftclipPctSupport), ";RSCPCT=%.1f", row$RightSoftclipPctSupport)
  info
}

#' Build VCF header lines
#' @export
build_vcf_header <- function(sample_name, chrom, genome_build = NULL, contig_length = NULL) {
  pkg_ver <- tryCatch(as.character(utils::packageVersion("TALOS")), error = function(e) "dev")
  lines   <- c(
    "##fileformat=VCFv4.2",
    paste0("##fileDate=",    format(Sys.Date(), "%Y%m%d")),
    paste0("##source=TALOS_v", pkg_ver),
    paste0("##SAMPLE=<ID=", sample_name, ">")
  )
  if (!is.null(genome_build))
    lines <- c(lines, paste0("##reference=", genome_build))
  
  if (!is.null(chrom)) {
    if (!is.null(contig_length) && is.numeric(contig_length) && contig_length > 0) {
      lines <- c(lines, paste0("##contig=<ID=", chrom, ",length=", contig_length, ">"))
    } else {
      lines <- c(lines, paste0("##contig=<ID=", chrom, ",length=?>"))
    }
  }
  
  lines <- c(lines,
    '##ALT=<ID=DUP,Description="Tandem duplication">',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of duplication (bp)">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of duplication">',
    '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">',
    '##INFO=<ID=CDS,Number=1,Type=String,Description="HGVS cDNA notation">',
    '##INFO=<ID=AA,Number=1,Type=String,Description="HGVS protein notation">',
    '##INFO=<ID=AF,Number=1,Type=Float,Description="Mutant allele frequency">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth at breakpoint">',
    '##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Bias-corrected fragment count">',
    '##INFO=<ID=WT,Number=1,Type=Integer,Description="Wildtype fragment count">',
    '##INFO=<ID=REFMATCH_OBS,Number=1,Type=Float,Description="% identity between observed ITD sequence and reference segment (0-100)">',
    '##INFO=<ID=REFMATCH_TOT,Number=1,Type=Float,Description="% identity between full reconstructed ITD sequence (observed+imputed) and reference segment (0-100)">',
    '##INFO=<ID=ITDCOV,Number=1,Type=Float,Description="% of ITD sequence positions covered by actual reads (0-100)">',
    '##INFO=<ID=CD,Number=1,Type=Float,Description="Coverage drop fold-change across breakpoint">',
    '##INFO=<ID=DRATIO,Number=1,Type=Float,Description="Discordant pair ratio spanning breakpoint">',
    '##INFO=<ID=MH,Number=1,Type=Integer,Description="Median microhomology length (bp)">',
    '##INFO=<ID=ENT,Number=1,Type=Float,Description="Dinucleotide Shannon entropy around breakpoint">',
    '##INFO=<ID=CONSIST,Number=1,Type=Float,Description="Support consistency (% of reads containing ITD sequence)">',
    '##INFO=<ID=HOTSPOT,Number=1,Type=String,Description="Known hotspot name if applicable">',
    '##INFO=<ID=ALNSCORE,Number=1,Type=Float,Description="Alignment score between ITD and reference (0-1)">',
    '##INFO=<ID=SUPBASES,Number=1,Type=Integer,Description="Total soft-clip bases from supporting reads">',
    '##INFO=<ID=SB,Number=1,Type=Float,Description="Strand bias: fraction of supporting reads on reverse strand">',
    '##INFO=<ID=MMAPQ,Number=1,Type=Float,Description="Mean MAPQ of supporting reads">',
    '##INFO=<ID=BPSPREAD,Number=1,Type=Integer,Description="Range of breakpoint positions within cluster (bp)">',
    '##INFO=<ID=SCFRAC,Number=1,Type=Float,Description="Fraction of supporting reads detected via soft-clip">',
    '##INFO=<ID=UBP,Number=1,Type=Integer,Description="Number of unique breakpoint positions in cluster">',
    '##INFO=<ID=ORIENT,Number=1,Type=String,Description="Duplication orientation: tandem or inverted">',
    '##INFO=<ID=LSC,Number=1,Type=Integer,Description="Number of supporting reads with left soft-clip">',
    '##INFO=<ID=RSC,Number=1,Type=Integer,Description="Number of supporting reads with right soft-clip">',
    '##INFO=<ID=LSCPCT,Number=1,Type=Float,Description="Percentage of supporting reads with left soft-clip">',
    '##INFO=<ID=RSCPCT,Number=1,Type=Float,Description="Percentage of supporting reads with right soft-clip">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name),
          collapse = "\t")
  )
  lines
}

#' Build VCF data lines from TALOS result data frame
#' @export
build_vcf_records <- function(itd_df, ref_dna, genomic_start, chrom) {
  if (nrow(itd_df) == 0L) return(character(0L))
  records <- character(nrow(itd_df))
  for (i in seq_len(nrow(itd_df))) {
    row <- itd_df[i, ]
    if (row$DepthAtBreakpoint == 0L) next
    pos      <- row$GenomicPosition
    ref_base <- .resolve_ref_base(pos, ref_dna, genomic_start)
    id_field <- if (!is.na(row$HGVS_cDNA)) row$HGVS_cDNA else "."
    info     <- .format_info(row)
    records[i] <- paste(chrom, pos, id_field, ref_base, "<DUP>", ".", "PASS",
                        info, "GT", "0/1", sep = "\t")
  }
  records[nzchar(records)]
}

#' Write a VCF file from TALOS results
#' @export
write_itd_vcf <- function(itd_df, ref_dna, genomic_start, chrom,
                           sample_name, genome_build = NULL,
                           vcf_path, overwrite = FALSE) {
  if (!overwrite && file.exists(vcf_path)) stop("VCF already exists: ", vcf_path)
  out_dir <- dirname(vcf_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  contig_length <- NULL
  if (!is.null(genome_build) && requireNamespace("BSgenome", quietly = TRUE)) {
    bs_pkg <- if (genome_build == "hg19") "BSgenome.Hsapiens.UCSC.hg19" else "BSgenome.Hsapiens.UCSC.hg38"
    if (requireNamespace(bs_pkg, quietly = TRUE)) {
      bs_obj <- BSgenome::getBSgenome(bs_pkg)
      if (chrom %in% names(bs_obj)) {
        contig_length <- length(bs_obj[[chrom]])
      }
    }
  }
  
  lines <- c(build_vcf_header(sample_name, chrom, genome_build, contig_length),
             build_vcf_records(itd_df, ref_dna, genomic_start, chrom))
  writeLines(lines, vcf_path)
  invisible(vcf_path)
}


.infer_sample_names <- function(bam_paths, part = 1L, sep = NULL) {
  basenames <- sub("\\.bam$", "", basename(bam_paths), ignore.case = TRUE)

  vapply(basenames, function(bn) {
    if (is.null(sep)) {
      sep_used <- if (grepl("\\.", bn)) "\\."
                  else if (grepl("_",  bn)) "_"
                  else if (grepl("-",  bn)) "-"
                  else NULL
    } else {
      sep_used <- sep
    }
    if (is.null(sep_used)) return(bn)
    parts   <- strsplit(bn, sep_used)[[1L]]
    n_parts <- length(parts)
    if (part > n_parts) {
      warning(sprintf("BAM '%s' has only %d part(s); using full name.", bn, n_parts))
      return(bn)
    }
    parts[part]
  }, character(1L), USE.NAMES = FALSE)
}

#' Batch TALOS: analyse every BAM in a directory
#' @export
talos_batch_dir <- function(
    bam_dir,
    genes,
    build            = "hg19",
    output_folder    = "./results",
    pattern          = "\\.bam$",
    sample_name_part = 1L,
    sample_name_sep  = NULL,
    recursive        = FALSE,
    combine_output   = TRUE,
    ...
) {
  if (!dir.exists(bam_dir)) stop("Directory not found: ", bam_dir)
  bam_files <- list.files(bam_dir, pattern = pattern,
                           full.names = TRUE, recursive = recursive)
  if (length(bam_files) == 0)
    stop("No BAM files found in '", bam_dir, "' matching pattern '", pattern, "'")

  sample_names <- .infer_sample_names(bam_files, part = sample_name_part,
                                       sep = sample_name_sep)
  message(sprintf("[TALOS batch_dir] Found %d BAM file(s) in '%s'.",
                  length(bam_files), bam_dir))
  for (i in seq_along(bam_files))
    message(sprintf("  %s  ->  sample: %s", basename(bam_files[i]), sample_names[i]))

  talos_batch(
    bam_paths     = bam_files,
    genes         = genes,
    build         = build,
    sample_names  = sample_names,
    output_folder = output_folder,
    combine_output = combine_output,
    ...
  )
}

#' Batch TALOS: analyse multiple genes and/or samples
#' @export
talos_batch <- function(
    bam_paths,
    genes,
    build          = "hg19",
    sample_names   = NULL,
    output_folder  = "./results",
    combine_output = TRUE,
    ...
) {
  .parse_input <- function(x) {
    if (length(x) == 1L && grepl(",", x)) trimws(strsplit(x, ",")[[1L]])
    else as.character(x)
  }

  bam_paths    <- .parse_input(bam_paths)
  genes        <- .parse_input(genes)
  sample_names <- if (!is.null(sample_names)) .parse_input(sample_names)
                  else .infer_sample_names(bam_paths)

  if (length(sample_names) != length(bam_paths))
    stop("talos_batch: length of sample_names must equal length of bam_paths.")

  batch_start <- Sys.time()
  message(sprintf("[TALOS batch] %d sample(s) x %d gene(s) = %d run(s)",
    length(bam_paths), length(genes), length(bam_paths) * length(genes)))

  all_results <- list()
  for (bi in seq_along(bam_paths)) {
    for (gi in seq_along(genes)) {
      bam <- bam_paths[bi]; gene <- genes[gi]; sn <- sample_names[bi]
      message(sprintf("[TALOS batch] Running: sample=%s  gene=%s", sn, gene))
      res <- tryCatch(
        talos(bam_path = bam, gene = gene, build = build,
              sample_name = sn, output_folder = output_folder, ...),
        error = function(e) {
          message(sprintf("[TALOS batch] ERROR for %s / %s: %s",
                          sn, gene, conditionMessage(e)))
          data.frame()
        }
      )
      if (is.data.frame(res) && nrow(res) > 0)
        all_results[[length(all_results) + 1L]] <- res
    }
  }

  combined <- if (length(all_results) > 0) do.call(rbind, all_results) else data.frame()

  if (combine_output && nrow(combined) > 0) {
    if (!dir.exists(output_folder))
      dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
    ts    <- format(Sys.time(), "%Y%m%d_%H%M%S")
    cfile <- file.path(output_folder, paste0("TALOS_batch_", ts, ".tsv"))
    write.table(combined, file = cfile, sep = "\t", quote = FALSE,
                row.names = FALSE, na = ".")
    message(sprintf("[TALOS batch] Combined results written to: %s", cfile))
  }

  elapsed <- as.numeric(difftime(Sys.time(), batch_start, units = "secs"))
  message(sprintf("[TALOS batch] Total duration: %.1f sec  |  %d variant(s) found",
    elapsed, nrow(combined)))
  invisible(combined)
}
#' Compute softclip-restricted paired-end support and orientation metrics
#'
#' Counts PE support only among QNAMEs with soft-clipped supporting reads, and
#' summarizes pair orientation for those paired QNAMEs.
#'
#' @keywords internal
compute_pe_softclip_metrics <- function(all_pairs, breakpoint, support_rows, best_len,
                                        flank = 5000L, min_mapq = 20L,
                                        max_insert = 5000L) {
  empty <- list(
    pe_softclip_support = 0L,
    pe_softclip_event_pairs = 0L,
    pe_softclip_long_pairs = 0L,
    pe_orientation_fr = 0L,
    pe_orientation_rf = 0L,
    pe_orientation_ff = 0L,
    pe_orientation_rr = 0L,
    pe_orientation_other = 0L,
    pe_orientation_dominant = NA_character_
  )
  if (is.null(all_pairs) || length(all_pairs) == 0L || is.null(support_rows) || nrow(support_rows) == 0L) return(empty)
  if (is.null(support_rows$cigar) || is.null(support_rows$read_name)) return(empty)

  has_sc <- grepl("^\\d+S", support_rows$cigar, perl = TRUE) |
            grepl("\\d+S$", support_rows$cigar, perl = TRUE)
  sc_qnames <- unique(support_rows$read_name[has_sc & !is.na(support_rows$read_name)])
  if (length(sc_qnames) == 0L) return(empty)

  r1 <- GenomicAlignments::first(all_pairs)
  r2 <- GenomicAlignments::last(all_pairs)
  pair_qnames <- names(all_pairs)
  if (is.null(pair_qnames) || length(pair_qnames) != length(all_pairs) || all(is.na(pair_qnames))) {
    pair_qnames <- as.character(S4Vectors::mcols(r1)$qname)
  }
  if (is.null(pair_qnames) || length(pair_qnames) != length(all_pairs)) return(empty)

  mq1 <- S4Vectors::mcols(r1)$mapq; if (is.null(mq1)) mq1 <- rep(0L, length(r1)); mq1[is.na(mq1)] <- 0L
  mq2 <- S4Vectors::mcols(r2)$mapq; if (is.null(mq2)) mq2 <- rep(0L, length(r2)); mq2[is.na(mq2)] <- 0L
  f1 <- S4Vectors::mcols(r1)$flag; if (is.null(f1)) f1 <- rep(0L, length(r1))
  f2 <- S4Vectors::mcols(r2)$flag; if (is.null(f2)) f2 <- rep(0L, length(r2))

  primary <- bitwAnd(f1, 0x100L) == 0L & bitwAnd(f1, 0x800L) == 0L &
             bitwAnd(f2, 0x100L) == 0L & bitwAnd(f2, 0x800L) == 0L
  valid <- primary & mq1 >= min_mapq & mq2 >= min_mapq

  s1 <- GenomicRanges::start(r1); e1 <- GenomicRanges::end(r1)
  s2 <- GenomicRanges::start(r2); e2 <- GenomicRanges::end(r2)
  pair_start <- pmin(s1, s2)
  pair_end <- pmax(e1, e2)
  outer_insert <- pair_end - pair_start + 1L
  size_ok <- outer_insert > 0L & outer_insert <= max_insert

  support_pair_mask <- valid & size_ok & pair_qnames %in% sc_qnames
  if (!any(support_pair_mask, na.rm = TRUE)) return(empty)

  span_bp <- pair_start <= breakpoint & pair_end >= breakpoint
  bg_mask <- valid & size_ok &
    pair_start >= (breakpoint - flank) & pair_end <= (breakpoint + flank) &
    !span_bp
  bg_med <- stats::median(outer_insert[bg_mask], na.rm = TRUE)
  bg_mad <- stats::mad(outer_insert[bg_mask], na.rm = TRUE)
  if (!is.finite(bg_mad) || bg_mad == 0) bg_mad <- stats::sd(outer_insert[bg_mask], na.rm = TRUE)
  if (!is.finite(bg_mad) || bg_mad == 0) bg_mad <- 50

  len_used <- if (is.na(best_len) || best_len <= 0L) 0L else as.numeric(best_len)
  event_size_span <- long_span <- rep(FALSE, length(all_pairs))
  if (is.finite(bg_med) && len_used > 0) {
    expected_insert <- bg_med + len_used
    insert_tol <- max(50, 2 * bg_mad)
    event_size_span <- span_bp & abs(outer_insert - expected_insert) <= insert_tol
    long_span <- span_bp & outer_insert >= (bg_med + 0.75 * len_used)
  }

  rev1 <- bitwAnd(f1, 0x10L) != 0L
  rev2 <- bitwAnd(f2, 0x10L) != 0L
  left_is_1 <- s1 <= s2
  left_rev <- ifelse(left_is_1, rev1, rev2)
  right_rev <- ifelse(left_is_1, rev2, rev1)
  orientation <- ifelse(!left_rev & right_rev, "FR",
                 ifelse(left_rev & !right_rev, "RF",
                 ifelse(!left_rev & !right_rev, "FF",
                 ifelse(left_rev & right_rev, "RR", "other"))))

  orient_tab <- table(factor(orientation[support_pair_mask], levels = c("FR", "RF", "FF", "RR", "other")))
  dominant <- if (sum(orient_tab) > 0L) names(sort(orient_tab, decreasing = TRUE))[1L] else NA_character_

  list(
    pe_softclip_support = length(unique(pair_qnames[support_pair_mask])),
    pe_softclip_event_pairs = as.integer(sum(support_pair_mask & event_size_span, na.rm = TRUE)),
    pe_softclip_long_pairs = as.integer(sum(support_pair_mask & long_span, na.rm = TRUE)),
    pe_orientation_fr = as.integer(orient_tab[["FR"]]),
    pe_orientation_rf = as.integer(orient_tab[["RF"]]),
    pe_orientation_ff = as.integer(orient_tab[["FF"]]),
    pe_orientation_rr = as.integer(orient_tab[["RR"]]),
    pe_orientation_other = as.integer(orient_tab[["other"]]),
    pe_orientation_dominant = dominant
  )
}

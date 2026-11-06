# ── Null-coalescing operator (used throughout the package) ───────────────────
`%||%` <- function(x, y) if (is.null(x)) y else x

# ── Shared run plumbing (used by detect_itd, detect_alu and the writers) ─────

#' Sample name derived from a BAM path (extension and everything after the
#' first dot removed), shared by detect_itd(), detect_alu() and the wrappers
#' so the three cannot disagree.
.sample_from_bam <- function(bam_path) sub("\\..*", "", basename(bam_path))

#' Timestamp used in output file names (one format everywhere).
.ts_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

#' Create a directory (recursively) and return it invisibly.
.ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

#' Flatten list columns to semicolon-joined strings so a data.frame can be
#' written as TSV (shared by .write_talos_output and .write_alu_output).
.flatten_list_cols <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- vapply(df[[col]], function(x) {
        if (length(x) == 0 || all(is.na(x))) NA_character_ else paste(x, collapse = ";")
      }, FUN.VALUE = character(1), USE.NAMES = FALSE)
    }
  }
  df
}

#' Write a results data.frame as TSV (list columns flattened, NA as ".").
.write_results_tsv <- function(df, path) {
  df <- .flatten_list_cols(df)
  utils::write.table(df, file = path, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = ".")
  invisible(path)
}

#' mcols column with a default, NA-filled (shared by .pair_info() and the
#' BAM evidence reader so both read record metadata the same way).
.mcol_or <- function(x, name, default) {
  v <- S4Vectors::mcols(x)[[name]]
  n <- if (is.null(x)) 0L else NROW(x)
  if (is.null(v)) v <- rep(default, n)
  v[is.na(v)] <- default
  v
}

# ── Shared helper for softclip extraction ───────────────────────────────────
.get_softclips <- function(cig, r_seq) {
  clips <- .get_softclip_lengths(cig)
  read_len <- if (is.na(r_seq)) 0L else nchar(r_seq)
  lead <- if (clips$lead > 0L && clips$lead <= read_len) {
    substr(r_seq, 1L, clips$lead)
  } else {
    NA_character_
  }
  trail <- if (clips$trail > 0L && clips$trail <= read_len) {
    substr(r_seq, read_len - clips$trail + 1L, read_len)
  } else {
    NA_character_
  }
  list(lead = lead, trail = trail)
}

.parse_cigar_ops <- function(cig) {
  if (length(cig) != 1L || is.na(cig) || !nzchar(cig)) {
    return(character())
  }
  regmatches(cig, gregexpr("[0-9]+[MIDNSHP=X]", cig, perl = TRUE))[[1L]]
}

.get_softclip_lengths <- function(cig) {
  empty <- list(lead = 0L, trail = 0L)
  ops <- .parse_cigar_ops(cig)
  if (length(ops) == 0L) {
    return(empty)
  }

  non_hard <- which(!grepl("H$", ops))
  if (length(non_hard) == 0L) {
    return(empty)
  }
  first_aligned <- non_hard[1L]
  last_aligned <- non_hard[length(non_hard)]
  lead <- if (grepl("S$", ops[first_aligned])) {
    as.integer(sub("S$", "", ops[first_aligned]))
  } else {
    0L
  }
  trail <- if (grepl("S$", ops[last_aligned])) {
    as.integer(sub("S$", "", ops[last_aligned]))
  } else {
    0L
  }
  if (first_aligned == last_aligned) {
    trail <- 0L
  }
  list(lead = lead, trail = trail)
}

.parse_cigar_insertions <- function(
  cig,
  read_start,
  genomic_start,
  ref_len,
  min_size = 1L
) {
  ops <- .parse_cigar_ops(cig)
  if (length(ops) == 0L) {
    return(data.frame())
  }

  ref_cursor <- as.integer(read_start)
  insertions <- list()
  n_insertions <- 0L
  for (op in ops) {
    op_len <- as.integer(sub("[A-Z=]$", "", op))
    op_type <- sub("^[0-9]+", "", op)
    if (op_type == "I") {
      bp <- ref_cursor - 1L
      local_bp <- bp - genomic_start + 1L
      if (op_len > 0L && op_len >= min_size && local_bp >= 1L && local_bp <= ref_len) {
        n_insertions <- n_insertions + 1L
        insertions[[n_insertions]] <- data.frame(
          local_breakpoint = as.integer(local_bp),
          length = op_len
        )
      }
    }
    if (op_type %in% c("M", "D", "N", "=", "X")) {
      ref_cursor <- ref_cursor + op_len
    }
  }

  if (n_insertions == 0L) {
    return(data.frame())
  }
  do.call(rbind, insertions)
}

.has_softclip <- function(cig) {
  clips <- lapply(cig, .get_softclip_lengths)
  vapply(clips, function(x) x$lead > 0L || x$trail > 0L, logical(1L))
}

.is_breakpoint_exonic <- function(position, chrom, exons) {
  if (is.null(exons) || length(exons) == 0L || is.na(position)) {
    return(FALSE)
  }
  breakpoint <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(start = position, end = position)
  )
  length(GenomicRanges::findOverlaps(breakpoint, exons, type = "within")) > 0L
}

.get_gene_exons <- function(gene_config) {
  gene_config$all_exons %||% gene_config$exons %||% gene_config$target_exons
}

# ── 3-letter amino acid lookup for HGVS protein notation ────────────────────
.aa1_to_aa3 <- c(
  A = "Ala",
  R = "Arg",
  N = "Asn",
  D = "Asp",
  C = "Cys",
  Q = "Gln",
  E = "Glu",
  G = "Gly",
  H = "His",
  I = "Ile",
  L = "Leu",
  K = "Lys",
  M = "Met",
  F = "Phe",
  P = "Pro",
  S = "Ser",
  T = "Thr",
  W = "Trp",
  Y = "Tyr",
  V = "Val",
  `*` = "Ter"
)

.to_three_letter <- function(aa1_str) {
  chars <- strsplit(aa1_str, "")[[1L]]
  paste(
    vapply(chars, function(ch) .aa1_to_aa3[[ch]] %||% ch, character(1L)),
    collapse = ""
  )
}


# ── Vectorised CIGAR helpers (base R; used by the split-read DUP parser) ─────
#' Parse CIGAR strings into reference width and end clips (soft + hard)
#' @param cig Character vector of CIGAR strings.
#' @return data.frame(width, lead, trail), one row per CIGAR. \code{width} is
#'   the span along the reference (M, D, N, =, X); \code{lead}/\code{trail}
#'   are the summed leading / trailing S and H operations.
#' @keywords internal
.cigar_layout <- function(cig) {
  one <- function(cg) {
    if (is.na(cg) || !nzchar(cg) || cg == "*") return(c(NA_integer_, 0L, 0L))
    ops <- regmatches(cg, gregexpr("[0-9]+[MIDNSHP=X]", cg, perl = TRUE))[[1L]]
    if (!length(ops)) return(c(NA_integer_, 0L, 0L))
    len <- as.integer(sub("[A-Z=]$", "", ops))
    op <- sub("^[0-9]+", "", ops)
    clip <- op %in% c("S", "H")
    if (all(clip)) return(c(NA_integer_, 0L, 0L))
    first_al <- which(!clip)[1L]
    last_al <- max(which(!clip))
    c(sum(len[op %in% c("M", "D", "N", "=", "X")]),
      sum(len[seq_len(first_al - 1L)]),
      if (last_al < length(ops)) sum(len[(last_al + 1L):length(ops)]) else 0L)
  }
  m <- vapply(as.character(cig), one, integer(3L), USE.NAMES = FALSE)
  if (!length(cig)) m <- matrix(integer(0), nrow = 3L)
  data.frame(width = m[1L, ], lead = m[2L, ], trail = m[3L, ])
}

#' Annotate ITDs with known hotspot information
#' @param itd_df Data frame with columns \code{Gene}, a position column and
#'   (optionally) \code{Length}.
#' @param db_path Optional path to CSV file (forces CSV mode).
#' @param genome_build \code{"hg19"} or \code{"hg38"}.
#' @param db_conn Optional DBI connection (overrides CSV).
#' @param fallback_csv_path Optional fallback CSV path.
#' @param pos_col Position column (default \code{"GenomicPosition"};
#'   \code{"InsertionSite"} for ALU results).
#' @return Data frame with added columns \code{Hotspot} and \code{HotspotName}.
#' @export
annotate_hotspots <- function(
  itd_df,
  db_path = NULL,
  genome_build = NULL,
  db_conn = NULL,
  fallback_csv_path = NULL,
  pos_col = "GenomicPosition"
) {
  itd_df$Hotspot <- rep(FALSE, nrow(itd_df))
  itd_df$HotspotName <- rep(NA_character_, nrow(itd_df))
  if (nrow(itd_df) == 0L) {
    return(itd_df)
  }
  if (!pos_col %in% names(itd_df)) {
    stop("annotate_hotspots(): position column '", pos_col, "' not found.")
  }

  hotspots <- NULL
  if (!is.null(db_conn)) {
    if (!requireNamespace("DBI", quietly = TRUE)) {
      stop("Package 'DBI' required for database mode.")
    }
    hs_list <- lapply(unique(itd_df$Gene), function(g) {
      build_clause <- if (!is.null(genome_build)) " AND Build = ?" else ""
      params <- if (!is.null(genome_build)) list(g, genome_build) else list(g)
      DBI::dbGetQuery(
        db_conn,
        paste0("SELECT Gene, Start, End, Name FROM hotspots WHERE Gene = ?", build_clause),
        params = params
      )
    })
    hotspots <- do.call(rbind, hs_list)
  } else {
    csv_path <- db_path %||% system.file("extdata", "hotspots.csv", package = "TALOS")
    if (nzchar(csv_path) && file.exists(csv_path)) {
      hotspots <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
      if (!is.null(genome_build) && "Build" %in% colnames(hotspots)) {
        hotspots <- hotspots[hotspots$Build == genome_build, , drop = FALSE]
      }
    } else if (!is.null(fallback_csv_path) && file.exists(fallback_csv_path)) {
      return(annotate_hotspots(itd_df, db_path = fallback_csv_path,
                               genome_build = genome_build, pos_col = pos_col))
    } else {
      message("No hotspot data available. Hotspot annotation skipped.")
      return(itd_df)
    }
  }
  if (is.null(hotspots) || nrow(hotspots) == 0L) {
    return(itd_df)
  }
  .overlap_hotspots(itd_df, hotspots, pos_col)
}

# Shared overlap step for the DB and CSV paths of annotate_hotspots().
.overlap_hotspots <- function(itd_df, hotspots, pos_col) {
  pos <- as.integer(itd_df[[pos_col]])
  len <- if ("Length" %in% names(itd_df)) as.integer(itd_df$Length) else rep(0L, nrow(itd_df))
  len[is.na(len)] <- 0L
  itd_end <- pmax(pos + len - 1L, pos)
  ok <- !is.na(pos) & !is.na(itd_df$Gene)
  if (!any(ok)) return(itd_df)
  hs_gr <- GenomicRanges::GRanges(
    seqnames = hotspots$Gene,
    ranges = IRanges::IRanges(start = hotspots$Start, end = hotspots$End)
  )
  itd_gr <- GenomicRanges::GRanges(
    seqnames = itd_df$Gene[ok],
    ranges = IRanges::IRanges(start = pos[ok], end = itd_end[ok])
  )
  hits <- suppressWarnings(GenomicRanges::findOverlaps(itd_gr, hs_gr, type = "any"))
  if (length(hits) == 0L) return(itd_df)
  rows <- which(ok)
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  for (qi in unique(qh)) {
    itd_df$Hotspot[rows[qi]] <- TRUE
    itd_df$HotspotName[rows[qi]] <- paste(unique(hotspots$Name[sh[qh == qi]]), collapse = ";")
  }
  itd_df
}

.annotate_exonic_region <- function(
  final_df,
  exons_gr,
  pos_col = "GenomicPosition"
) {
  if (is.null(exons_gr) || length(exons_gr) == 0L || nrow(final_df) == 0L) {
    final_df$Region <- rep(NA_character_, nrow(final_df))
    final_df$ExonNumber <- rep(NA_integer_, nrow(final_df))
    return(final_df)
  }

  if (!pos_col %in% names(final_df)) {
    stop("Position column not found in result data: ", pos_col)
  }

  exon_num_vec <- NULL
  if (!is.null(exons_gr$exon_number) && length(exons_gr$exon_number) > 0) {
    exon_num_vec <- as.integer(unlist(exons_gr$exon_number))
  } else if (!is.null(exons_gr$exon_idx) && length(exons_gr$exon_idx) > 0) {
    exon_num_vec <- as.integer(unlist(exons_gr$exon_idx))
  } else if (!is.null(names(exons_gr))) {
    exon_num_vec <- suppressWarnings(as.integer(names(exons_gr)))
  }

  bp_gr <- GenomicRanges::GRanges(
    seqnames = GenomeInfoDb::seqnames(exons_gr)[1L],
    ranges = IRanges::IRanges(
      start = final_df[[pos_col]],
      end = final_df[[pos_col]]
    )
  )
  hits <- GenomicRanges::findOverlaps(bp_gr, exons_gr, type = "within")
  final_df$Region <- "intronic"
  final_df$ExonNumber <- NA_integer_

  if (length(hits) > 0) {
    q_idx <- S4Vectors::queryHits(hits)
    s_idx <- S4Vectors::subjectHits(hits)
    keep <- !duplicated(q_idx)
    q_idx <- q_idx[keep]
    s_idx <- s_idx[keep]
    final_df$Region[q_idx] <- "exonic"
    if (!is.null(exon_num_vec)) {
      final_df$ExonNumber[q_idx] <- exon_num_vec[s_idx]
    }
  }
  final_df
}


#' HGVS c./p. annotation of a tandem duplication
#'
#' @param gene_config Resolved gene configuration.
#' @param genomic_pos Leftmost genomic base of the duplicated segment.
#' @param dup_seq Duplicated sequence (optional, used for its length).
#' @param dup_len Duplication length (bp).
#' @param debug Verbose debugging.
#' @return list(c_notation, p_notation)
#' @export
compute_hgvs_annotations <- function(
  gene_config,
  genomic_pos,
  dup_seq = NA_character_,
  dup_len = NA_integer_,
  debug = FALSE
) {
  genomic_pos <- suppressWarnings(as.numeric(genomic_pos))
  if (is.na(genomic_pos)) {
    return(list(c_notation = NA_character_, p_notation = NA_character_))
  }

  # ---- Check if transcript reference is available ----
  if (
    is.null(gene_config$transcript) ||
      is.na(gene_config$transcript) ||
      length(gene_config$transcript) == 0 ||
      gene_config$transcript == ""
  ) {
    return(list(c_notation = NA_character_, p_notation = NA_character_))
  }

  full_cdna <- gene_config$full_cdna
  if (is.null(full_cdna) || is.na(full_cdna) || nchar(full_cdna) == 0) {
    if (!is.null(gene_config$all_exons) && !is.null(gene_config$bsgenome_obj)) {
      if (debug) {
        message("[HGVS] Rebuilding full_cdna from exons and BSgenome")
      }
      exons_gr <- gene_config$all_exons
      strand <- gene_config$strand
      bs_obj <- gene_config$bsgenome_obj
      exon_seq_list <- lapply(seq_along(exons_gr), function(i) {
        seq <- as.character(BSgenome::getSeq(bs_obj, exons_gr[i]))
        if (strand == "-") {
          seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(
            seq
          )))
        }
        seq
      })
      full_cdna <- paste(unlist(exon_seq_list), collapse = "")
      if (nchar(full_cdna) == 0) {
        warning(
          "Could not rebuild full_cdna. cDNA annotation may be incomplete."
        )
        full_cdna <- NA_character_
      }
    } else {
      warning(
        "full_cdna missing and cannot be rebuilt. cDNA annotation may be incomplete."
      )
      full_cdna <- NA_character_
    }
  }

  seq_len <- if (!is.na(dup_seq) && nchar(dup_seq) > 0) {
    nchar(dup_seq)
  } else {
    NA_integer_
  }
  dup_len <- if (!is.na(dup_len) && dup_len > 0) {
    as.integer(dup_len)
  } else if (!is.na(seq_len)) {
    seq_len
  } else {
    NA_integer_
  }
  if (is.na(dup_len) || dup_len == 0L) {
    return(list(c_notation = NA_character_, p_notation = NA_character_))
  }

  exons_gr <- gene_config$all_exons
  if (is.null(exons_gr) || length(exons_gr) == 0) {
    return(list(c_notation = NA_character_, p_notation = NA_character_))
  }

  strand <- gene_config$strand

  if (!is.null(exons_gr$exon_idx)) {
    exons_gr <- exons_gr[order(exons_gr$exon_idx)]
  } else {
    if (strand == "+") {
      exons_gr <- exons_gr[order(GenomicRanges::start(exons_gr))]
    } else {
      exons_gr <- exons_gr[order(
        GenomicRanges::start(exons_gr),
        decreasing = TRUE
      )]
    }
  }

  exons_df <- data.frame(
    start = GenomicRanges::start(exons_gr),
    end = GenomicRanges::end(exons_gr),
    width = GenomicRanges::width(exons_gr),
    exon_idx = if (!is.null(exons_gr$exon_idx)) {
      exons_gr$exon_idx
    } else {
      seq_along(exons_gr)
    },
    stringsAsFactors = FALSE
  )
  exons_df$cum_len_before <- c(0L, cumsum(exons_df$width)[-nrow(exons_df)])

  cds_offset <- if (is.null(gene_config$cds_offset)) {
    0L
  } else {
    gene_config$cds_offset
  }

  # The duplicated segment is always [genomic_pos, genomic_pos + dup_len - 1]
  # on the reference (+) strand (GenomicPosition is its leftmost base). On a
  # minus-strand transcript the 5' end of the segment in cDNA space is its
  # rightmost genomic base, so c.start / c.end are taken from opposite ends.
  # (Previous versions used genomic_pos - dup_len + 1, i.e. the wrong
  # segment, for minus-strand genes such as FLT3.)
  seg_lo <- genomic_pos
  seg_hi <- genomic_pos + dup_len - 1L
  cdna_first <- if (strand == "-") seg_hi else seg_lo
  cdna_last <- if (strand == "-") seg_lo else seg_hi

  .genomic_to_cdna_hgvs <- function(pos, exons_df, strand, cds_offset) {
    exon_hit <- which(exons_df$start <= pos & exons_df$end >= pos)
    if (length(exon_hit) > 0) {
      idx <- exon_hit[1L]
      pos_in_exon <- if (strand == "+") {
        pos - exons_df$start[idx] + 1L
      } else {
        exons_df$end[idx] - pos + 1L
      }
      c_pos <- exons_df$cum_len_before[idx] + pos_in_exon - cds_offset
      # 5' UTR: there is no c.0, the base before c.1 is c.-1.
      if (c_pos <= 0L) {
        return(as.character(c_pos - 1L))
      }
      return(as.character(c_pos))
    }

    if (strand == "+") {
      up_mask <- exons_df$end < pos
      dn_mask <- exons_df$start > pos
    } else {
      up_mask <- exons_df$start > pos
      dn_mask <- exons_df$end < pos
    }

    up_note <- NULL
    up_off <- Inf
    dn_note <- NULL
    dn_off <- Inf

    if (any(up_mask)) {
      up_idx <- tail(which(up_mask), 1L)
      up_c <- exons_df$cum_len_before[up_idx] +
        exons_df$width[up_idx] -
        cds_offset
      up_off <- if (strand == "+") {
        pos - exons_df$end[up_idx]
      } else {
        exons_df$start[up_idx] - pos
      }
      up_note <- paste0(up_c, "+", up_off)
    }
    if (any(dn_mask)) {
      dn_idx <- head(which(dn_mask), 1L)
      dn_c <- exons_df$cum_len_before[dn_idx] + 1L - cds_offset
      dn_off <- if (strand == "+") {
        exons_df$start[dn_idx] - pos
      } else {
        pos - exons_df$end[dn_idx]
      }
      dn_note <- paste0(dn_c, "-", dn_off)
    }

    if (is.null(up_note) && is.null(dn_note)) {
      return(paste0("?+", pos))
    }
    if (is.null(up_note)) {
      return(dn_note)
    }
    if (is.null(dn_note)) {
      return(up_note)
    }
    if (up_off <= dn_off) up_note else dn_note
  }

  start_notation <- .genomic_to_cdna_hgvs(
    cdna_first,
    exons_df,
    strand,
    cds_offset
  )
  end_notation <- .genomic_to_cdna_hgvs(
    cdna_last,
    exons_df,
    strand,
    cds_offset
  )

  if (is.na(start_notation) || is.na(end_notation)) {
    return(list(c_notation = "c.?+?", p_notation = "intronic"))
  }

  c_notation <- paste0("c.", start_notation, "_", end_notation, "dup")

  p_notation <- NA_character_
  if (
    !is.na(c_notation) &&
      !is.na(dup_len) &&
      !is.na(full_cdna) &&
      nchar(full_cdna) > cds_offset + 3L
  ) {
    coding_seq <- substr(full_cdna, cds_offset + 1L, nchar(full_cdna))
    start_c <- suppressWarnings(as.numeric(start_notation))
    end_c <- suppressWarnings(as.numeric(end_notation))

    if (!is.na(start_c) && !is.na(end_c) && start_c >= 1L && end_c >= 1L) {
      c_start <- as.integer(start_c)
      c_end <- as.integer(end_c)

      if (c_start >= 1L && c_end <= nchar(coding_seq)) {
        dup_cds <- substr(coding_seq, c_start, c_end)

        if (nchar(dup_cds) >= 3L) {
          aa_start_pos <- ceiling(c_start / 3L)
          aa_end_pos <- ceiling(c_end / 3L)

          cds_trim <- substr(
            coding_seq,
            1L,
            nchar(coding_seq) - nchar(coding_seq) %% 3L
          )
          if (nchar(cds_trim) >= 3L) {
            full_aa <- suppressWarnings(
              as.character(Biostrings::translate(Biostrings::DNAString(
                cds_trim
              )))
            )
            if (aa_end_pos <= nchar(full_aa)) {
              aa1_1 <- substr(full_aa, aa_start_pos, aa_start_pos)
              aa1_2 <- substr(full_aa, aa_end_pos, aa_end_pos)
              aa3_1 <- .to_three_letter(aa1_1)
              aa3_2 <- .to_three_letter(aa1_2)

              if (nchar(dup_cds) %% 3L != 0L) {
                p_notation <- paste0("p.", aa3_1, aa_start_pos, "fsTer?")
              } else if (aa_start_pos == aa_end_pos) {
                p_notation <- paste0("p.", aa3_1, aa_start_pos, "dup")
              } else {
                p_notation <- paste0(
                  "p.",
                  aa3_1,
                  aa_start_pos,
                  "_",
                  aa3_2,
                  aa_end_pos,
                  "dup"
                )
              }
            }
          }
        }
      }
    }
  }

  # "123+45" / "124-30" are intronic offsets; a leading "-" is the 5' UTR.
  intronic_pat <- "[0-9][+-][0-9]"
  if (is.na(p_notation)) {
    if (grepl(intronic_pat, start_notation) || grepl(intronic_pat, end_notation)) {
      p_notation <- "intronic"
    } else if (startsWith(start_notation, "-")) {
      p_notation <- "5'UTR"
    }
  }

  list(c_notation = c_notation, p_notation = p_notation)
}



.get_coverage_rle <- function(cov_all, chrom) {
  if (is.null(cov_all)) {
    return(NULL)
  }
  if (inherits(cov_all, "Rle") || is.numeric(cov_all)) {
    return(cov_all)
  }
  if (is.list(cov_all) && length(cov_all) == 1L && is.null(names(cov_all))) {
    return(cov_all[[1L]])
  }
  cov_names <- names(cov_all)
  if (!is.null(cov_names) && chrom %in% cov_names) {
    return(cov_all[[chrom]])
  }
  NULL
}

compute_coverage_drop <- function(
  cov_all,
  chrom,
  breakpoint,
  flank = 200,
  min_depth = 5,
  min_flank_bp = 5,
  verbose = FALSE
) {
  if (is.null(cov_all)) {
    if (verbose) {
      message("Coverage drop: cov_all is NULL")
    }
    return(NA_real_)
  }

  cov_rle <- .get_coverage_rle(cov_all, chrom)
  if (is.null(cov_rle)) {
    if (verbose) {
      message("Coverage drop: chrom missing in cov_all list")
    }
    return(NA_real_)
  }

  cov_len <- length(cov_rle)
  start_pos <- max(1L, breakpoint - flank)
  end_pos <- min(cov_len, breakpoint + flank)
  if (start_pos >= end_pos) {
    if (verbose) {
      message("Coverage drop: start_pos >= end_pos")
    }
    return(NA_real_)
  }

  cov_vec <- as.numeric(cov_rle[start_pos:end_pos])
  bp_relative <- breakpoint - start_pos + 1L
  left_end_rel <- max(1L, bp_relative - 10L)
  right_start_rel <- min(length(cov_vec), bp_relative + 10L)

  if (
    left_end_rel < min_flank_bp ||
      (length(cov_vec) - right_start_rel) < min_flank_bp
  ) {
    if (verbose) {
      message(
        "Coverage drop: insufficient flanking bases (left=",
        left_end_rel,
        ", right=",
        length(cov_vec) - right_start_rel,
        ")"
      )
    }
    return(NA_real_)
  }

  left_cov <- mean(cov_vec[1L:left_end_rel], na.rm = TRUE)
  right_cov <- mean(cov_vec[right_start_rel:length(cov_vec)], na.rm = TRUE)

  if (left_cov < min_depth || right_cov < min_depth) {
    if (verbose) {
      message(
        "Coverage drop: left_cov=",
        left_cov,
        " right_cov=",
        right_cov,
        " below min_depth=",
        min_depth
      )
    }
    return(NA_real_)
  }

  ratio <- max(left_cov, right_cov) / min(left_cov, right_cov)
  if (verbose) {
    message("Coverage drop ratio = ", ratio)
  }
  ratio
}

compute_local_coverage <- function(
  cov_all,
  chrom,
  breakpoint,
  buffer = 10L,
  verbose = FALSE
) {
  if (is.null(cov_all)) {
    if (verbose) {
      message("Local coverage: cov_all is NULL")
    }
    return(NA_real_)
  }

  cov_rle <- .get_coverage_rle(cov_all, chrom)
  if (is.null(cov_rle)) {
    if (verbose) {
      message("Local coverage: chrom missing in cov_all list")
    }
    return(NA_real_)
  }

  cov_len <- length(cov_rle)
  start_pos <- breakpoint - buffer
  end_pos <- breakpoint + buffer
  if (start_pos < 1L || end_pos > cov_len || start_pos > end_pos) {
    if (verbose) {
      message("Local coverage: window out of range")
    }
    return(NA_real_)
  }

  min_cov <- min(as.numeric(cov_rle[start_pos:end_pos]))
  if (verbose) {
    message("Local coverage (min over bp+/-", buffer, ") = ", min_cov)
  }
  min_cov
}

compute_span_coverage <- function(
  cov_all,
  chrom,
  start,
  span_len,
  verbose = FALSE
) {
  summary <- .compute_span_coverage_summary(cov_all, chrom, start, span_len)
  if (verbose && !is.na(summary$mean_cov)) {
    message(
      "Span coverage over ",
      span_len,
      "bp: min=",
      summary$min_cov,
      " mean=",
      round(summary$mean_cov, 1)
    )
  }
  summary[c("min_cov", "mean_cov")]
}

.compute_span_coverage_summary <- function(
  cov_all,
  chrom,
  start,
  span_len,
  flank = 300L
) {
  empty_result <- list(
    min_cov = NA_real_,
    mean_cov = NA_real_,
    depth_fold_change = NA_real_
  )
  if (is.null(cov_all) || is.na(start) || is.na(span_len) || span_len < 1L) {
    return(empty_result)
  }

  cov_rle <- .get_coverage_rle(cov_all, chrom)
  if (is.null(cov_rle)) {
    return(empty_result)
  }

  cov_len <- length(cov_rle)
  start_pos <- start
  end_pos <- start + span_len - 1L
  if (start_pos < 1L || end_pos > cov_len || start_pos > end_pos) {
    return(empty_result)
  }

  cov_vec <- as.numeric(cov_rle[start_pos:end_pos])
  span_mean <- mean(cov_vec)
  left_start <- max(1L, start_pos - flank)
  left_end <- start_pos - 1L
  right_start <- end_pos + 1L
  right_end <- min(cov_len, end_pos + flank)
  flank_vals <- c(
    if (left_end >= left_start) {
      as.numeric(cov_rle[left_start:left_end])
    } else {
      numeric()
    },
    if (right_end >= right_start) {
      as.numeric(cov_rle[right_start:right_end])
    } else {
      numeric()
    }
  )
  background_mean <- if (length(flank_vals)) mean(flank_vals) else NA_real_
  depth_fold_change <- if (is.finite(background_mean) && background_mean > 0) {
    span_mean / background_mean
  } else {
    NA_real_
  }

  list(
    min_cov = min(cov_vec),
    mean_cov = span_mean,
    depth_fold_change = depth_fold_change
  )
}

compute_span_depth_fold_change <- function(
  cov_all,
  chrom,
  start,
  span_len,
  flank = 300L,
  verbose = FALSE
) {
  if (is.null(cov_all) || is.na(start) || is.na(span_len) || span_len < 1L) {
    return(NA_real_)
  }

  .compute_span_coverage_summary(
    cov_all,
    chrom,
    start,
    span_len,
    flank = flank
  )$depth_fold_change
}



#' Shared per-pair summary for the paired-end metrics (DRY helper)
#'
#' Extracts first/last mates, MAPQ, flags, coordinates, primary/MAPQ validity
#' and pair orientation once, so compute_discordant_ratio(),
#' compute_pe_itd_length() and compute_pe_softclip_metrics() read pairs the
#' same way.
#' @param all_pairs A \code{GAlignmentPairs} object.
#' @param min_mapq Minimum MAPQ for both mates.
#' @return list of per-pair vectors, or NULL when there are no pairs.
#' @keywords internal
.pair_info <- function(all_pairs, min_mapq = 20L) {
  if (is.null(all_pairs) || length(all_pairs) == 0L) return(NULL)
  r1 <- GenomicAlignments::first(all_pairs)
  r2 <- GenomicAlignments::last(all_pairs)
  n <- length(r1)
  mq1 <- .mcol_or(r1, "mapq", 0L); mq2 <- .mcol_or(r2, "mapq", 0L)
  f1 <- .mcol_or(r1, "flag", 0L);  f2 <- .mcol_or(r2, "flag", 0L)
  primary <- bitwAnd(f1, 0x100L) == 0L & bitwAnd(f1, 0x800L) == 0L &
    bitwAnd(f2, 0x100L) == 0L & bitwAnd(f2, 0x800L) == 0L
  s1 <- GenomicRanges::start(r1); e1 <- GenomicRanges::end(r1)
  s2 <- GenomicRanges::start(r2); e2 <- GenomicRanges::end(r2)
  rev1 <- bitwAnd(f1, 0x10L) != 0L
  rev2 <- bitwAnd(f2, 0x10L) != 0L
  left_is_1 <- s1 <= s2
  left_rev <- ifelse(left_is_1, rev1, rev2)
  right_rev <- ifelse(left_is_1, rev2, rev1)
  orientation <- ifelse(!left_rev & right_rev, "FR",
                 ifelse(left_rev & !right_rev, "RF",
                 ifelse(!left_rev & !right_rev, "FF", "RR")))
  qn <- names(all_pairs)
  if (is.null(qn) || length(qn) != n || all(is.na(qn))) {
    qn <- as.character(S4Vectors::mcols(r1)$qname)
  }
  pair_start <- pmin(s1, s2)
  pair_end <- pmax(e1, e2)
  list(
    n = n, qname = qn, s1 = s1, e1 = e1, s2 = s2, e2 = e2,
    rev1 = rev1, rev2 = rev2, left_rev = left_rev, right_rev = right_rev,
    orientation = orientation, primary = primary,
    valid = primary & mq1 >= min_mapq & mq2 >= min_mapq,
    pair_start = pair_start, pair_end = pair_end,
    outer_insert = pair_end - pair_start + 1L
  )
}

compute_discordant_ratio <- function(
  all_pairs,
  breakpoint,
  flank = 5000,
  min_mapq = 20,
  insert_size_factor = 2.0,
  min_everted_separation = 100L
) {
  pi <- .pair_info(all_pairs, min_mapq)
  if (is.null(pi)) return(NA_real_)
  in_window <- pi$s1 >= max(1L, breakpoint - flank) & pi$e1 <= breakpoint + flank
  keep <- in_window & pi$valid
  if (!any(in_window) || sum(keep) == 0L) return(NA_real_)

  start1 <- pi$s1[keep]; start2 <- pi$s2[keep]
  mapped_isize <- pmax(pi$e1[keep], pi$e2[keep]) - pmin(start1, start2)
  med_insert <- stats::median(mapped_isize, na.rm = TRUE)
  sd_insert <- stats::sd(mapped_isize, na.rm = TRUE)
  if (is.na(sd_insert)) sd_insert <- 0
  large_insert <- mapped_isize > (med_insert + insert_size_factor * sd_insert)
  spanning <- (start1 <= breakpoint & start2 > breakpoint) |
    (start2 <= breakpoint & start1 > breakpoint)
  everted <- pi$left_rev[keep] & !pi$right_rev[keep] &
    abs(start1 - start2) >= min_everted_separation
  sum(spanning & (large_insert | everted)) / sum(keep)
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
compute_pe_itd_length <- function(
  all_pairs,
  breakpoint,
  nominal_read_len = 150L,
  flank = 5000L,
  min_mapq = 20L,
  min_spanning = 5L,
  max_insert = 5000L
) {
  empty <- list(length_pe = NA_real_, n_spanning = 0L, n_background = 0L)
  pi <- .pair_info(all_pairs, min_mapq)
  if (is.null(pi)) return(empty)
  isize <- pi$outer_insert
  size_mask <- isize > 0L & isize <= max_insert
  bg_mask <- pi$valid & size_mask &
    pi$pair_start >= (breakpoint - flank) & pi$pair_end <= (breakpoint + flank) &
    !(pi$pair_start < breakpoint & pi$pair_end > breakpoint)
  span_mask <- pi$valid & size_mask &
    pi$pair_start <= (breakpoint - nominal_read_len) &
    pi$pair_end >= (breakpoint + nominal_read_len)
  n_span <- sum(span_mask, na.rm = TRUE)
  n_bg <- sum(bg_mask, na.rm = TRUE)
  if (n_span < min_spanning || n_bg < 10L) return(empty)
  est <- round(stats::median(isize[span_mask], na.rm = TRUE) -
                 stats::median(isize[bg_mask], na.rm = TRUE))
  if (is.na(est) || est <= 0L) return(empty)
  list(length_pe = est, n_spanning = n_span, n_background = n_bg)
}

compute_microhomology <- function(
  support_rows,
  ref_dna,
  breakpoint,
  genomic_start,
  debug = FALSE
) {
  if (nrow(support_rows) == 0) {
    return(NA_integer_)
  }
  mh_lengths <- integer()
  local_bp <- breakpoint - genomic_start + 1L
  ref_len_total <- nchar(ref_dna)

  has_biostrings <- .has_biostrings

  for (i in seq_len(nrow(support_rows))) {
    cig <- support_rows$cigar[i]
    read_seq <- support_rows$read_seq[i]

    sc <- .get_softclips(cig, read_seq)

    if (!is.na(sc$lead)) {
      len <- nchar(sc$lead)
      ref_start <- local_bp + 1L
      ref_end <- min(ref_len_total, local_bp + len)
      ref_down <- substr(ref_dna, ref_start, ref_end)
      if (nchar(ref_down) >= len) {
        mm <- if (has_biostrings) {
          Biostrings::lcprefix(sc$lead, ref_down)
        } else {
          min_len <- min(nchar(sc$lead), nchar(ref_down))
          if (substr(sc$lead, 1, min_len) == substr(ref_down, 1, min_len)) {
            min_len
          } else {
            0
          }
        }
        if (mm > 0) mh_lengths <- c(mh_lengths, mm)
      }
    }

    if (!is.na(sc$trail)) {
      len <- nchar(sc$trail)
      ref_start <- max(1L, local_bp - len)
      ref_end <- local_bp - 1L
      ref_up <- substr(ref_dna, ref_start, ref_end)
      if (nchar(ref_up) >= len) {
        mm <- if (has_biostrings) {
          Biostrings::lcprefix(sc$trail, ref_up)
        } else {
          min_len <- min(nchar(sc$trail), nchar(ref_up))
          if (substr(sc$trail, 1, min_len) == substr(ref_up, 1, min_len)) {
            min_len
          } else {
            0
          }
        }
        if (mm > 0) mh_lengths <- c(mh_lengths, mm)
      }
    }
  }
  if (length(mh_lengths) == 0) {
    return(NA_integer_)
  }
  as.integer(median(mh_lengths))
}

compute_repeat_entropy <- function(
  ref_dna,
  breakpoint,
  genomic_start,
  window = 50
) {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    warning("Biostrings not available; repeat entropy set to 0")
    return(0)
  }
  local_pos <- breakpoint - genomic_start + 1L
  start_pos <- max(1L, local_pos - window)
  end_pos <- min(nchar(ref_dna), local_pos + window)
  seq_window <- substr(ref_dna, start_pos, end_pos)
  if (nchar(seq_window) < 10L) {
    return(0)
  }
  freq <- Biostrings::oligonucleotideFrequency(
    Biostrings::DNAString(seq_window),
    2L
  )
  freq <- freq / sum(freq)
  -sum(freq * log2(freq + 1e-12))
}

#' Detect orientation of a duplication (tandem vs inverted) using alignment
#' @param itd_seq Duplicated sequence (character)
#' @param ref_seq Reference segment (character)
#' @param min_pid Minimum percent identity to call (0-1)
#' @export
detect_orientation <- function(itd_seq, ref_seq, min_pid = 0.90) {
  if (
    is.na(itd_seq) ||
      is.na(ref_seq) ||
      nchar(itd_seq) == 0 ||
      nchar(ref_seq) == 0
  ) {
    return(NA_character_)
  }
  itd_seq <- toupper(itd_seq)
  ref_seq <- toupper(ref_seq)

  has_pwalign <- requireNamespace("pwalign", quietly = TRUE)
  has_biostrings <- requireNamespace("Biostrings", quietly = TRUE)

  if (nchar(itd_seq) < 10 || nchar(ref_seq) < 10) {
    if (itd_seq == ref_seq) {
      return("tandem")
    }
    if (has_biostrings) {
      rc <- as.character(Biostrings::reverseComplement(Biostrings::DNAString(
        ref_seq
      )))
      if (itd_seq == rc) return("inverted")
    }
    return(NA_character_)
  }

  if (has_pwalign && has_biostrings) {
    aln_fwd <- pwalign::pairwiseAlignment(itd_seq, ref_seq, type = "local")
    if (pwalign::pid(aln_fwd) >= min_pid * 100) {
      return("tandem")
    }
    rc <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAString(ref_seq))
    )
    aln_rev <- pwalign::pairwiseAlignment(itd_seq, rc, type = "local")
    if (pwalign::pid(aln_rev) >= min_pid * 100) return("inverted")
  } else if (has_biostrings) {
    if (itd_seq == ref_seq) {
      return("tandem")
    }
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
  if (!is.na(local_pos) && local_pos > 0L && local_pos <= nchar(ref_dna)) {
    substr(ref_dna, local_pos, local_pos)
  } else {
    "N"
  }
}

# VCF INFO values must not contain ';', '=', whitespace or commas-in-names.
.vcf_safe <- function(x) gsub("[;=[:space:]]", "_", as.character(x))

# Optional INFO fields: column, key, sprintf format, coercion.
.VCF_OPTIONAL_INFO <- list(
  list("StrandBias", "SB", "%.4f", as.numeric),
  list("MeanSupportMAPQ", "MMAPQ", "%.1f", as.numeric),
  list("BreakpointSpread", "BPSPREAD", "%d", as.integer),
  list("SoftclipFraction", "SCFRAC", "%.4f", as.numeric),
  list("UniqueBreakpoints", "UBP", "%d", as.integer),
  list("CoverageDrop", "CD", "%.2f", as.numeric),
  list("DiscordantRatio", "DRATIO", "%.4f", as.numeric),
  list("MedianMicrohomology", "MH", "%d", as.integer),
  list("RepeatEntropy", "ENT", "%.2f", as.numeric),
  list("Orientation", "ORIENT", "%s", .vcf_safe),
  list("SupportConsistency", "CONSIST", "%.4f", as.numeric),
  list("RefMatch_Observed", "REFMATCH_OBS", "%.1f", as.numeric),
  list("RefMatch_Total", "REFMATCH_TOT", "%.1f", as.numeric),
  list("ITDReadCoverage", "ITDCOV", "%.1f", as.numeric),
  list("HotspotName", "HOTSPOT", "%s", function(v) .vcf_safe(gsub(";[[:space:]]*", ",", v))),
  list("AlignmentScore", "ALNSCORE", "%.4f", as.numeric),
  list("TotalSupportBases", "SUPBASES", "%d", as.integer),
  list("LeftSoftclipCount", "LSC", "%d", as.integer),
  list("RightSoftclipCount", "RSC", "%d", as.integer),
  list("LeftSoftclipPctSupport", "LSCPCT", "%.1f", as.numeric),
  list("RightSoftclipPctSupport", "RSCPCT", "%.1f", as.numeric)
)

.format_info <- function(row) {
  len <- suppressWarnings(as.integer(row$Length))
  pos <- as.integer(row$GenomicPosition)
  end <- if (!is.na(len) && len > 0L) pos + len - 1L else pos
  cds_val <- if (is.na(row$HGVS_cDNA)) "." else .vcf_safe(row$HGVS_cDNA)
  info <- sprintf(
    "SVTYPE=DUP;SVLEN=%s;END=%d;GENE=%s;CDS=%s;AF=%.4f;DP=%d;SUPPORT=%d;WT=%d",
    if (is.na(len)) "." else as.character(len), end, .vcf_safe(row$Gene), cds_val,
    as.numeric(row$AlleleFrequency), as.integer(row$DepthAtBreakpoint),
    as.integer(row$SupportingReads), as.integer(row$WildtypeReads)
  )
  if (!is.na(row$HGVS_Protein) && row$HGVS_Protein != "") {
    info <- paste0(info, ";AA=", .vcf_safe(row$HGVS_Protein))
  }
  for (spec in .VCF_OPTIONAL_INFO) {
    col <- spec[[1L]]
    v <- row[[col]]
    if (is.null(v) || length(v) == 0L || is.na(v)) next
    if (col == "HotspotName" && !isTRUE(row$Hotspot)) next
    info <- paste0(info, ";", spec[[2L]], "=", sprintf(spec[[3L]], spec[[4L]](v)))
  }
  info
}

#' Build VCF header lines
#' @export
build_vcf_header <- function(
  sample_name,
  chrom,
  genome_build = NULL,
  contig_length = NULL
) {
  pkg_ver <- tryCatch(
    as.character(utils::packageVersion("TALOS")),
    error = function(e) "dev"
  )
  lines <- c(
    "##fileformat=VCFv4.2",
    paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
    paste0("##source=TALOS_v", pkg_ver),
    paste0("##SAMPLE=<ID=", sample_name, ">")
  )
  if (!is.null(genome_build)) {
    lines <- c(lines, paste0("##reference=", genome_build))
  }

  if (!is.null(chrom)) {
    if (
      !is.null(contig_length) && is.numeric(contig_length) && contig_length > 0
    ) {
      lines <- c(
        lines,
        paste0("##contig=<ID=", chrom, ",length=", contig_length, ">")
      )
    } else {
      lines <- c(lines, paste0("##contig=<ID=", chrom, ",length=?>"))
    }
  }

  lines <- c(
    lines,
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
    paste(
      c(
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
        sample_name
      ),
      collapse = "\t"
    )
  )
  lines
}

#' Build VCF data lines from TALOS result data frame
#' @export
build_vcf_records <- function(itd_df, ref_dna, genomic_start, chrom) {
  if (nrow(itd_df) == 0L) {
    return(character(0L))
  }
  records <- character(nrow(itd_df))
  for (i in seq_len(nrow(itd_df))) {
    row <- itd_df[i, ]
    # A missing position or depth would emit "NA" fields and produce an
    # invalid VCF row; such rows are skipped instead.
    if (!isTRUE(row$DepthAtBreakpoint > 0)) next
    if (is.na(row$GenomicPosition)) next
    pos <- row$GenomicPosition
    id_field <- if (!is.na(row$HGVS_cDNA)) .vcf_safe(row$HGVS_cDNA) else "."
    records[i] <- paste(
      chrom, pos, id_field, .resolve_ref_base(pos, ref_dna, genomic_start),
      "<DUP>", ".", "PASS", .format_info(row), "GT", "0/1",
      sep = "\t"
    )
  }
  records[nzchar(records)]
}

#' Write a VCF file from TALOS results
#' @export
write_itd_vcf <- function(
  itd_df,
  ref_dna,
  genomic_start,
  chrom,
  sample_name,
  genome_build = NULL,
  vcf_path,
  overwrite = FALSE
) {
  if (!overwrite && file.exists(vcf_path)) {
    stop("VCF already exists: ", vcf_path)
  }
  out_dir <- dirname(vcf_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  contig_length <- NULL
  if (!is.null(genome_build) && requireNamespace("BSgenome", quietly = TRUE)) {
    bs_pkg <- if (genome_build == "hg19") {
      "BSgenome.Hsapiens.UCSC.hg19"
    } else {
      "BSgenome.Hsapiens.UCSC.hg38"
    }
    if (requireNamespace(bs_pkg, quietly = TRUE)) {
      bs_obj <- BSgenome::getBSgenome(bs_pkg)
      if (chrom %in% names(bs_obj)) {
        contig_length <- length(bs_obj[[chrom]])
      }
    }
  }

  lines <- c(
    build_vcf_header(sample_name, chrom, genome_build, contig_length),
    build_vcf_records(itd_df, ref_dna, genomic_start, chrom)
  )
  writeLines(lines, vcf_path)
  invisible(vcf_path)
}


.infer_sample_names <- function(bam_paths, part = 1L, sep = NULL) {
  basenames <- sub("\\.bam$", "", basename(bam_paths), ignore.case = TRUE)

  vapply(
    basenames,
    function(bn) {
      if (is.null(sep)) {
        sep_used <- if (grepl("\\.", bn)) {
          "\\."
        } else if (grepl("_", bn)) {
          "_"
        } else if (grepl("-", bn)) {
          "-"
        } else {
          NULL
        }
      } else {
        sep_used <- sep
      }
      if (is.null(sep_used)) {
        return(bn)
      }
      parts <- strsplit(bn, sep_used)[[1L]]
      n_parts <- length(parts)
      if (part > n_parts) {
        warning(sprintf(
          "BAM '%s' has only %d part(s); using full name.",
          bn,
          n_parts
        ))
        return(bn)
      }
      parts[part]
    },
    character(1L),
    USE.NAMES = FALSE
  )
}

#' Batch TALOS: analyse every BAM in a directory
#' @export
talos_batch_dir <- function(
  bam_dir,
  genes,
  build = "hg19",
  output_folder = "./results",
  pattern = "\\.bam$",
  sample_name_part = 1L,
  sample_name_sep = NULL,
  recursive = FALSE,
  combine_output = TRUE,
  ...
) {
  if (!dir.exists(bam_dir)) {
    stop("Directory not found: ", bam_dir)
  }
  bam_files <- list.files(
    bam_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = recursive
  )
  if (length(bam_files) == 0) {
    stop(
      "No BAM files found in '",
      bam_dir,
      "' matching pattern '",
      pattern,
      "'"
    )
  }

  sample_names <- .infer_sample_names(
    bam_files,
    part = sample_name_part,
    sep = sample_name_sep
  )
  message(sprintf(
    "[TALOS batch_dir] Found %d BAM file(s) in '%s'.",
    length(bam_files),
    bam_dir
  ))
  for (i in seq_along(bam_files)) {
    message(sprintf(
      "  %s  ->  sample: %s",
      basename(bam_files[i]),
      sample_names[i]
    ))
  }

  talos_batch(
    bam_paths = bam_files,
    genes = genes,
    build = build,
    sample_names = sample_names,
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
  build = "hg19",
  sample_names = NULL,
  output_folder = "./results",
  combine_output = TRUE,
  ...
) {
  .parse_input <- function(x) {
    if (length(x) == 1L && grepl(",", x)) {
      trimws(strsplit(x, ",")[[1L]])
    } else {
      as.character(x)
    }
  }

  bam_paths <- .parse_input(bam_paths)
  genes <- .parse_input(genes)
  sample_names <- if (!is.null(sample_names)) {
    .parse_input(sample_names)
  } else {
    .infer_sample_names(bam_paths)
  }

  if (length(sample_names) != length(bam_paths)) {
    stop("talos_batch: length of sample_names must equal length of bam_paths.")
  }

  batch_start <- Sys.time()
  message(sprintf(
    "[TALOS batch] %d sample(s) x %d gene(s) = %d run(s)",
    length(bam_paths),
    length(genes),
    length(bam_paths) * length(genes)
  ))

  all_results <- list()
  for (bi in seq_along(bam_paths)) {
    for (gi in seq_along(genes)) {
      bam <- bam_paths[bi]
      gene <- genes[gi]
      sn <- sample_names[bi]
      message(sprintf("[TALOS batch] Running: sample=%s  gene=%s", sn, gene))
      res <- tryCatch(
        talos(
          bam_path = bam,
          gene = gene,
          build = build,
          sample_name = sn,
          output_folder = output_folder,
          ...
        ),
        error = function(e) {
          message(sprintf(
            "[TALOS batch] ERROR for %s / %s: %s",
            sn,
            gene,
            conditionMessage(e)
          ))
          data.frame()
        }
      )
      if (is.data.frame(res) && nrow(res) > 0) {
        all_results[[length(all_results) + 1L]] <- res
      }
    }
  }

  combined <- if (length(all_results) > 0) {
    do.call(rbind, all_results)
  } else {
    data.frame()
  }

  if (combine_output && nrow(combined) > 0) {
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
    }
    ts <- .ts_stamp()
    cfile <- file.path(output_folder, paste0("TALOS_batch_", ts, ".tsv"))
    write.table(
      combined,
      file = cfile,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      na = "."
    )
    message(sprintf("[TALOS batch] Combined results written to: %s", cfile))
  }

  elapsed <- as.numeric(difftime(Sys.time(), batch_start, units = "secs"))
  message(sprintf(
    "[TALOS batch] Total duration: %.1f sec  |  %d variant(s) found",
    elapsed,
    nrow(combined)
  ))
  invisible(combined)
}
#' Compute softclip-restricted paired-end support and orientation metrics
#'
#' Counts PE support only among QNAMEs with soft-clipped supporting reads, and
#' summarizes pair orientation for those paired QNAMEs.
#'
#' @keywords internal
compute_pe_softclip_metrics <- function(
  all_pairs,
  breakpoint,
  support_rows,
  best_len,
  flank = 5000L,
  min_mapq = 20L,
  max_insert = 5000L
) {
  orient_levels <- c("FR", "RF", "FF", "RR", "other")
  empty <- list(
    pe_softclip_support = 0L, pe_softclip_event_pairs = 0L, pe_softclip_long_pairs = 0L,
    pe_orientation_fr = 0L, pe_orientation_rf = 0L, pe_orientation_ff = 0L,
    pe_orientation_rr = 0L, pe_orientation_other = 0L,
    pe_orientation_dominant = NA_character_
  )
  if (is.null(support_rows) || nrow(support_rows) == 0L ||
      is.null(support_rows$cigar) || is.null(support_rows$read_name)) {
    return(empty)
  }
  has_sc <- grepl("^\\d+S", support_rows$cigar, perl = TRUE) |
    grepl("\\d+S$", support_rows$cigar, perl = TRUE)
  sc_qnames <- unique(support_rows$read_name[has_sc & !is.na(support_rows$read_name)])
  if (length(sc_qnames) == 0L) return(empty)

  pi <- .pair_info(all_pairs, min_mapq)
  if (is.null(pi) || is.null(pi$qname) || length(pi$qname) != pi$n) return(empty)

  outer_insert <- pi$outer_insert
  size_ok <- outer_insert > 0L & outer_insert <= max_insert
  support_pair_mask <- pi$valid & size_ok & pi$qname %in% sc_qnames
  if (!any(support_pair_mask, na.rm = TRUE)) return(empty)

  span_bp <- pi$pair_start <= breakpoint & pi$pair_end >= breakpoint
  bg_mask <- pi$valid & size_ok &
    pi$pair_start >= (breakpoint - flank) & pi$pair_end <= (breakpoint + flank) & !span_bp
  bg_med <- stats::median(outer_insert[bg_mask], na.rm = TRUE)
  bg_mad <- stats::mad(outer_insert[bg_mask], na.rm = TRUE)
  if (!is.finite(bg_mad) || bg_mad == 0) bg_mad <- stats::sd(outer_insert[bg_mask], na.rm = TRUE)
  if (!is.finite(bg_mad) || bg_mad == 0) bg_mad <- 50

  len_used <- if (is.na(best_len) || best_len <= 0L) 0 else as.numeric(best_len)
  event_size_span <- long_span <- rep(FALSE, pi$n)
  if (is.finite(bg_med) && len_used > 0) {
    insert_tol <- max(50, 2 * bg_mad)
    event_size_span <- span_bp & abs(outer_insert - (bg_med + len_used)) <= insert_tol
    long_span <- span_bp & outer_insert >= (bg_med + 0.75 * len_used)
  }

  orient_tab <- table(factor(pi$orientation[support_pair_mask], levels = orient_levels))
  dominant <- if (sum(orient_tab) > 0L) names(sort(orient_tab, decreasing = TRUE))[1L] else NA_character_

  list(
    pe_softclip_support = length(unique(pi$qname[support_pair_mask])),
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

#' Plain data.frame view of an exon GRanges (DRY helper)
#'
#' One place resolves the exon numbering column (\code{exon_idx}, else
#' \code{exon_number}, else order) so the graph, plots and report agree.
#' @param exons_gr GRanges of exons, or a data.frame with start/end.
#' @return data.frame(start, end, exon_idx) in transcript order, or NULL.
#' @keywords internal
.exons_df <- function(exons_gr) {
  if (is.null(exons_gr) || length(exons_gr) == 0L) return(NULL)
  if (is.data.frame(exons_gr)) {
    df <- exons_gr
    if (is.null(df$exon_idx)) df$exon_idx <- seq_len(nrow(df))
  } else {
    mc <- S4Vectors::mcols(exons_gr)
    idx <- if (!is.null(mc$exon_idx)) mc$exon_idx else if (!is.null(mc$exon_number)) mc$exon_number else seq_along(exons_gr)
    df <- data.frame(start = GenomicRanges::start(exons_gr), end = GenomicRanges::end(exons_gr),
                     exon_idx = as.integer(unlist(idx)))
  }
  df <- df[order(df$exon_idx), c("start", "end", "exon_idx"), drop = FALSE]
  rownames(df) <- NULL
  df
}

#' Resolve run settings: explicit argument > YAML gene settings > default
#'
#' Shared by talos() and talos_alu() (DRY). Defaults come straight from the
#' engine's formals, so the wrapper and the engine can no longer drift apart.
#' @param env Calling environment holding the user's arguments (NULL = unset).
#' @param defaults Named list / pairlist of default expressions.
#' @param yaml_vals Named list from the gene YAML (may be NULL).
#' @return Named list of resolved values.
#' @keywords internal
.resolve_settings <- function(env, defaults, yaml_vals = NULL) {
  yaml_vals <- yaml_vals %||% list()
  out <- lapply(names(defaults), function(nm) {
    user_val <- if (exists(nm, envir = env, inherits = FALSE)) base::get(nm, envir = env, inherits = FALSE) else NULL
    if (!is.null(user_val)) return(user_val)
    if (!is.null(yaml_vals[[nm]])) return(yaml_vals[[nm]])
    d <- defaults[[nm]]
    if (is.null(d)) NULL else eval(d, envir = baseenv())
  })
  names(out) <- names(defaults)
  out
}

.has_fastmatch  <- requireNamespace("fastmatch",  quietly = TRUE)
.has_cigarillo  <- requireNamespace("cigarillo",  quietly = TRUE)
.has_pwalign    <- requireNamespace("pwalign",    quietly = TRUE)
.has_biostrings <- requireNamespace("Biostrings", quietly = TRUE)


.prepare_kmers <- function(ref_dna, k) {
  ref_len <- nchar(ref_dna)
  if (ref_len < k) stop("Reference shorter than k-mer length.")

  if (.has_biostrings) {
    dna_obj <- Biostrings::DNAString(ref_dna)
    v       <- Biostrings::Views(dna_obj, start = seq_len(ref_len - k + 1L), width = k)
    kmers   <- as.character(v)
  } else {
    kmer_starts <- seq_len(ref_len - k + 1L)
    kmers       <- substring(ref_dna, kmer_starts, kmer_starts + k - 1L)
  }
  if (.has_fastmatch) fastmatch::fmatch("", kmers)  # prime hash table
  kmers
}


.kmer_match <- function(query, table) {
  if (.has_fastmatch) fastmatch::fmatch(query, table) else match(query, table)
}
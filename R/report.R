# ============================================================================
# TALOS reporting: shared variant/plot data, PDF (Gviz), plotly widget and the
# self-contained HTML report (no R Markdown; SPARKS-style single page).
# ============================================================================

# ── Shared plot geometry ─────────────────────────────────────────────────────

#' Exons and genomic window to plot
#' @param gene_config Gene configuration.
#' @param flank_exons Exons added on each side of the targeted exons.
#' @param flank_bp Bases added on each side of the plotted exons.
#' @param span Optional c(start, end) that must also be shown (e.g. a PTD
#'   reaching beyond the targeted exons).
#' @return list(exons_plot = data.frame(start, end, exon_idx) or NULL,
#'   plot_start, plot_end)
#' @keywords internal
.build_plot_exon_context <- function(gene_config, flank_exons = 1L, flank_bp = 50L, span = NULL) {
  all_ex <- .exons_df(gene_config$all_exons)
  if (is.null(all_ex)) {
    lo <- min(c(gene_config$genomic_start, span[1L]), na.rm = TRUE)
    hi <- max(c(gene_config$genomic_end, span[2L]), na.rm = TRUE)
    return(list(exons_plot = NULL, plot_start = lo, plot_end = hi))
  }
  tgt <- .exons_df(gene_config$target_exons)
  keep <- rep(TRUE, nrow(all_ex))
  if (!is.null(tgt) && nrow(tgt) > 0L) {
    lo_i <- max(min(all_ex$exon_idx), min(tgt$exon_idx) - flank_exons)
    hi_i <- min(max(all_ex$exon_idx), max(tgt$exon_idx) + flank_exons)
    keep <- all_ex$exon_idx >= lo_i & all_ex$exon_idx <= hi_i
  }
  if (!is.null(span) && all(!is.na(span))) {
    keep <- keep | (all_ex$end >= span[1L] & all_ex$start <= span[2L])
  }
  ex <- all_ex[keep, , drop = FALSE]
  ex <- ex[order(ex$start), , drop = FALSE]
  list(
    exons_plot = ex,
    plot_start = max(1L, min(c(ex$start, span[1L]), na.rm = TRUE) - flank_bp),
    plot_end = max(c(ex$end, span[2L]), na.rm = TRUE) + flank_bp
  )
}

.parse_itd_coverage_rle <- function(rle_str) {
  if (is.null(rle_str) || length(rle_str) == 0L || is.na(rle_str) || !nzchar(rle_str)) {
    return(numeric(0))
  }
  sp <- strsplit(strsplit(rle_str, ",", fixed = TRUE)[[1L]], ":", fixed = TRUE)
  sp <- sp[lengths(sp) == 2L]
  v <- suppressWarnings(as.numeric(vapply(sp, `[`, character(1L), 1L)))
  n <- suppressWarnings(as.integer(vapply(sp, `[`, character(1L), 2L)))
  ok <- !is.na(v) & !is.na(n) & n > 0L
  rep(v[ok], n[ok])
}

#' Compressed genomic axis: exons at full scale, introns/flanks shrunk
#' @param exons_plot data.frame(start, end, exon_idx) (or GRanges), or NULL.
#' @keywords internal
.build_compressed_axis <- function(exons_plot, raw_plot_start, raw_plot_end,
                                   intron_display = 40L,
                                   flank_display = 25L,
                                   outer_margin = 20L) {
  ex <- .exons_df(exons_plot)
  if (is.null(ex) || nrow(ex) == 0L) {
    segments <- data.frame(
      raw_start = raw_plot_start, raw_end = raw_plot_end,
      disp_start = 1 + outer_margin,
      disp_end = max(1L, raw_plot_end - raw_plot_start + 1L) + outer_margin,
      type = "linear", exon_label = NA_character_, stringsAsFactors = FALSE
    )
  } else {
    ex <- ex[order(ex$start), , drop = FALSE]
    rs <- re <- dl <- numeric(0); ty <- lab <- character(0)
    add_seg <- function(a, b, len, type, label = NA_character_) {
      if (b < a || len <= 0L) return(invisible())
      rs <<- c(rs, a); re <<- c(re, b); dl <<- c(dl, len); ty <<- c(ty, type); lab <<- c(lab, label)
    }
    if (raw_plot_start < ex$start[1L]) {
      add_seg(raw_plot_start, ex$start[1L] - 1L, min(flank_display, ex$start[1L] - raw_plot_start), "flank")
    }
    for (i in seq_len(nrow(ex))) {
      add_seg(ex$start[i], ex$end[i], ex$end[i] - ex$start[i] + 1L, "exon", paste0("E", ex$exon_idx[i]))
      if (i < nrow(ex)) {
        gs <- ex$end[i] + 1L
        ge <- ex$start[i + 1L] - 1L
        if (ge >= gs) add_seg(gs, ge, min(intron_display, ge - gs + 1L), "intron")
      }
    }
    last_end <- max(ex$end)
    if (raw_plot_end > last_end) {
      add_seg(last_end + 1L, raw_plot_end, min(flank_display, raw_plot_end - last_end), "flank")
    }
    disp_start <- 1 + outer_margin + c(0, cumsum(dl)[-length(dl)])
    segments <- data.frame(raw_start = rs, raw_end = re, disp_start = disp_start,
                           disp_end = disp_start + dl - 1L, type = ty, exon_label = lab,
                           stringsAsFactors = FALSE)
  }
  # Vectorised: findInterval instead of a per-position which().
  transform_pos <- function(pos) {
    pos <- as.numeric(pos)
    i <- findInterval(pos, segments$raw_start)
    out <- rep(NA_real_, length(pos))
    before <- i == 0L
    out[before] <- segments$disp_start[1L]
    j <- pmax(i, 1L)
    inside <- !before & pos <= segments$raw_end[j]
    w <- segments$raw_end[j] - segments$raw_start[j]
    frac <- ifelse(w > 0, (pos - segments$raw_start[j]) / w, 0)
    out[inside] <- (segments$disp_start[j] + frac * (segments$disp_end[j] - segments$disp_start[j]))[inside]
    # Gaps between segments (should not happen) and positions past the end.
    gap <- !before & !inside
    out[gap] <- ifelse(i[gap] >= nrow(segments), segments$disp_end[nrow(segments)], segments$disp_end[j[gap]])
    out
  }
  exon_rows <- which(segments$type == "exon")
  list(
    segments = segments,
    transform_pos = transform_pos,
    disp_from = 0,
    disp_to = max(segments$disp_end) + outer_margin,
    exon_centers = (segments$disp_start[exon_rows] + segments$disp_end[exon_rows]) / 2,
    exon_labels = segments$exon_label[exon_rows]
  )
}

.sample_compressed_coverage <- function(cov, axis_ctx, max_gap_points = 30L,
                                        coverage_start = 1L) {
  seg <- axis_ctx$segments
  out <- lapply(seq_len(nrow(seg)), function(i) {
    raw_n <- as.integer(seg$raw_end[i] - seg$raw_start[i] + 1L)
    raw_pos <- if (seg$type[i] == "exon" || raw_n <= max_gap_points) {
      seq.int(as.integer(seg$raw_start[i]), as.integer(seg$raw_end[i]))
    } else {
      as.integer(unique(round(seq(seg$raw_start[i], seg$raw_end[i], length.out = max_gap_points))))
    }
    depth <- if (is.null(cov)) rep(NA_real_, length(raw_pos)) else {
      # Coverage queried for a report window is stored window-relative to
      # avoid constructing a chromosome-length Rle. Legacy caller-supplied
      # precomputed Rles retain genomic-coordinate indexing via start=1.
      idx <- as.integer(raw_pos - as.numeric(coverage_start) + 1L)
      inr <- !is.na(idx) & idx >= 1L & idx <= length(cov)
      d <- rep(0, length(raw_pos))
      d[inr] <- as.numeric(cov[idx[inr]])
      d
    }
    data.frame(raw_pos = raw_pos, disp_pos = axis_ctx$transform_pos(raw_pos), depth = depth,
               type = seg$type[i], stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Resolve the per-base ITD coverage vector used by report tracks.
#'
#' Prefer the already folded allele-specific display RLE. Older / merged result
#' rows can have that field unset while retaining source-specific RLE columns;
#' reconstruct from direct read-sequence and k-mer evidence only. PE coordinate
#' projection is diagnostic and must never fill an ITD allele coverage gap.
.itd_coverage_vector_for_row <- function(row) {
  read_col <- function(nm) {
    if (!nm %in% names(row)) return(integer(0))
    .parse_itd_coverage_rle(row[[nm]][1L])
  }
  combined <- read_col("ITDCoverageRLE")
  if (length(combined)) return(as.integer(combined))

  direct <- read_col("ITDCoverageRLE_Direct")
  kmer <- read_col("ITDCoverageRLE_Kmer")
  lens <- c(length(direct), length(kmer))
  n <- max(lens)
  if (n == 0L) return(integer(0))
  # Do not align vectors of different lengths by stretching them: only
  # same-coordinate vectors can be safely combined. Keep the best available
  # complete allele-specific vector when an older result is inconsistent.
  direct_ok <- length(direct) == n
  kmer_ok <- length(kmer) == n
  if (length(unique(lens[lens > 0L])) > 1L) {
    if (direct_ok) return(as.integer(direct))
    return(as.integer(kmer))
  }
  .fold_display_coverage(
    if (direct_ok) direct else integer(0),
    integer(0),
    if (kmer_ok) kmer else integer(0)
  )
}

#' Build per-variant ITD-coverage display rows from stored per-base coverage.
#' Missing combined RLEs fall back to source-specific RLEs when available.
#'
#' @param itd_df Data frame with GenomicPosition and coverage RLE columns.
#' @param dup_ends Integer vector of duplication end positions.
#' @param axis_ctx Compressed axis context from .build_compressed_axis().
#' @return data.frame(variant, disp_pos, count, raw_idx); possibly empty.
.build_itd_cov_display <- function(itd_df, dup_ends, axis_ctx) {
  if (is.null(itd_df) || nrow(itd_df) == 0L || !"GenomicPosition" %in% names(itd_df)) {
    return(data.frame())
  }
  coverage_cols <- c("ITDCoverageRLE", "ITDCoverageRLE_Direct",
                     "ITDCoverageRLE_PE", "ITDCoverageRLE_Kmer")
  if (!any(coverage_cols %in% names(itd_df))) return(data.frame())
  rows <- lapply(seq_len(nrow(itd_df)), function(i) {
    cov_vec <- .itd_coverage_vector_for_row(itd_df[i, , drop = FALSE])
    if (length(cov_vec) == 0L) return(NULL)
    x0 <- axis_ctx$transform_pos(as.numeric(itd_df$GenomicPosition[i]))
    x1 <- axis_ctx$transform_pos(as.numeric(dup_ends[i]))
    data.frame(variant = i,
               disp_pos = if (length(cov_vec) == 1L) x0 else seq(x0, x1, length.out = length(cov_vec)),
               count = as.numeric(cov_vec), raw_idx = seq_along(cov_vec))
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) data.frame() else do.call(rbind, rows)
}

# Convert genomic alignment blocks into a coverage Rle whose index 1
# represents `from`. This is deliberately window-relative: report windows can
# include distal exons (e.g. FLT3 E16) that are several Mb from the detection
# interval, and indexing a short target-only Rle by genomic coordinates silently
# rendered those exons as zero.
.coverage_window_from_ranges <- function(gr, bam_chrom, from, to) {
  from <- suppressWarnings(as.integer(from)[1L])
  to <- suppressWarnings(as.integer(to)[1L])
  if (is.na(from) || is.na(to) || from < 1L || to < from) {
    stop("Invalid genomic coverage window.", call. = FALSE)
  }
  width <- as.integer(to - from + 1L)
  if (is.null(gr) || length(gr) == 0L) {
    return(IRanges::coverage(IRanges::IRanges(), width = width))
  }
  if (methods::is(gr, "GRangesList")) gr <- unlist(gr, use.names = FALSE)
  if (!methods::is(gr, "GRanges")) {
    stop("'gr' must be a GRanges or GRangesList object.", call. = FALSE)
  }
  keep_seq <- as.character(GenomicRanges::seqnames(gr)) == as.character(bam_chrom)
  gr <- gr[keep_seq]
  if (length(gr) == 0L) {
    return(IRanges::coverage(IRanges::IRanges(), width = width))
  }
  local_start <- as.integer(GenomicRanges::start(gr)) - from + 1L
  local_end <- as.integer(GenomicRanges::end(gr)) - from + 1L
  keep <- !is.na(local_start) & !is.na(local_end) &
    local_end >= 1L & local_start <= width & local_end >= local_start
  if (!any(keep)) {
    return(IRanges::coverage(IRanges::IRanges(), width = width))
  }
  local_ranges <- IRanges::IRanges(
    start = pmax(1L, local_start[keep]),
    end = pmin(width, local_end[keep])
  )
  IRanges::coverage(local_ranges, width = width)
}

# Coverage over [from, to]: precomputed genomic-coordinate-indexed Rle if
# supplied, otherwise query the BAM and construct an explicitly window-relative
# Rle. The latter avoids relying on the output width/coordinate origin chosen
# by coverage(GAlignments) for a region-restricted read set.
.plot_coverage <- function(bam_path, gene_config, from, to, precomputed_cov = NULL) {
  if (!is.null(precomputed_cov)) {
    finite_precomputed <- tryCatch({
      if (methods::is(precomputed_cov, "Rle")) {
        any(is.finite(S4Vectors::runValue(precomputed_cov)))
      } else {
        any(is.finite(as.numeric(precomputed_cov)))
      }
    }, error = function(e) FALSE)
    # Public precomputed coverage vectors are genomic-coordinate indexed
    # (index 1 = genomic position 1). A shorter vector cannot cover this window.
    covers_window <- length(precomputed_cov) >= as.numeric(to)
    if (isTRUE(finite_precomputed) && isTRUE(covers_window)) {
      return(list(cov = precomputed_cov, start = 1L))
    }
    message("[HTML report] Precomputed genomic coverage is empty, non-finite, or does not span the plotted window; falling back to BAM coverage.")
  }
  if (is.null(bam_path) || length(bam_path) != 1L || is.na(bam_path) ||
      !nzchar(bam_path) || !file.exists(bam_path)) {
    message("[HTML report] Coverage unavailable: no readable BAM path and no usable precomputed coverage.")
    return(list(cov = NULL, start = from))
  }

  tryCatch({
    bam_chrom <- .resolve_bam_seqname(bam_path, gene_config$chrom)
    param <- Rsamtools::ScanBamParam(
      which = GenomicRanges::GRanges(bam_chrom, IRanges::IRanges(from, to)),
      flag = Rsamtools::scanBamFlag(
        isSecondaryAlignment = FALSE,
        isSupplementaryAlignment = FALSE,
        isDuplicate = FALSE,
        isNotPassingQualityControls = FALSE
      ),
      mapqFilter = as.integer(gene_config$min_mapq %||% 20L)
    )
    gal <- GenomicAlignments::readGAlignments(bam_path, param = param)
    gr <- GenomicAlignments::grglist(gal, drop.D.ranges = FALSE)
    list(cov = .coverage_window_from_ranges(gr, bam_chrom, from, to), start = from)
  }, error = function(e) {
    warning(
      "Could not calculate coverage across plotted interval ", from, "-", to,
      " in ", basename(bam_path), ": ", conditionMessage(e),
      call. = FALSE
    )
    list(cov = NULL, start = from)
  })
}

.genome_build_label <- function(gene_config) {
  b <- gene_config$build
  if (is.null(b) || is.na(b)) "custom" else as.character(b)
}

#' Everything a genome-view plot needs for one sample / gene (shared by the
#' PDF, the plotly widget and the HTML report)
#' @param itd_df Results for one sample and gene.
#' @param include_ptd Keep zero-length (clip-only PTD) calls.
#' @return list(itd_df, dup_ends, axis, cov_df, exons, itd_cov) or NULL.
#' @keywords internal
.variant_plot_data <- function(itd_df, gene_config, bam_path = NULL, precomputed_cov = NULL,
                               include_ptd = FALSE, pos_col = "GenomicPosition") {
  if (is.null(itd_df) || nrow(itd_df) == 0L) return(NULL)
  pos <- as.numeric(itd_df[[pos_col]])
  len <- if ("Length" %in% names(itd_df)) as.numeric(itd_df$Length) else rep(0, nrow(itd_df))
  len[is.na(len)] <- 0
  keep <- if (include_ptd) rep(TRUE, nrow(itd_df)) else len > 0
  if (!any(keep)) return(NULL)
  itd_df <- itd_df[keep, , drop = FALSE]
  pos <- pos[keep]; len <- len[keep]
  dup_ends <- pos + pmax(len, 1) - 1
  ctx <- .build_plot_exon_context(gene_config, flank_exons = 1L, flank_bp = 50L,
                                  span = c(min(pos), max(dup_ends)))
  cov_res <- .plot_coverage(bam_path, gene_config, ctx$plot_start, ctx$plot_end, precomputed_cov)
  cov <- cov_res$cov
  cov_start <- cov_res$start
  if (length(cov_start) != 1L || !is.finite(cov_start)) cov_start <- ctx$plot_start
  axis <- .build_compressed_axis(ctx$exons_plot, ctx$plot_start, ctx$plot_end,
                                 intron_display = 40L, flank_display = 30L, outer_margin = 20L)
  list(
    itd_df = itd_df, pos = pos, len = len, dup_ends = dup_ends, axis = axis,
    cov_df = .sample_compressed_coverage(cov, axis, max_gap_points = 30L,
                                        coverage_start = cov_start),
    exons = ctx$exons_plot,
    itd_cov = .build_itd_cov_display(itd_df, dup_ends, axis)
  )
}

# ── Shared variant summary (PDF table, HTML details) ─────────────────────────

.fmt_na <- function(x, fmt = NULL, suffix = "") {
  if (is.null(x) || length(x) == 0L || is.na(x[[1L]]) || identical(as.character(x[[1L]]), "")) return("N/A")
  x <- x[[1L]]
  paste0(if (is.null(fmt)) as.character(x) else sprintf(fmt, x), suffix)
}

#' Human-readable fields of one result row (single source for every report)
#' @param row One-row data.frame (ITD/PTD or ALU result).
#' @param genome_build Build label.
#' @param mode "itd" or "alu".
#' @return Named character vector (label -> value).
#' @keywords internal
.variant_fields <- function(row, genome_build = "custom", mode = c("itd", "alu")) {
  mode <- match.arg(mode)
  g <- function(col) if (col %in% names(row)) row[[col]] else NULL
  vaf <- suppressWarnings(as.numeric(g("AlleleFrequency")))
  hotspot <- if (isTRUE(g("Hotspot"))) .fmt_na(g("HotspotName")) else "No"
  if (mode == "alu") {
    return(c(
      "Sample" = .fmt_na(g("Sample")), "Gene" = .fmt_na(g("Gene")), "Genome Build" = genome_build,
      "Insertion site" = .fmt_na(g("InsertionSite")), "ALU subtype" = .fmt_na(g("ALUSubtype")),
      "Orientation" = .fmt_na(g("ALUOrientation")), "TSD length" = .fmt_na(g("TSD_Length"), suffix = " bp"),
      "TSD sequence" = .fmt_na(g("TSD_Sequence")), "Poly-A tail" = .fmt_na(g("PolyALength"), suffix = " nt"),
      "5' truncation" = .fmt_na(g("ALU5pTruncation"), suffix = " bp"),
      "Supporting reads" = .fmt_na(g("SupportingReads")), "Wildtype reads" = .fmt_na(g("WildtypeReads")),
      "VAF (%)" = .fmt_na(vaf * 100, "%.1f"), "Mean MAPQ" = .fmt_na(g("MeanSupportMAPQ"), "%.1f"),
      "Strand bias" = .fmt_na(g("StrandBias"), "%.3f"), "Region" = .fmt_na(g("Region")),
      "Exon number" = .fmt_na(g("ExonNumber")), "Hotspot" = hotspot
    ))
  }
  pe_len <- if (!is.null(g("LengthPE")) && !is.na(g("LengthPE"))) {
    sprintf("%d bp (%d pairs)", as.integer(g("LengthPE")), as.integer(g("LengthPE_NSpanning") %||% 0L))
  } else "N/A"
  pe_support_len <- if (!is.null(g("LengthPE_SupportQNAME")) &&
                        !is.na(g("LengthPE_SupportQNAME"))) {
    sprintf("%d bp (%d pairs)", as.integer(g("LengthPE_SupportQNAME")),
            as.integer(g("LengthPE_SupportQNAME_NSpanning") %||% 0L))
  } else "N/A"
  n_merged <- suppressWarnings(as.integer(g("MergedFromN") %||% 1L))
  if (is.na(n_merged)) n_merged <- 1L
  merged <- if (n_merged <= 1L) "1 (not merged)" else
    sprintf("%d%s", n_merged, if (isTRUE(g("MergeFragmentWarning"))) " (!! unusually high !!)" else "")
  pe_orient <- if (is.null(g("PEOrientationDominant")) || is.na(g("PEOrientationDominant"))) "N/A" else
    sprintf("%s (FR=%s RF=%s FF=%s RR=%s)", g("PEOrientationDominant"), .fmt_na(g("PEOrientationFR")),
            .fmt_na(g("PEOrientationRF")), .fmt_na(g("PEOrientationFF")), .fmt_na(g("PEOrientationRR")))
  c(
    "Sample" = .fmt_na(g("Sample")), "Gene" = .fmt_na(g("Gene")), "Genome Build" = genome_build,
    "Genomic Position" = .fmt_na(g("GenomicPosition")), "Length (bp)" = .fmt_na(g("Length")),
    "Duplicated exons" = .fmt_na(g("DupExons")),
    "Length PE (all pairs) (bp)" = pe_len,
    "Length PE (support QNAME pairs) (bp)" = pe_support_len,
    "PE all-pair spanning pairs" = .fmt_na(g("LengthPE_AllPairs_NSpanning")),
    "PE support-QNAME spanning pairs" = .fmt_na(g("LengthPE_SupportQNAME_NSpanning")),
    "PE all-pair median spanning insert (bp)" = .fmt_na(g("LengthPE_AllPairs_MedianSpanningInsert")),
    "PE background median insert (bp)" = .fmt_na(g("LengthPE_AllPairs_MedianBackgroundInsert")),
    "PE support-QNAME matched pairs" = .fmt_na(g("LengthPE_SupportQNAME_NMatched")),
    "PE support-QNAME median spanning insert (bp)" = .fmt_na(g("LengthPE_SupportQNAME_MedianSpanningInsert")),
    "PE support-insert estimate (bp; relaxed span)" = .fmt_na(g("LengthPE_SupportInsert")),
    "PE support-insert candidate pairs" = .fmt_na(g("LengthPE_SupportInsert_NPairs")),
    "PE support-insert median insert (bp)" = .fmt_na(g("LengthPE_SupportInsert_MedianInsert")),
    "PE support-insert status" = .fmt_na(g("LengthPE_SupportInsert_Status")),
    "PE support-QNAME estimate status" = .fmt_na(g("LengthPE_SupportQNAME_Status")),
    "PE all-pair estimate status" = .fmt_na(g("LengthPE_AllPairs_Status")),
    "PE all-pair eligibility" = .fmt_na(g("LengthPE_AllPairs_Eligibility")),
    "PE support-QNAME eligibility" = .fmt_na(g("LengthPE_SupportQNAME_Eligibility")),
    "PE support-insert eligibility" = .fmt_na(g("LengthPE_SupportInsert_Eligibility")),
    "PE support-QNAME unavailability reason" = .fmt_na(g("LengthPE_SupportQNAME_Status")),
    "Length Ext (bp)" = .fmt_na(g("LengthExt")),
    "PE Length Delta (PE - call)" = .fmt_na(g("LengthPE_Delta"), "%+.0f", " bp"),
    "PE Length Concordant" = .fmt_na(g("LengthPE_Concordant")),
    "VAF (%) - raw" = .fmt_na(vaf * 100, "%.1f"),
    "VAF (%) - size-corrected diagnostic" = .fmt_na(suppressWarnings(as.numeric(g("CorrectedAlleleFrequency"))) * 100, "%.1f"),
    "In frame (diagnostic)" = .fmt_na(g("InFrame")),
    "Supporting Reads - raw" = .fmt_na(g("SupportingReads")),
    "Supporting Reads - size-corrected diagnostic" = .fmt_na(g("AdjustedSupportingReads")),
    "Size Correction Factor" = .fmt_na(g("SizeCorrectionFactor"), "%.3f"),
    "Independent Supporting QNAMEs" = .fmt_na(g("IndependentSupportQNAMEs")),
    "Anchor-confirmed Supporting QNAMEs" = .fmt_na(g("AnchorConfirmedQNAMEs")),
    "Anchor-confirmation status" = .fmt_na(g("AnchorConfirmedQNAMEsStatus")),
    "Wildtype Reads" = .fmt_na(if (is.null(g("WildtypeReads"))) NULL else round(g("WildtypeReads"))),
    "Depth at Breakpoint - raw" = .fmt_na(g("DepthAtBreakpoint")),
    "Depth at Breakpoint - corrected diagnostic" = .fmt_na(g("CorrectedDepthAtBreakpoint")),
    "Local Coverage" = .fmt_na(g("LocalCoverage")),
    "Span Min Coverage" = .fmt_na(g("SpanMinCoverage")),
    "Span Mean Coverage" = .fmt_na(g("SpanMeanCoverage")),
    "Span Depth Fold Change" = .fmt_na(g("SpanDepthFoldChange"), "%.2f", "x"),
    "HGVS cDNA" = .fmt_na(g("HGVS_cDNA")), "HGVS Protein" = .fmt_na(g("HGVS_Protein")),
    "Region" = .fmt_na(g("Region")), "Exon Number" = .fmt_na(g("ExonNumber")),
    "Orientation" = .fmt_na(g("Orientation")),
    "Observed Fragment Orientation" = .fmt_na(g("ITDObservedFragmentOrientation")),
    "Observed Fragment Identity (%)" = .fmt_na(g("ITDObservedFragmentPID"), "%.1f", "%"),
    "Observed Fragment Query Coverage (%)" = .fmt_na(g("ITDObservedFragmentQueryCoveragePct"), "%.1f", "%"),
    "Observed Fragment Aligned Bases" = .fmt_na(g("ITDObservedFragmentAlignedBases")),
    "Observed Fragment Template Start" = .fmt_na(g("ITDObservedFragmentTemplateStart")),
    "Observed Fragment Template End" = .fmt_na(g("ITDObservedFragmentTemplateEnd")),
    "Observed Fragment Supporting Reads" = .fmt_na(g("ITDObservedFragmentSupportReads")),
    "Observed Fragment Unique Sequences" = .fmt_na(g("ITDObservedFragmentUniqueSequences")),
    "Observed Fragment Status" = .fmt_na(g("ITDObservedFragmentStatus")),
    "RefMatch Observed (%)" = .fmt_na(g("RefMatch_Observed"), "%.1f", "%"),
    "RefMatch Total (%)" = .fmt_na(g("RefMatch_Total"), "%.1f", "%"),
    "Alignment Score" = .fmt_na(g("AlignmentScore"), "%.3f"),
    "Sequence Source" = .fmt_na(g("SequenceSource")),
    "Sequence-dependent Metrics Status" = .fmt_na(g("SequenceDependentMetricsStatus")),
    "Coverage Template Source" = .fmt_na(g("CoverageTemplateSource")),
    "Coverage Template Length" = .fmt_na(g("CoverageTemplateLength")),
    "Adapter Screening Status" = .fmt_na(g("ITDAdapterScreeningStatus")),
    "Adapter-Positive Support Reads" = .fmt_na(g("ITDAdapterPositiveSupportReads")),
    "Adapter-Positive Reads (region)" = .fmt_na(g("ITDAdapterPositiveReads")),
    "Base Quality Status" = .fmt_na(g("ITDBaseQualityStatus")),
    "Mean Base Quality (support reads)" = .fmt_na(g("ITDMeanBaseQuality"), "%.1f"),
    "Min Base Quality (support reads)" = .fmt_na(g("ITDMinBaseQuality"), "%.1f"),
    "Terminal-Window Base Quality" = .fmt_na(g("ITDTerminalMeanBaseQuality"), "%.1f"),
    "ITD Direct Coverage (%)" = .fmt_na(g("ITDReadCoverage_Direct"), "%.1f", "%"),
    "ITD PE Coordinate Projection (%) - diagnostic only" = .fmt_na(g("ITDReadCoverage_PE"), "%.1f", "%"),
    "ITD K-mer Coverage (%)" = .fmt_na(g("ITDReadCoverage_Kmer"), "%.1f", "%"),
    "ITD Combined Coverage (%) - direct + k-mer" = .fmt_na(g("ITDReadCoverage_Combined"), "%.1f", "%"),
    "ITD Combined Depth Method" = .fmt_na(g("ITDCoverageCombinedDepthMethod")),
    "Direct Certain Coverage (%)" = .fmt_na(g("ITDReadCoverage_DirectCertain"), "%.1f", "%"),
    "Direct Possible Coverage (%)" = .fmt_na(g("ITDReadCoverage_DirectPossible"), "%.1f", "%"),
    "Combined Certain Coverage (%)" = .fmt_na(g("ITDReadCoverage_CombinedCertain"), "%.1f", "%"),
    "Combined Possible Coverage (%)" = .fmt_na(g("ITDReadCoverage_CombinedPossible"), "%.1f", "%"),
    "Direct Ambiguous Coverage (%)" = .fmt_na(g("ITDReadCoverage_DirectAmbiguous"), "%.1f", "%"),
    "Direct Covered Bases" = .fmt_na(g("ITDCoverageDirectCoveredBases")),
    "PE Projected Bases (diagnostic only)" = .fmt_na(g("ITDCoveragePECoveredBases")),
    "K-mer Covered Bases" = .fmt_na(g("ITDCoverageKmerCoveredBases")),
    "Combined Covered Bases" = .fmt_na(g("ITDCoverageCombinedCoveredBases")),
    "Direct Mean Depth" = .fmt_na(g("ITDCoverageDirectMeanDepth"), "%.2f"),
    "PE Projected Mean Depth (diagnostic only)" = .fmt_na(g("ITDCoveragePEMeanDepth"), "%.2f"),
    "K-mer Mean Depth" = .fmt_na(g("ITDCoverageKmerMeanDepth"), "%.2f"),
    "Combined Mean Depth" = .fmt_na(g("ITDCoverageCombinedMeanDepth"), "%.2f"),
    "Direct Uncovered Bases" = .fmt_na(g("ITDUncoveredBases_Direct")),
    "PE Projection Uncovered Bases (diagnostic only)" = .fmt_na(g("ITDUncoveredBases_PE")),
    "K-mer Uncovered Bases" = .fmt_na(g("ITDUncoveredBases_Kmer")),
    "Combined Uncovered Bases" = .fmt_na(g("ITDUncoveredBases_Combined")),
    "Combined Uncovered Intervals" = .fmt_na(g("ITDUncoveredIntervals_Combined")),
    "Coverage Read Alignments Processed" = .fmt_na(g("ITDCoverageReadAlignmentsProcessed")),
    "Coverage Read Alignments Matched" = .fmt_na(g("ITDCoverageReadsMatched")),
    "Coverage Fragments Processed" = .fmt_na(g("ITDCoverageFragmentsProcessed")),
    "Coverage Fragments with Direct Alignment" = .fmt_na(g("ITDCoverageFragmentsWithDirectAlignment")),
    "Direct Alignment Rate (%)" = .fmt_na(g("ITDCoverageDirectAlignmentRatePct"), "%.1f", "%"),
    "Direct Median PID (%)" = .fmt_na(g("ITDCoverageDirectMedianPID"), "%.1f", "%"),
    "Direct Mean PID (%)" = .fmt_na(g("ITDCoverageDirectMeanPID"), "%.1f", "%"),
    "Direct Median Aligned Bases" = .fmt_na(g("ITDCoverageDirectMedianAlignedBases")),
    "Direct Median Query Coverage (%)" = .fmt_na(g("ITDCoverageDirectMedianQueryCoveragePct"), "%.1f", "%"),
    "Any Direct Alignment Reads" = .fmt_na(g("ITDCoverageDirectReadsAnyAlignment")),
    "Qualifying Direct Alignment Reads" = .fmt_na(g("ITDCoverageDirectReadsQualifying")),
    "PE Spanning Pairs" = .fmt_na(g("ITDPE_SpanningPairs")),
    "PE QNAME Matched Pairs" = .fmt_na(g("ITDPE_QnameMatchedPairs")),
    "PE Coverage Status (detail)" = .fmt_na(g("ITDCoveragePEStatus")),
    "K-mer Matched Read Ends" = .fmt_na(g("ITDKmerReads")),
    "K-mer Tested Read Ends" = .fmt_na(g("ITDKmerTestedReads")),
    "K-mer Matched Read Ends (%)" = .fmt_na(g("ITDKmerMatchedReadPct"), "%.1f", "%"),
    "K-mer Tested Fragments" = .fmt_na(g("ITDKmerTestedFragments")),
    "K-mer Matched Fragments" = .fmt_na(g("ITDKmerMatchedFragments")),
    "K-mer Matched Fragments (%)" = .fmt_na(g("ITDKmerMatchedFragmentPct"), "%.1f", "%"),
    "K-mer 7 Coverage (%)" = .fmt_na(g("ITDReadCoverage_Kmer7"), "%.1f", "%"),
    "K-mer 9 Coverage (%)" = .fmt_na(g("ITDReadCoverage_Kmer9"), "%.1f", "%"),
    "K-mer 11 Coverage (%)" = .fmt_na(g("ITDReadCoverage_Kmer11"), "%.1f", "%"),
    "K-mer 13 Coverage (%)" = .fmt_na(g("ITDReadCoverage_Kmer13"), "%.1f", "%"),
    "K-mer 15 Coverage (%)" = .fmt_na(g("ITDReadCoverage_Kmer15"), "%.1f", "%"),
    "Direct Ambiguous Read Ends" = .fmt_na(g("ITDDirectAmbiguousReadEnds")),
    "Direct Multi-placement Read Ends" = .fmt_na(g("ITDDirectMultiPlacementReadEnds")),
    "Direct Ends at Alternative Search Limit" = .fmt_na(g("ITDDirectAlternativeSearchLimitReadEnds")),
    "Direct Uncovered Intervals" = .fmt_na(g("ITDUncoveredIntervals_Direct")),
    "PE Uncovered Intervals" = .fmt_na(g("ITDUncoveredIntervals_PE")),
    "K-mer Uncovered Intervals" = .fmt_na(g("ITDUncoveredIntervals_Kmer")),
    "Length PE Status" = .fmt_na(g("LengthPEStatus")),
    "RefMatch Observed Status" = .fmt_na(g("RefMatchObservedStatus")),
    "RefMatch Total Status" = .fmt_na(g("RefMatchTotalStatus")),
    "Alignment Score Status" = .fmt_na(g("AlignmentScoreStatus")),
    "Orientation Status" = .fmt_na(g("OrientationStatus")),
    "Support Consistency (%)" = .fmt_na(g("SupportConsistency"), "%.1f", "%"),
    "Support Consistency Status" = .fmt_na(g("SupportConsistencyStatus")),
    "Coverage Status" = .fmt_na(g("ITDCoverageStatus")),
    "Direct Coverage Status" = .fmt_na(g("ITDCoverageDirectStatus")),
    "PE Coverage Status" = .fmt_na(g("ITDCoveragePEStatus")),
    "K-mer Coverage Status" = .fmt_na(g("ITDCoverageKmerStatus")),
    "Strand Bias" = .fmt_na(g("StrandBias"), "%.3f"),
    "Mean MAPQ" = .fmt_na(g("MeanSupportMAPQ"), "%.1f"),
    "Hotspot" = hotspot, "Evidence" = .fmt_na(g("CandidateTypes")),
    "Merged From N Fragment(s)" = merged,
    "Merged Support Semantics" = .fmt_na(g("MergedSupportSemantics")),
    "PE Orientation" = pe_orient,
    "Left Softclip Count" = .fmt_na(g("LeftSoftclipCount")),
    "Right Softclip Count" = .fmt_na(g("RightSoftclipCount")),
    "Left Softclip % (Support)" = .fmt_na(g("LeftSoftclipPctSupport"), "%.1f"),
    "Right Softclip % (Support)" = .fmt_na(g("RightSoftclipPctSupport"), "%.1f")
  )
}

# Compact variant-detail field set for the PDF report (canonical labels of
# .variant_fields()). The full ~150-metric matrix caused multi-page tables
# whose wrapped rows escaped the page frame; the PDF now prints this curated
# one-page digest instead. The complete matrix remains available in the HTML
# report (Results -> All result values) and in the results TSV.
.COMPACT_VARIANT_FIELDS <- c(
  "Sample", "Gene", "Genome Build", "Genomic Position", "Length (bp)",
  "Duplicated exons", "Length Ext (bp)", "Length PE (all pairs) (bp)",
  "PE Length Concordant", "PE Length Delta (PE - call)",
  "VAF (%) - raw", "VAF (%) - size-corrected diagnostic", "In frame (diagnostic)",
  "Supporting Reads - raw", "Independent Supporting QNAMEs",
  "Anchor-confirmed Supporting QNAMEs", "Wildtype Reads", "Depth at Breakpoint - raw",
  "Local Coverage", "Span Mean Coverage", "ITD Direct Coverage (%)",
  "ITD K-mer Coverage (%)", "ITD Combined Coverage (%) - direct + k-mer", "HGVS cDNA", "HGVS Protein",
  "Region", "Exon Number", "Orientation", "Strand Bias", "Mean MAPQ", "Hotspot",
  "Evidence", "Merged From N Fragment(s)"
)
.COMPACT_ALU_FIELDS <- c(
  "Sample", "Gene", "Genome Build", "Insertion site", "ALU subtype", "Orientation",
  "TSD length", "Poly-A tail", "5' truncation", "Supporting reads", "Wildtype reads",
  "VAF (%)", "Mean MAPQ", "Strand bias", "Region", "Exon number", "Hotspot"
)

#' Curated short field list of one result row for the PDF detail page
#' @inheritParams .variant_fields
#' @return Named character vector (label -> value), subset of .variant_fields().
#' @keywords internal
.compact_variant_fields <- function(row, genome_build = "custom",
                                    mode = c("itd", "alu")) {
  mode <- match.arg(mode)
  fields <- .variant_fields(row, genome_build = genome_build, mode = mode)
  keep <- intersect(if (mode == "alu") .COMPACT_ALU_FIELDS else .COMPACT_VARIANT_FIELDS,
                    names(fields))
  fields[keep]
}

# Configuration summary (PDF page, plotly annotation, HTML tab).
.config_fields <- function(gene_config) {
  gs <- gene_config$gene_settings
  tgt <- .exons_df(gene_config$target_exons)
  c(
    "Gene" = .fmt_na(gene_config$gene), "Build" = .genome_build_label(gene_config),
    "Chromosome" = .fmt_na(gene_config$chrom), "Strand" = .fmt_na(gene_config$strand),
    "Transcript" = if (is.null(gene_config$transcript) || length(gene_config$transcript) != 1L ||
                        is.na(gene_config$transcript) || !nzchar(gene_config$transcript)) "none" else gene_config$transcript,
    "CDS offset" = .fmt_na(gene_config$cds_offset),
    "Targeted exons" = if (is.null(tgt)) "N/A" else paste(unique(tgt$exon_idx), collapse = ", "),
    "Genomic window" = sprintf("%s-%s", gene_config$genomic_start, gene_config$genomic_end),
    "min_support" = as.character(gs$min_support %||% "default"),
    "cluster_tolerance" = as.character(gs$cluster_tolerance %||% "default"),
    "min_mapq" = as.character(gs$min_mapq %||% "default")
  )
}

# Variant-details pagination for the PDF report.
#
# History of this function: gridExtra::tableGrob() draws at its natural size,
# so a summary taller than one page is split into page-sized chunks. The first
# implementation drew the full ~150-metric table in one fixed 0.88-height
# viewport and silently clipped every row outside the page; the second
# implementation measured the row pitch from 1-row / 3-row probe grobs. That
# probe under-reports the real drawn pitch by ~50% (cell padding is only
# resolved when the table is drawn, and long eligibility text wraps to extra
# lines the probe never sees), so chunks were again taller than the page and
# the table ran past the top and bottom page edges ("second and third page
# not in frame").
#
# The current design is deterministic instead of trusting probe measurements:
#   * every cell value is pre-wrapped with strwrap() to a fixed character
#     width, so the number of drawn text lines of every row is known exactly;
#   * rows are packed per page from a line-count height model that slightly
#     over-estimates the drawn pitch, so a chunk can only under-fill a page,
#     never overflow it;
#   * the table carries the compact variant-detail field set
#     (.compact_variant_fields()), which fits a single page for any sample.
# The width is still measured directly (widths are reported reliably) and the
# font shrunk when a wide multi-call table would overflow the page width.
# The no-gridExtra fallback paginates the printed text the same way.
.draw_summary_pages <- function(summary_df, title, width, height) {
  page_h_pt <- height * 72
  page_w_pt <- width * 72
  title_h_pt <- 30
  vp_h_frac <- 0.88
  vp_y_frac <- 0.455
  avail_h_pt <- page_h_pt * vp_h_frac - 6
  avail_w_pt <- page_w_pt * 0.95
  n <- nrow(summary_df)
  title_gp <- grid::gpar(fontsize = 13, fontface = "bold", col = "#2c3e50")

  # Pre-wrap cell values and count the drawn lines per row. tableGrob renders
  # one text line per strwrap segment, so the pagination below can compute the
  # height of every row exactly instead of probing the grob (whose reported
  # height under-estimates the drawn pitch; see the header comment). Line
  # counts are taken from the original values BEFORE the wrapped copy is
  # built: strwrap() would collapse the inserted newlines back into one
  # segment and under-count again.
  cell_lines <- function(v) {
    v <- as.character(v)
    if (is.na(v) || !nzchar(v)) return(1L)
    max(1L, length(strwrap(v, width = 64)))
  }
  row_lines <- pmax(1L, apply(summary_df, 1L, function(r) max(vapply(r, cell_lines, integer(1L)))))
  wrap_cell <- function(v) {
    v <- as.character(v)
    if (is.na(v) || !nzchar(v)) return(v)
    parts <- strwrap(v, width = 64)
    if (length(parts) <= 1L) return(v)
    paste(parts, collapse = "\n")
  }
  summary_df <- as.data.frame(
    lapply(summary_df, function(col_chr) vapply(col_chr, wrap_cell, character(1L))),
    stringsAsFactors = FALSE
  )

  page_title <- function(page_i, npages, continued) {
    ttl <- title
    if (npages > 1L) {
      ttl <- paste0(ttl, sprintf("  (%d/%d)", page_i, npages))
      if (continued) ttl <- paste0(ttl, ", continued")
    }
    grid::grid.text(ttl, x = 0.5, y = grid::unit(0.975, "npc"),
                    just = c("centre", "top"), gp = title_gp)
  }

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    theme_for <- function(core_cex, head_cex, nr) gridExtra::ttheme_minimal(
      core = list(fg_params = list(cex = core_cex, fontface = "plain"),
                  bg_params = list(fill = rep(c("white", "#f0f4f8"),
                                               length.out = max(1L, nr)), col = NA)),
      colhead = list(fg_params = list(cex = head_cex, fontface = "bold", col = "white"),
                     bg_params = list(fill = "#2c3e50", col = NA)))
    grob_width_pt <- function(df, th) {
      g <- gridExtra::tableGrob(df, rows = NULL, theme = th)
      w <- tryCatch(as.numeric(grid::convertUnit(grid::grobWidth(g), "pt", valueOnly = TRUE)),
                    error = function(e) NA_real_)
      w
    }

    th <- theme_for(0.7, 0.8, n)
    full_w <- grob_width_pt(summary_df, th)
    if (is.finite(full_w) && full_w > avail_w_pt && full_w > 0) {
      scale <- max(0.45, avail_w_pt / full_w * 0.98)
      th <- theme_for(0.7 * scale, 0.8 * scale, n)
    }

    # Deterministic pagination from the line counts: a cex-0.7 tableGrob row
    # draws ~17.3 pt per text line and the header ~17 pt, so 18 pt/line with a
    # 20 pt header and 6 pt slack over-estimates the drawn height and can only
    # under-fill a page, never overflow it.
    line_pitch_pt <- 18
    head_pt <- 20
    budget_pt <- avail_h_pt - head_pt - 6
    chunks <- list()
    start_i <- 1L
    while (start_i <= n) {
      used <- 0
      end_i <- start_i - 1L
      while (end_i < n) {
        next_h <- row_lines[end_i + 1L] * line_pitch_pt
        if (used + next_h > budget_pt && end_i >= start_i) break
        used <- used + next_h
        end_i <- end_i + 1L
      }
      chunks[[length(chunks) + 1L]] <- start_i:end_i
      start_i <- end_i + 1L
    }
    npages <- length(chunks)
    for (p in seq_along(chunks)) {
      chunk <- summary_df[chunks[[p]], , drop = FALSE]
      tbl <- gridExtra::tableGrob(chunk, rows = NULL, theme = th)
      grid::grid.newpage()
      page_title(p, npages, p > 1L)
      grid::pushViewport(grid::viewport(x = 0.5, y = vp_y_frac, width = 0.96, height = vp_h_frac))
      grid::grid.draw(tbl)
      grid::popViewport()
    }
  } else {
    lns <- utils::capture.output(print(summary_df))
    line_h_pt <- 12 * 0.7 * 1.35
    per_page <- max(8L, floor((page_h_pt * vp_h_frac - title_h_pt) / line_h_pt))
    chunks <- split(lns, ceiling(seq_along(lns) / per_page))
    npages <- length(chunks)
    for (p in seq_along(chunks)) {
      grid::grid.newpage()
      page_title(p, npages, p > 1L)
      grid::grid.text(paste(chunks[[p]], collapse = "\n"), x = 0.04, y = 0.90,
                      just = c("left", "top"),
                      gp = grid::gpar(cex = 0.7, fontfamily = "mono"))
    }
  }
  invisible(NULL)
}

# ── PDF (Gviz) ───────────────────────────────────────────────────────────────

#' Generate a publication-ready PDF plot (Gviz + compact variant details)
#' Only the region restricted directly by targeted_exons context is shown
#' (extended to include the plotted variants). Page 1 carries the Gviz
#' tracks; the following pages carry a compact variant-details table
#' (.compact_variant_fields()). The full ~150-metric matrix is not printed
#' any more: it caused multi-page tables whose rows escaped the page frame;
#' it remains available in the HTML report and the results TSV.
#'
#' @param itd_row Data frame of TALOS results for one sample/gene.
#' @param gene_config Gene configuration list.
#' @param bam_path Path to BAM file.
#' @param sample_name Sample name.
#' @param output_pdf Output PDF file path.
#' @param width Plot dimension (inches).
#' @param height Plot dimension (inches).
#' @param show_config If TRUE, display configuration summary on a separate page.
#' @param precomputed_cov Optional precomputed coverage Rle indexed by genomic
#'   coordinate (index 1 is genomic position 1) and long enough to include the
#'   full plotted interval, including flanking exons.
#' @param include_ptd If TRUE, include zero-length (PTD) events (default FALSE).
#' @param gap Gap between Gviz tracks.
#' @return Invisibly, output_pdf.
#' @export
plot_talos_report <- function(itd_row, gene_config, bam_path,
                              sample_name, output_pdf,
                              width = 14, height = 14,
                              show_config = FALSE,
                              precomputed_cov = NULL,
                              include_ptd = FALSE,
                              gap = 0.05) {
  if (!requireNamespace("Gviz", quietly = TRUE)) {
    warning("Gviz not installed - PDF plot skipped.")
    return(invisible(NULL))
  }
  if (nrow(itd_row) == 0) {
    warning("No ITDs in itd_df.")
    return(invisible(NULL))
  }
  pd <- .variant_plot_data(itd_row, gene_config, bam_path, precomputed_cov, include_ptd)
  if (is.null(pd)) {
    warning("No ITDs with positive length to plot (PTDs excluded). Set include_ptd=TRUE to include them.")
    return(invisible(NULL))
  }
  itd_df <- pd$itd_df
  axis_ctx <- pd$axis
  chrom <- gene_config$chrom
  genome_build <- .genome_build_label(gene_config)
  breakpoints <- pd$pos
  dup_lens <- pd$len
  vafs <- sprintf("%.1f", itd_df$AlleleFrequency * 100)
  orient <- if ("Orientation" %in% names(itd_df)) itd_df$Orientation else rep(NA_character_, nrow(itd_df))
  orient <- ifelse(is.na(orient) | orient == "?", "", orient)

  cov_df <- pd$cov_df
  y_max <- max(cov_df$depth, na.rm = TRUE)
  if (!is.finite(y_max) || y_max == 0) y_max <- 1L
  title_args <- list(background.title = "transparent", col.border.title = NA, fontcolor.title = "black")

  cov_track <- do.call(Gviz::DataTrack, c(list(
    start = cov_df$disp_pos, end = cov_df$disp_pos, data = cov_df$depth, chromosome = chrom,
    genome = genome_build, type = "histogram", name = "Coverage",
    col.histogram = "steelblue", fill.histogram = "#dfeefa",
    ylab = "Read depth", ylim = c(0, y_max * 1.05),
    grid = TRUE, lwd.grid = 0.4, col.grid = "#d0d0d0",
    cex.axis = 0.8, cex.title = 0.9, labelPos = "below"
  ), title_args))

  grtrack <- NULL
  if (!is.null(pd$exons) && nrow(pd$exons) > 0) {
    exon_gr <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges = IRanges::IRanges(start = round(axis_ctx$transform_pos(pd$exons$start)),
                                end = round(axis_ctx$transform_pos(pd$exons$end))),
      strand = if (gene_config$strand %in% c("+", "-")) gene_config$strand else "*"
    )
    exon_gr$gene <- gene_config$gene
    exon_gr$transcript <- paste0(gene_config$gene, "_tx")
    exon_gr$symbol <- gene_config$gene
    exon_gr$exon <- paste0("E", pd$exons$exon_idx)
    grtrack <- do.call(Gviz::GeneRegionTrack, c(list(
      exon_gr, genome = genome_build, chromosome = chrom,
      name = paste0(gene_config$gene, " exons (+ padding)"),
      fill = "#2ecc71", col = "#1a8a4a", transcriptAnnotation = "none",
      exonAnnotation = "exon", showExonId = TRUE, fontsize = 11, cex.feature = 1.0,
      fontcolor.feature = "black", shape = "box", stackHeight = 0.45
    ), title_args))
  }

  itd_cov_track <- NULL
  if (nrow(pd$itd_cov) > 0L) {
    itd_cov_ok <- is.finite(pd$itd_cov$disp_pos) & is.finite(pd$itd_cov$count)
    itd_cov_plot <- pd$itd_cov[itd_cov_ok, , drop = FALSE]
    if (nrow(itd_cov_plot) > 0L) {
      # Gviz drops / mis-scales a DataTrack with ylim=c(0, 0). Keep a visible
      # zero-depth baseline and coerce display coordinates to finite integers.
      itd_cov_max <- max(1, itd_cov_plot$count, na.rm = TRUE)
      itd_cov_track <- tryCatch(do.call(Gviz::DataTrack, c(list(
        start = pmax(1L, as.integer(round(itd_cov_plot$disp_pos))),
        end = pmax(1L, as.integer(round(itd_cov_plot$disp_pos))),
        data = as.numeric(itd_cov_plot$count),
        chromosome = chrom, genome = genome_build, type = "histogram", name = "ITD cov reads",
        col.histogram = "#f39c12", fill.histogram = "#fdebd0", lwd = 2,
        ylim = c(0, itd_cov_max * 1.05), ylab = "ITD cov"
      ), title_args)), error = function(e) {
        warning("Could not construct the ITD coverage PDF track: ", conditionMessage(e))
        NULL
      })
    }
  }

  itd_labels <- ifelse(nzchar(orient),
                       sprintf("%d bp | VAF %s%% | %s", as.integer(dup_lens), vafs, orient),
                       sprintf("%d bp | VAF %s%%", as.integer(dup_lens), vafs))
  itd_track <- do.call(Gviz::AnnotationTrack, c(list(
    start = round(axis_ctx$transform_pos(breakpoints)), end = round(axis_ctx$transform_pos(pd$dup_ends)),
    chromosome = chrom, genome = genome_build, name = if (length(breakpoints) > 1) "ITDs" else "ITD",
    id = itd_labels, group = factor(seq_along(breakpoints)), stacking = "full", stackHeight = 0.45,
    showFeatureId = TRUE, just = "right", fill = "#f39c12", col = "#c0392b",
    fontcolor.feature = "#4a1c00", cex.feature = max(0.55, 0.8 - 0.04 * length(breakpoints)),
    cex.title = 0.9, rotation.title = 0
  ), title_args))

  axis_track <- Gviz::GenomeAxisTrack(cex = 0.7, labelPos = "below", add35 = FALSE, add53 = FALSE,
                                      col = "black", at = axis_ctx$exon_centers, labels = axis_ctx$exon_labels)
  itd_h <- max(1.0, 0.4 * length(breakpoints))
  tracks <- list(axis_track, grtrack, cov_track, itd_cov_track, itd_track)
  sizes <- c(0.5, 0.65, 2.2, 0.75, itd_h)
  present <- !vapply(tracks, is.null, logical(1L))

  grDevices::pdf(output_pdf, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  main_title <- sprintf("Sample: %s  |  %s  |  %d ITD(s)  |  %s", sample_name, gene_config$gene,
                        length(breakpoints), genome_build)
  subtitle <- sprintf("Positions: %s  |  Lengths: %s  |  VAFs: %s%%",
                      paste(breakpoints, collapse = ", "), paste(dup_lens, collapse = ", "),
                      paste(vafs, collapse = ", "))
  # Keep the subtitle inside Gviz's title margin. Drawing it later at a fixed
  # page coordinate (0.885 npc) placed it over the genome tracks on tall plots.
  main_title <- paste(main_title, subtitle, sep = "\n")
  Gviz::plotTracks(tracks[present], from = axis_ctx$disp_from, to = axis_ctx$disp_to,
                   sizes = sizes[present], collapse = FALSE, main = main_title, cex.main = 1.05,
                   fontcolor.title = "#2c3e50", col.axis = "black", cex.axis = 0.8,
                   margins = c(18, 14, 6.5, 14), gap = gap)

  # Compact variant-details page: the PDF carries a curated digest of every
  # call; the full result matrix stays in the HTML report and the TSV.
  fields <- lapply(seq_len(nrow(itd_df)), function(i)
    .compact_variant_fields(itd_df[i, , drop = FALSE], genome_build))
  summary_df <- data.frame(Metric = names(fields[[1L]]), stringsAsFactors = FALSE)
  for (i in seq_along(fields)) summary_df[[paste0("ITD_", i)]] <- unname(fields[[i]])
  summary_df[summary_df$Metric == "Sample", -1L] <- as.character(sample_name)

  .draw_summary_pages(
    summary_df,
    sprintf("Variant Details (compact) - Sample: %s | %s (%d variants)", sample_name,
            gene_config$gene, nrow(itd_df)),
    width, height
  )

  if (show_config) {
    grid::grid.newpage()
    grid::grid.text("Analysis Configuration", x = 0.5, y = 0.95,
                    gp = grid::gpar(fontsize = 14, fontface = "bold", col = "#2c3e50"))
    cf <- .config_fields(gene_config)
    lines <- paste0(names(cf), ": ", cf)
    for (k in seq_along(lines)) {
      grid::grid.text(lines[k], x = 0.05, y = 0.85 - 0.04 * (k - 1), just = "left",
                      gp = grid::gpar(fontsize = 10, fontfamily = "mono"))
    }
  }
  message("PDF plot saved to: ", output_pdf)
  invisible(output_pdf)
}

# ── plotly htmlwidget (kept for interactive R sessions) ──────────────────────

.TALOS_PALETTE <- c("#e67e22", "#e74c3c", "#9b59b6", "#1abc9c", "#3498db",
                    "#f39c12", "#c0392b", "#8e44ad", "#16a085", "#2980b9")

# Short hover text for one variant (plotly widget and HTML report).
.variant_hover <- function(row, i, vaf_str) {
  f <- .variant_fields(row)
  # Resolve fields by their current canonical labels. Older labels such as
  # "RefMatch (%)" and "ITD Read Coverage (%)" are not emitted by
  # .variant_fields(); indexing those names with [[ ]] raises subscript-out-of-
  # bounds and caused the entire genome panel to be discarded after coverage
  # had already been calculated.
  field <- function(name, fallback = "N/A") {
    if (name %in% names(f) && length(f[[name]]) == 1L && !is.na(f[[name]])) {
      as.character(f[[name]])
    } else fallback
  }
  direct_cov <- field("ITD Direct Coverage (%)")
  combined_cov <- field("ITD Combined Coverage (%)")
  cov_line <- if (combined_cov != "N/A" && direct_cov != "N/A" &&
                  combined_cov != direct_cov) {
    paste0(direct_cov, " (combined ", combined_cov, ")")
  } else if (combined_cov != "N/A") combined_cov else direct_cov

  paste0("<b>ITD ", i, "</b><br>",
         "Breakpoint: ", field("Genomic Position"), "<br>",
         "Length: ", field("Length (bp)"), " bp<br>",
         "VAF: ", vaf_str, "<br>",
         "HGVS cDNA: ", field("HGVS cDNA"), "<br>",
         "HGVS Protein: ", field("HGVS Protein"), "<br>",
         "Region: ", field("Region"), "<br>",
         if (field("Duplicated exons") != "N/A") paste0("Exons: ", field("Duplicated exons"), "<br>") else "",
         "Adjusted support: ", field("Supporting Reads"), " reads<br>",
         "Independent QNAMEs: ", field("Independent Supporting QNAMEs"), "<br>",
         "Anchor-confirmed QNAMEs: ", field("Anchor-confirmed Supporting QNAMEs"), "<br>",
         "RefMatch: ", field("RefMatch Observed (%)"), "<br>",
         "ITD Cov: ", cov_line, "<br>",
         "Orientation: ", sub("N/A", "?", field("Orientation")), "<br>",
         "Left Softclip: ", field("Left Softclip Count"), " (",
         field("Left Softclip % (Support)"), "% of support)<br>",
         "Right Softclip: ", field("Right Softclip Count"), " (",
         field("Right Softclip % (Support)"), "% of support)")
}

#' Generate an interactive plotly visualisation for TALOS results
#' Only the region restricted by targeted exons is shown.
#'
#' @param itd_df Data frame of TALOS results.
#' @param gene_config Gene configuration list.
#' @param bam_path Path to indexed BAM file.
#' @param sample_name Sample identifier.
#' @param show_config Add a config annotation string.
#' @param precomputed_cov Optional precomputed coverage Rle indexed by genomic
#'   coordinate (index 1 is genomic position 1) and long enough to include the
#'   full plotted interval, including flanking exons.
#' @param include_ptd If TRUE, include zero-length (PTD) events (default FALSE).
#' @return A plotly htmlwidget.
#' @export
plot_talos_interactive <- function(itd_df, gene_config, bam_path,
                                   sample_name = NULL,
                                   show_config = FALSE,
                                   precomputed_cov = NULL,
                                   include_ptd = FALSE) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    warning("plotly not installed.")
    return(invisible(NULL))
  }
  if (nrow(itd_df) == 0) return(invisible(NULL))
  pd <- .variant_plot_data(itd_df, gene_config, bam_path, precomputed_cov, include_ptd)
  if (is.null(pd)) {
    warning("No ITDs with positive length to plot interactively (PTDs excluded). Set include_ptd=TRUE to include them.")
    return(plotly::plot_ly())
  }
  itd_df <- pd$itd_df
  ax <- pd$axis
  xr <- c(ax$disp_from, ax$disp_to)
  genome_build <- .genome_build_label(gene_config)
  hide_x <- list(title = "", showticklabels = FALSE, range = xr, rangeslider = list(visible = FALSE))

  p_cov <- plotly::plot_ly(x = pd$cov_df$disp_pos, y = pd$cov_df$depth, type = "scatter", mode = "lines",
                           fill = "tozeroy", name = "Coverage", line = list(color = "#4a90d9", width = 1.2),
                           fillcolor = "rgba(74,144,217,0.22)", hovertemplate = "Depth: %{y}<extra></extra>")
  p_cov <- plotly::layout(p_cov, yaxis = list(title = "Depth", gridcolor = "#e8e8e8"), xaxis = hide_x)

  p_exon <- plotly::plot_ly()
  if (!is.null(pd$exons) && nrow(pd$exons) > 0) {
    for (i in seq_len(nrow(pd$exons))) {
      es <- ax$transform_pos(pd$exons$start[i]); ee <- ax$transform_pos(pd$exons$end[i])
      p_exon <- plotly::add_trace(p_exon, x = c(es, ee, ee, es, es), y = c(0.1, 0.1, 0.9, 0.9, 0.1),
                                  type = "scatter", mode = "lines", fill = "toself",
                                  fillcolor = "rgba(46,204,113,0.7)", line = list(color = "#1a8a4a", width = 1),
                                  name = paste0("Exon ", pd$exons$exon_idx[i]), showlegend = FALSE,
                                  hovertemplate = sprintf("<b>Exon %s</b><br>%d-%d<extra></extra>",
                                                          pd$exons$exon_idx[i], pd$exons$start[i], pd$exons$end[i]))
      p_exon <- plotly::add_annotations(p_exon, x = (es + ee) / 2, y = 0.5, xref = "x", yref = "y",
                                        text = paste0("E", pd$exons$exon_idx[i]), showarrow = FALSE,
                                        font = list(size = 10, color = "#1a5e30"))
    }
  }
  p_exon <- plotly::layout(p_exon, xaxis = hide_x,
                           yaxis = list(title = "Exons", showticklabels = FALSE, range = c(0, 1), fixedrange = TRUE))

  p_itd_cov <- plotly::plot_ly()
  for (vid in unique(pd$itd_cov$variant)) {
    sub <- pd$itd_cov[pd$itd_cov$variant == vid, , drop = FALSE]
    p_itd_cov <- plotly::add_trace(p_itd_cov, x = sub$disp_pos, y = sub$count, type = "scatter", mode = "lines",
                                   fill = "tozeroy", fillcolor = "rgba(243,156,18,0.22)",
                                   line = list(color = .TALOS_PALETTE[((vid - 1L) %% 10L) + 1L], width = 2),
                                   name = paste0("ITD cov ", vid), text = sub$raw_idx, showlegend = FALSE,
                                   hovertemplate = "ITD pos: %{text}<br>Coverage reads: %{y}<extra></extra>")
  }
  ymax <- if (nrow(pd$itd_cov) > 0) max(pd$itd_cov$count, na.rm = TRUE) else 1
  p_itd_cov <- plotly::layout(p_itd_cov, xaxis = hide_x,
                              yaxis = list(title = "ITD cov", range = c(0, max(1, ymax)), fixedrange = TRUE))

  y_step <- max(0.10, min(0.18, 0.85 / max(1L, nrow(itd_df))))
  p_itd <- plotly::plot_ly()
  for (i in seq_len(nrow(itd_df))) {
    row <- itd_df[i, , drop = FALSE]
    vaf_str <- sprintf("%.1f%%", row$AlleleFrequency * 100)
    col_i <- .TALOS_PALETTE[((i - 1L) %% 10L) + 1L]
    y_pos <- 0.05 + (i - 1L) * y_step
    label <- sprintf("%d bp  VAF %s", as.integer(pd$len[i]), vaf_str)
    if (!is.na(row$HGVS_cDNA)) label <- paste0(label, "  ", row$HGVS_cDNA)
    x0 <- ax$transform_pos(pd$pos[i]); x1 <- ax$transform_pos(pd$dup_ends[i])
    p_itd <- plotly::add_trace(p_itd, x = c(x0, x1), y = c(y_pos, y_pos), type = "scatter", mode = "lines",
                               line = list(color = col_i, width = 14), name = label,
                               text = .variant_hover(row, i, vaf_str), hoverinfo = "text", showlegend = TRUE)
    p_itd <- plotly::add_trace(p_itd, x = x0, y = y_pos, type = "scatter", mode = "markers",
                               marker = list(symbol = "line-ns", size = 14, color = col_i,
                                             line = list(color = "white", width = 2)),
                               hoverinfo = "skip", showlegend = FALSE)
  }
  p_itd <- plotly::layout(
    p_itd,
    xaxis = list(title = paste0(gene_config$chrom, " (compressed ", genome_build, ")"),
                 gridcolor = "#e8e8e8", range = xr, rangeslider = list(visible = FALSE)),
    yaxis = list(title = "ITDs", showticklabels = FALSE,
                 range = c(0, max(1.1, nrow(itd_df) * y_step + 0.30)), fixedrange = TRUE),
    legend = list(orientation = "h", x = 0, y = -0.18, font = list(size = 9), tracegroupgap = 0)
  )

  title_text <- paste0("TALOS | ", if (!is.null(sample_name)) paste0(sample_name, "  |  "),
                       gene_config$gene, "  |  ", nrow(itd_df), " ITD(s)  |  ", genome_build)
  fig <- plotly::subplot(p_cov, p_exon, p_itd_cov, p_itd, nrows = 4, shareX = TRUE,
                         heights = c(0.26, 0.10, 0.09, 0.55))
  fig <- plotly::layout(fig, title = list(text = title_text, x = 0.5, font = list(size = 14)),
                        hovermode = "x unified", showlegend = TRUE,
                        margin = list(t = 60, l = 60, r = 40, b = 75),
                        plot_bgcolor = "white", paper_bgcolor = "white")
  if (show_config) {
    cf <- .config_fields(gene_config)
    config_text <- paste("Config:", paste0("Targeted exons: ", cf[["Targeted exons"]]),
                         paste0("Transcript: ", cf[["Transcript"]]), paste0("CDS offset: ", cf[["CDS offset"]]),
                         sep = " | ")
    fig <- plotly::layout(fig, annotations = list(text = config_text, x = 0, y = 1.05, xref = "paper",
                                                  yref = "paper", showarrow = FALSE, font = list(size = 9)))
  }
  fig
}

# ── Self-contained HTML report (no R Markdown) ───────────────────────────────

# Internal columns kept in the downloadable table but hidden by default.
# "Genome" is redundant with GenomeBuild (which the report appends next to
# Gene); the TSD_RawSeq_* columns never existed in the ALU output.
.HTML_HIDDEN_COLS <- c("ITD_Sequence", "ITDCoverageRLE", "SequenceImputed", "SequencePartial",
                       "SupportConsistency", "SpanMeanCoverage", "Genome",
                       "ALU_RawScore", "ALU_Sequence")

# Stylesheet / app script of the HTML report. Embedded in R (report_assets.R);
# options(TALOS.report_dir = "<dir>") swaps in external copies while editing.
.report_asset <- function(name) {
  dev_dir <- getOption("TALOS.report_dir", "")
  dev_file <- file.path(dev_dir, name)
  if (nzchar(dev_dir) && file.exists(dev_file)) {
    return(paste(readLines(dev_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  }
  switch(name,
         talos_report.css = .TALOS_REPORT_CSS,
         talos_report.js = .TALOS_REPORT_JS,
         stop("Unknown report asset: ", name, call. = FALSE))
}


# Browser-side coverage enhancement. Installed before the report app so it can
# decorate each Plotly coverage trace at creation time without changing the
# underlying coverage values or the variant-calling algorithm.
.talos_coverage_enhancement_js <- function() {
  r"---((function () {
  'use strict';
  if (!window.Plotly) return;
  ['newPlot', 'react'].forEach(function (method) {
    var original = window.Plotly[method];
    if (typeof original !== 'function' || original.__talosCoverageWrapped) return;
    var wrapped = function (gd, data, layout, config) {
      if (Array.isArray(data)) {
        var cov = data.find(function (trace) {
          return trace && trace.name === 'Coverage' && Array.isArray(trace.y) &&
            Array.isArray(trace.x) && trace.y.length === trace.x.length;
        });
        var already = data.some(function (trace) {
          return trace && trace.meta && trace.meta.talosCoverageTrend;
        });
        if (cov && !already && cov.y.length >= 7) {
          var ys = cov.y.map(function (v) {
            return (typeof v === 'number' && Number.isFinite(v)) ? v : null;
          });
          var half = Math.max(3, Math.min(10, Math.floor(ys.length / 30)));
          var smoothed = ys.map(function (_, i) {
            var win = ys.slice(Math.max(0, i - half), Math.min(ys.length, i + half + 1))
              .filter(function (v) { return v !== null; }).sort(function (a, b) { return a - b; });
            if (!win.length) return null;
            var mid = Math.floor(win.length / 2);
            return win.length % 2 ? win[mid] : (win[mid - 1] + win[mid]) / 2;
          });
          data = data.concat([{
            type: 'scatter', mode: 'lines', x: cov.x.slice(), y: smoothed,
            name: 'Coverage trend (rolling median)', hoverinfo: 'skip',
            line: { color: '#d97706', width: 2, dash: 'dot' },
            connectgaps: false, showlegend: true, meta: { talosCoverageTrend: true },
            yaxis: cov.yaxis || 'y'
          }]);
        }
      }
      return original.call(this, gd, data, layout, config);
    };
    wrapped.__talosCoverageWrapped = true;
    window.Plotly[method] = wrapped;
  });
})();)---"
}

# Locate a local Plotly JavaScript bundle without requiring the optional R
# package `plotly`. The HTML report can always fall back to the CDN, so a
# missing bundle must never abort report generation.
.local_plotly_js <- function() {
  is_js_file <- function(x) {
    length(x) == 1L && !is.na(x) && nzchar(x) &&
      file.exists(x) && !dir.exists(x) && file.access(x, 4L) == 0L
  }

  configured <- getOption("TALOS.plotly_js", NULL)
  if (is.null(configured) || !length(configured)) {
    configured <- Sys.getenv("TALOS_PLOTLY_JS", unset = "")
  }
  if (length(configured) == 1L && !is.na(configured) && nzchar(configured)) {
    configured <- path.expand(configured)
    if (dir.exists(configured)) {
      configured <- file.path(configured, c("plotly.min.js", "plotly-latest.min.js"))
    }
    configured <- as.character(configured)
    hit <- configured[vapply(configured, is_js_file, logical(1L))]
    if (length(hit)) return(normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE))
  }

  # Plotly has used both names across releases. `system.file()` is safe when
  # the optional package is not installed (it returns an empty string).
  candidates <- unlist(lapply(c("plotly", "htmlwidgets"), function(pkg) {
    root <- tryCatch(system.file(package = pkg), error = function(e) "")
    if (!length(root) || is.na(root) || !nzchar(root)) return(character(0))
    file.path(root, c(
      "htmlwidgets/lib/plotlyjs/plotly.min.js",
      "htmlwidgets/lib/plotlyjs/plotly-latest.min.js",
      "htmlwidgets/lib/plotlyjs/plotly.js"
    ))
  }), use.names = FALSE)
  hit <- candidates[vapply(candidates, is_js_file, logical(1L))]
  if (!length(hit)) return(NA_character_)
  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

# plotly.js <script> tag.
# Inline mode embeds the library as a base64 data: URI. The previous version
# pasted readLines() output between <script> tags: readLines() drops what
# follows an embedded NUL on a line, and plotly's minified bundle has binary
# string literals (PNG signature), so the inlined copy was truncated, threw a
# SyntaxError and window.Plotly stayed undefined ("Interactive plots need
# plotly.js"). Bytes are now read with readBin() and never re-encoded, which
# also makes "</script>" escaping unnecessary. The page still falls back to the
# CDN at runtime if no copy is available.
.plotly_script_tag <- function(mode = c("auto", "inline", "cdn")) {
  mode <- match.arg(mode)
  cdn <- "<script src=\"https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js\"></script>"
  if (mode == "cdn") return(list(tag = cdn, inline = FALSE))
  path <- tryCatch(.local_plotly_js(), error = function(e) NA_character_)
  if (!length(path) || is.na(path[[1L]]) || !nzchar(path[[1L]])) path <- NA_character_
  if (is.na(path)) {
    if (mode == "inline") warning("plotly.js not found locally (install the R 'plotly' package or set options(TALOS.plotly_js=)); using the CDN.")
    return(list(tag = cdn, inline = FALSE))
  }
  inline <- tryCatch({
    size <- file.info(path)$size
    if (length(size) != 1L || is.na(size) || size < 1) stop("empty or unreadable Plotly bundle")
    bytes <- readBin(path, what = "raw", n = size)
    if (!length(bytes)) stop("Plotly bundle could not be read")
    b64 <- jsonlite::base64_enc(bytes)
    b64 <- gsub("[\r\n]", "", b64)
    list(tag = paste0("<script src=\"data:text/javascript;base64,", b64, "\"></script>"),
         inline = TRUE)
  }, error = function(e) {
    message(.html_error_diagnostic(e, context = "HTML report / Plotly bundle"))
    NULL
  })
  if (!is.null(inline)) return(inline)
  warning("Could not inline the local plotly.js bundle; using the CDN.", call. = FALSE)
  list(tag = cdn, inline = FALSE)
}

.arr <- function(x) I(as.vector(x))

# data.frame -> column-oriented JSON-ready list.
# Serialize every source result column for the browser-side Variant Details
# table. Values become stable display strings while preserving column order;
# hidden/internal fields remain selectable from the HTML column chooser.
.raw_variant_row_values <- function(row) {
  if (is.null(row) || !is.data.frame(row) || nrow(row) != 1L) {
    stop("'row' must be a one-row data.frame.", call. = FALSE)
  }
  keys <- names(row)
  values <- vapply(keys, function(nm) {
    value <- row[[nm]]
    if (is.factor(value)) value <- as.character(value)
    if (is.list(value)) {
      if (!length(value) || all(is.na(value))) return(NA_character_)
      return(paste(unlist(value, use.names = FALSE), collapse = ";"))
    }
    if (!length(value) || all(is.na(value))) return(NA_character_)
    paste(as.character(value), collapse = ";")
  }, character(1L), USE.NAMES = FALSE)
  list(keys = keys, values = values)
}

.table_payload <- function(df) {
  data <- lapply(df, function(col) {
    if (is.list(col)) col <- vapply(col, function(x) paste(x, collapse = ";"), character(1L))
    if (is.factor(col)) col <- as.character(col)
    # signif(, 15) keeps doubles exact (genomic positions, depths) while still
    # cleaning representation noise; signif(, 6) used to round e.g.
    # 28,057,643 -> 28,057,600 in the HTML table.
    if (is.numeric(col)) .arr(signif(col, 15L)) else .arr(col)
  })
  list(cols = .arr(names(df)), data = data, n = nrow(df))
}

# Genome-view data for one sample / gene (all calls; zero-length ones as markers).
# BUG FIX (NA subtype / NA VAF): the previous version used `row$ALUSubtype %||%
# "ALU"`, which passes NA through (only NULL triggers the fallback) and produced
# "NA%"-style labels for zero-VAF calls. NA/length-0 values are now normalised
# explicitly; hotspot/flag metadata is carried through so the details tab can
# draw badges without re-querying.
# Capture actionable diagnostics without coupling optional genome panels to
# the main report. Keep the message compact enough for batch logs while
# retaining the failing call and a bounded tail of the R call stack.
.html_error_diagnostic <- function(e, context = NULL, calls = sys.calls()) {
  call <- conditionCall(e)
  call_txt <- if (is.null(call)) "<no condition call>" else
    paste(deparse(call, width.cutoff = 120L), collapse = " ")
  stack_txt <- if (length(calls)) {
    paste(vapply(tail(calls, 10L), function(cl) {
      paste(deparse(cl, width.cutoff = 120L), collapse = "")
    }, character(1L)), collapse = " -> ")
  } else "<no call stack>"
  cls <- paste(class(e), collapse = "/")
  prefix <- if (is.null(context) || !nzchar(context)) "HTML report" else context
  sprintf("[%s] error=%s; class=%s; call=%s; stack=%s",
          prefix, conditionMessage(e), cls, call_txt, stack_txt)
}

.panel_payload <- function(sub_df, gene_config, bam_path, cov, sample, gene, mode, pos_col, ids) {
  pd <- tryCatch(
    .variant_plot_data(sub_df, gene_config, bam_path, cov, include_ptd = TRUE, pos_col = pos_col),
    error = function(e) {
      message(.html_error_diagnostic(
        e, context = sprintf("HTML report / genome payload %s/%s", sample, gene)
      ))
      NULL
    }
  )
  if (is.null(pd)) return(NULL)
  ax <- pd$axis
  cov_df <- pd$cov_df
  cov_n <- if (is.data.frame(cov_df)) nrow(cov_df) else 0L
  cov_depth <- if (is.data.frame(cov_df) && "depth" %in% names(cov_df)) {
    suppressWarnings(as.numeric(cov_df$depth))
  } else numeric(0)
  finite_cov <- is.finite(cov_depth)
  has_cov <- cov_n > 0L && any(finite_cov)
  message(sprintf("[HTML report] Genome panel %s/%s: coverage_points=%d; finite=%d; max_depth=%s",
                  sample, gene, cov_n, sum(finite_cov),
                  if (any(finite_cov)) as.character(max(cov_depth[finite_cov])) else "NA"))
  variants <- lapply(seq_len(nrow(pd$itd_df)), function(i) {
    row <- pd$itd_df[i, , drop = FALSE]
    vaf_raw <- suppressWarnings(as.numeric(row$AlleleFrequency))
    if (length(vaf_raw) != 1L || is.na(vaf_raw)) vaf_raw <- NA_real_
    vaf <- vaf_raw
    vaf_str <- if (is.na(vaf)) "N/A" else sprintf("%.1f%%", vaf * 100)
    label <- if (mode == "alu") {
      subtype <- if (is.null(row$ALUSubtype) || is.na(row$ALUSubtype)) "ALU" else as.character(row$ALUSubtype)
      sprintf("%s %s", subtype, vaf_str)
    } else {
      len_label <- if (pd$len[i] > 0) sprintf("%d bp", as.integer(pd$len[i])) else "PTD (clip)"
      extra <- if ("HGVS_cDNA" %in% names(row) && !is.na(row$HGVS_cDNA) && nzchar(row$HGVS_cDNA)) {
        paste0("  ", row$HGVS_cDNA)
      } else ""
      paste0(len_label, "  VAF ", vaf_str, extra)
    }
    hover <- if (mode == "alu") {
      f <- .variant_fields(row, mode = "alu")
      paste(names(f), f, sep = ": ", collapse = "<br>")
    } else {
      .variant_hover(row, i, vaf_str)
    }
    list(
      vid = ids[i], start = pd$pos[i], end = pd$dup_ends[i], len = pd$len[i],
      x0 = ax$transform_pos(pd$pos[i]), x1 = ax$transform_pos(pd$dup_ends[i]),
      vaf = vaf, label = label, hover = hover
    )
  })
  itd_cov <- list()
  if (is.data.frame(pd$itd_cov) && nrow(pd$itd_cov) > 0L) {
    itd_cov <- lapply(split(pd$itd_cov, pd$itd_cov$variant), function(d) {
      list(vid = ids[d$variant[1L]], x = .arr(d$disp_pos), y = .arr(d$count), pos = .arr(d$raw_idx))
    })
  }
  list(
    id = paste0("p_", make.names(sample), "_", make.names(gene)),
    sample = sample, gene = gene, chrom = gene_config$chrom,
    build = .genome_build_label(gene_config), strand = gene_config$strand %||% "",
    from = ax$disp_from, to = ax$disp_to,
    ticks = .arr(ax$exon_centers), tick_labels = .arr(ax$exon_labels),
    cov = if (has_cov) list(x = .arr(pd$cov_df$disp_pos), y = .arr(pd$cov_df$depth),
                            pos = .arr(pd$cov_df$raw_pos)) else NULL,
    exons = if (is.null(pd$exons)) list() else lapply(seq_len(nrow(pd$exons)), function(i) list(
      label = paste0("E", pd$exons$exon_idx[i]), start = pd$exons$start[i], end = pd$exons$end[i],
      x0 = ax$transform_pos(pd$exons$start[i]), x1 = ax$transform_pos(pd$exons$end[i]))),
    variants = variants,
    itd_cov = unname(itd_cov)
  )
}

#' Generate an HTML report for TALOS ITD/PTD or ALU results
#'
#' Writes \strong{one self-contained HTML file} (no R Markdown, no pandoc):
#' results are embedded as JSON and drawn client side, SPARKS style, with
#' tabs for an overview, a sortable / filterable results table with CSV
#' export, an interactive genome view per sample and gene (coverage, exons,
#' ITD read coverage, calls), per-variant details, the configuration and the
#' run information. Light / dark theme follows the system and can be toggled.
#'
#' @param result_df   Data frame from \code{\link{talos}} or \code{\link{detect_alu}}.
#' @param bam_paths   Named character vector of BAM paths (optional; used to
#'   draw coverage). Names must match \code{result_df$Sample}.
#' @param gene_configs Named list of resolved gene config lists.
#' @param output_file Path to the output HTML file.
#' @param title       Report title string.
#' @param show_config Include the per-gene configuration tab.
#' @param mode        \code{"itd"} (default) or \code{"alu"}.
#' @param coverages   Optional named list of precomputed genomic-coordinate-indexed
#'   coverage Rles (index 1 is genomic position 1), keyed
#'   \code{"<sample>_<gene>"} (avoids re-reading the BAM).
#' @param plotly_js   \code{"auto"} (inline plotly.js when the R plotly
#'   package is installed, else CDN), \code{"inline"} or \code{"cdn"}.
#' @return Invisibly, the normalised path to the written HTML file.
#' @export
talos_html_report <- function(result_df,
                              bam_paths = NULL,
                              gene_configs = NULL,
                              output_file = "TALOS_report.html",
                              title = "TALOS Analysis Report",
                              show_config = FALSE,
                              mode = c("itd", "alu"),
                              coverages = NULL,
                              plotly_js = c("auto", "inline", "cdn")) {
  mode <- match.arg(mode)
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Please install: jsonlite")
  if (is.null(result_df)) result_df <- data.frame()
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pos_col <- if (mode == "alu") "InsertionSite" else "GenomicPosition"

  # Results table (Genome build next to Gene, as before).
  df <- as.data.frame(result_df, stringsAsFactors = FALSE)
  get_config <- function(g) {
    if (is.null(gene_configs) || !is.list(gene_configs) || !length(gene_configs)) return(NULL)
    nm <- names(gene_configs)
    if (is.null(nm) || !g %in% nm) return(NULL)
    gene_configs[[match(g, nm)]]
  }
  build_of <- function(g) {
    cfg <- get_config(g)
    if (is.null(cfg)) "custom" else .genome_build_label(cfg)
  }
  if (nrow(df) && "Gene" %in% names(df)) {
    df$GenomeBuild <- vapply(as.character(df$Gene), build_of, character(1L))
    ord <- setdiff(names(df), "GenomeBuild")
    df <- df[, append(ord, "GenomeBuild", after = match("Gene", ord)), drop = FALSE]
  }

  # Variants (details + sequence), panels (genome view).
  groups <- if (nrow(df)) unique(df[, c("Sample", "Gene"), drop = FALSE]) else data.frame()
  variants <- list()
  panels <- list()
  for (gi in seq_len(nrow(groups))) {
    sn <- as.character(groups$Sample[gi]); gn <- as.character(groups$Gene[gi])
    rows <- which(df$Sample == sn & df$Gene == gn)
    sub_df <- df[rows, , drop = FALSE]
    gb <- build_of(gn)
    for (k in seq_along(rows)) {
      row <- sub_df[k, , drop = FALSE]
      fields <- .variant_fields(row, gb, mode)
      # Keep the complete raw result row alongside the curated clinical
      # details. The browser-side Details table uses the same column chooser
      # state as the Results table, so every result variable is selectable.
      raw_row <- .raw_variant_row_values(row)
      all_keys <- raw_row$keys
      all_values <- raw_row$values
      seq_col <- if (mode == "alu") "ALU_Sequence" else "ITD_Sequence"
      sq <- if (seq_col %in% names(row)) as.character(row[[seq_col]]) else NA_character_
      # Hotspot metadata carries through so the details tab can draw a badge
      # without re-reading the CSV/DB.
      hotspot_flag <- isTRUE(row$Hotspot)
      hotspot_name <- if (hotspot_flag && "HotspotName" %in% names(row) && !is.na(row$HotspotName)) {
        as.character(row$HotspotName)
      } else NULL
      variants[[length(variants) + 1L]] <- list(
        id = rows[k], sample = sn, gene = gn, n_in_group = k,
        title = if (mode == "alu") sprintf("ALU insertion %d @ %s (%s)", k, fields[["Insertion site"]], fields[["ALU subtype"]])
                else sprintf("Variant %d: %s", k, if (fields[["HGVS cDNA"]] != "N/A") fields[["HGVS cDNA"]] else
                  sprintf("%s bp @ %s", fields[["Length (bp)"]], fields[["Genomic Position"]])),
        keys = .arr(names(fields)), values = .arr(unname(fields)),
        all_keys = .arr(all_keys), all_values = .arr(all_values),
        sequence = if (!is.na(sq) && nzchar(sq)) sq else NULL,
        flag = isTRUE(row$MergeFragmentWarning) || isTRUE(row$ArtifactSuspect),
        hotspot = hotspot_flag,
        hotspot_name = hotspot_name
      )
    }
    cfg <- get_config(gn)
    if (is.null(cfg)) {
      message(sprintf("[HTML report] Genome panel skipped for %s/%s: no gene configuration.", sn, gn))
    }
    if (!is.null(cfg)) {
      key <- paste(sn, gn, sep = "_")
      bam <- if (!is.null(bam_paths) && !is.null(names(bam_paths)) && sn %in% names(bam_paths)) {
        bam_paths[[match(sn, names(bam_paths))]]
      } else NULL
      cov <- if (is.list(coverages) && !is.null(names(coverages)) && key %in% names(coverages)) {
        coverages[[match(key, names(coverages))]]
      } else NULL
      if (!is.null(cov) || (!is.null(bam) && length(bam) == 1L && !is.na(bam) && file.exists(bam))) {
        message(sprintf("[HTML report] Building genome panel %s/%s (BAM=%s; precomputed_coverage=%s)",
                        sn, gn,
                        if (!is.null(bam) && length(bam) == 1L && !is.na(bam) && file.exists(bam)) "available" else "unavailable",
                        if (!is.null(cov) && length(cov) > 0L) "available" else "unavailable"))
        # A malformed optional genome panel must not suppress the complete
        # results/details HTML report. Keep the panel error in the run log.
        pnl <- tryCatch(
          .panel_payload(sub_df, cfg, bam, cov, sn, gn, mode, pos_col, ids = rows),
          error = function(e) {
            message(.html_error_diagnostic(
              e, context = sprintf("HTML report / genome panel %s/%s", sn, gn)
            ))
            NULL
          }
        )
        if (!is.null(pnl)) panels[[length(panels) + 1L]] <- pnl
        else message(sprintf("[HTML report] Genome panel omitted for %s/%s; see preceding diagnostic.", sn, gn))
      } else {
        message(sprintf("[HTML report] Genome panel skipped for %s/%s: neither BAM nor precomputed coverage is available.", sn, gn))
      }
    }
  }

  configs <- if (isTRUE(show_config) && length(gene_configs)) {
    unname(lapply(names(gene_configs), function(gn) {
      cf <- .config_fields(gene_configs[[gn]])
      list(gene = gn, keys = .arr(names(cf)), values = .arr(unname(cf)))
    }))
  } else list()

  pl <- .plotly_script_tag(match.arg(plotly_js))
  is_alu <- mode == "alu"
  payload <- list(
    meta = list(
      title = title, generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      version = tryCatch(as.character(utils::packageVersion("TALOS")), error = function(e) "dev"),
      r_version = R.version.string, mode = mode,
      event = if (is_alu) "insertion" else "ITD", events = if (is_alu) "ALU insertions" else "ITDs / PTDs",
      plotly_inline = pl$inline
    ),
    hidden = .arr(intersect(.HTML_HIDDEN_COLS, names(df))),
    table = .table_payload(df),
    variants = variants,
    panels = panels,
    configs = configs,
    run = list(
      keys = .arr(c("Report generated", "TALOS version", "R version", "Total variants",
                    "Samples", "Genes", "Mode")),
      values = .arr(c(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                      tryCatch(as.character(utils::packageVersion("TALOS")), error = function(e) "dev"),
                      R.version.string, as.character(nrow(df)),
                      paste(unique(df$Sample), collapse = ", "), paste(unique(df$Gene), collapse = ", "),
                      toupper(mode)))
    )
  )
  message("[HTML report] Serializing report payload: variants=", length(variants),
          "; panels=", length(panels), "; rows=", nrow(df))
  json <- tryCatch(
    as.character(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", na = "null",
                                  digits = NA, force = TRUE)),
    error = function(e) {
      message(.html_error_diagnostic(e, context = "HTML report / JSON serialization"))
      stop(e)
    }
  )
  # Keep the JSON inert inside its <script> element.
  json <- gsub("<!--", "<\\!--", gsub("</", "<\\/", json, fixed = TRUE), fixed = TRUE)

  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    gsub(">", "&gt;", x, fixed = TRUE)
  }
  html <- paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", esc(title), "</title>\n",
    "<style>\n", .report_asset("talos_report.css"), "\n</style>\n",
    pl$tag, "\n</head>\n<body>\n",
    "<header class=\"top\"><div class=\"top-inner\">\n",
    "<div class=\"title-row\"><h1>", esc(title), "</h1>",
    "<span class=\"badge\">", toupper(mode), "</span>",
    "<span class=\"spacer\"></span>",
    "<button id=\"print-btn\" class=\"btn\" type=\"button\" aria-label=\"Print or save the report as PDF\">Print</button>",
    "<button id=\"theme-btn\" class=\"btn\" type=\"button\" aria-label=\"Toggle colour theme\">Theme</button></div>\n",
    "<div id=\"sample-pick\" class=\"sample-pick\" hidden></div>\n",
    "<nav id=\"tabs\" class=\"tabs\" role=\"tablist\" aria-label=\"Report sections\"></nav>\n",
    "</div></header>\n<main id=\"panels\"></main>\n",
    "<noscript><p style=\"padding:16px\">This report needs JavaScript enabled.</p></noscript>\n",
    "<script id=\"talos-data\" type=\"application/json\">", json, "</script>\n",
    "<script>\n", .talos_coverage_enhancement_js(), "\n</script>\n",
    "<script>\n", .report_asset("talos_report.js"), "\n</script>\n",
    "</body>\n</html>\n"
  )
  con <- tryCatch(file(output_file, open = "wb"), error = function(e) {
    message(.html_error_diagnostic(e, context = "HTML report / output file open"))
    stop(e)
  })
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  tryCatch(writeBin(charToRaw(enc2utf8(html)), con), error = function(e) {
    message(.html_error_diagnostic(e, context = "HTML report / output write"))
    stop(e)
  })
  message("[HTML report] HTML report written to: ", normalizePath(output_file))
  invisible(normalizePath(output_file))
}

.write_talos_output <- function(final_df, base_name, output_folder, sample_name,
                                gene_config, ref_dna, bam_path, write_vcf, plot,
                                html_report, verbose, add_config_to_report = FALSE,
                                precomputed_cov = NULL) {
  if (is.null(base_name)) return(invisible(NULL))
  sample_folder <- .ensure_dir(file.path(output_folder, sample_name))
  tsv_file <- file.path(sample_folder, paste0(base_name, ".tsv"))

  # Negative runs still produce a TSV (when there are columns to write) and an
  # HTML report with an explicit "No variants detected" state, so every
  # analysis has an audit artefact. VCF/PDF are skipped: they carry no
  # information for zero calls.
  if (ncol(final_df) > 0L) {
    .write_results_tsv(final_df, tsv_file)
    if (verbose) message(sprintf("Results written to: %s", tsv_file))
  }

  if (write_vcf && nrow(final_df) > 0L) {
    write_itd_vcf(itd_df = final_df, ref_dna = ref_dna, genomic_start = gene_config$genomic_start,
                  chrom = gene_config$chrom, sample_name = sample_name, genome_build = gene_config$build,
                  vcf_path = file.path(sample_folder, paste0(base_name, ".vcf")), overwrite = TRUE)
  }

  pdf_path <- NULL
  if (plot && nrow(final_df) > 0L) {
    pdf_path <- file.path(sample_folder, paste0(base_name, ".pdf"))
    plot_talos_report(itd_row = final_df, gene_config = gene_config, bam_path = bam_path,
                      sample_name = sample_name, output_pdf = pdf_path,
                      show_config = add_config_to_report, precomputed_cov = precomputed_cov)
  }

  report_path <- NULL
  if (html_report) {
    report_path <- file.path(sample_folder, paste0(base_name, "_report.html"))
    key <- paste(sample_name, gene_config$gene, sep = "_")
    message(sprintf("[HTML report] Starting report generation for %s/%s", sample_name, gene_config$gene))
    ok <- tryCatch({
      talos_html_report(
        result_df = final_df,
        bam_paths = if (is.null(bam_path)) NULL else stats::setNames(bam_path, sample_name),
        gene_configs = stats::setNames(list(gene_config), gene_config$gene),
        output_file = report_path,
        title = paste("TALOS Report -", gene_config$gene, sample_name),
        show_config = add_config_to_report,
        coverages = if (!is.null(precomputed_cov)) stats::setNames(list(precomputed_cov), key) else NULL
      )
      message(sprintf("[HTML report] Report generation completed for %s/%s", sample_name, gene_config$gene))
      TRUE
    }, error = function(e) {
      diagnostic <- .html_error_diagnostic(
        e, context = sprintf("HTML report / report generation %s/%s",
                             sample_name, gene_config$gene), calls = sys.calls()
      )
      warning("HTML report skipped: ", diagnostic, call. = FALSE)
      FALSE
    })
    if (!ok) report_path <- NULL
  }
  invisible(list(tsv = tsv_file, report = report_path, pdf = pdf_path))
}

.log_duration <- function(gene_name, sample_name, start_time, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  msg <- if (elapsed < 60) sprintf("%.1f sec", elapsed) else sprintf("%.1f min", elapsed / 60)
  message(sprintf("[TALOS] Analysis complete | Gene: %s | Sample: %s | Duration: %s",
                  gene_name, sample_name, msg))
}

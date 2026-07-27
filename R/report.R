# Helper: build padded exon display set and plot range
.build_plot_exon_context <- function(gene_config, flank_exons = 1L, flank_bp = 50L) {
  all_exons <- gene_config$all_exons
  target_exons <- gene_config$target_exons

  if (is.null(all_exons) || length(all_exons) == 0L) {
    return(list(
      exons_plot = NULL,
      plot_start = gene_config$genomic_start,
      plot_end = gene_config$genomic_end
    ))
  }

  all_idx <- S4Vectors::mcols(all_exons)$exon_idx
  if (is.null(all_idx) || length(all_idx) == 0L) {
    all_idx <- if (!is.null(S4Vectors::mcols(all_exons)$exon_number) && length(S4Vectors::mcols(all_exons)$exon_number) > 0L)
      as.integer(S4Vectors::mcols(all_exons)$exon_number)
    else seq_along(all_exons)
    S4Vectors::mcols(all_exons)$exon_idx <- all_idx
  }

  if (!is.null(target_exons) && length(target_exons) > 0L) {
    target_idx <- unique(S4Vectors::mcols(target_exons)$exon_idx)
    if (is.null(target_idx) || length(target_idx) == 0L)
      target_idx <- unique(S4Vectors::mcols(target_exons)$exon_number)
  } else {
    target_idx <- integer(0)
  }

  if (length(target_idx) > 0L) {
    min_t <- min(target_idx)
    max_t <- max(target_idx)
    left_idx  <- max(min(all_idx), min_t - flank_exons)
    right_idx <- min(max(all_idx), max_t + flank_exons)
    keep_idx  <- which(all_idx >= left_idx & all_idx <= right_idx)
    exons_plot <- all_exons[keep_idx]
  } else {
    exons_plot <- all_exons
  }

  exons_plot <- exons_plot[order(GenomicRanges::start(exons_plot))]
  plot_start <- max(1L, min(GenomicRanges::start(exons_plot)) - flank_bp)
  plot_end   <- max(GenomicRanges::end(exons_plot)) + flank_bp

  list(exons_plot = exons_plot, plot_start = plot_start, plot_end = plot_end)
}

.parse_itd_coverage_rle <- function(rle_str) {
  if (is.null(rle_str) || length(rle_str) == 0L || is.na(rle_str) || !nzchar(rle_str))
    return(numeric(0))
  parts <- strsplit(rle_str, ",", fixed = TRUE)[[1L]]
  vals <- numeric(0)
  lens <- integer(0)
  for (pt in parts) {
    sp <- strsplit(pt, ":", fixed = TRUE)[[1L]]
    if (length(sp) != 2L) next
    v <- suppressWarnings(as.numeric(sp[1L]))
    n <- suppressWarnings(as.integer(sp[2L]))
    if (is.na(v) || is.na(n) || n <= 0L) next
    vals <- c(vals, v)
    lens <- c(lens, n)
  }
  if (length(vals) == 0L) return(numeric(0))
  rep(vals, lens)
}

.build_compressed_axis <- function(exons_plot, raw_plot_start, raw_plot_end,
                                   intron_display = 40L,
                                   flank_display = 25L,
                                   outer_margin = 20L) {
  if (is.null(exons_plot) || length(exons_plot) == 0L) {
    raw_len <- max(1L, raw_plot_end - raw_plot_start + 1L)
    seg <- data.frame(
      raw_start = raw_plot_start, raw_end = raw_plot_end,
      disp_start = 1 + outer_margin, disp_end = raw_len + outer_margin,
      type = "linear", exon_label = NA_character_, stringsAsFactors = FALSE
    )
    transform_pos <- function(pos) {
      pos <- as.numeric(pos)
      seg$disp_start + (pos - seg$raw_start) / max(1, seg$raw_end - seg$raw_start) * (seg$disp_end - seg$disp_start)
    }
    return(list(
      segments = seg,
      transform_pos = transform_pos,
      disp_from = 0,
      disp_to = seg$disp_end + outer_margin,
      exon_centers = numeric(0),
      exon_labels = character(0)
    ))
  }

  exons_plot <- exons_plot[order(GenomicRanges::start(exons_plot))]
  exon_labels <- if (!is.null(S4Vectors::mcols(exons_plot)$exon_idx)) {
    paste0("E", S4Vectors::mcols(exons_plot)$exon_idx)
  } else if (!is.null(S4Vectors::mcols(exons_plot)$exon_number)) {
    paste0("E", S4Vectors::mcols(exons_plot)$exon_number)
  } else {
    paste0("E", seq_along(exons_plot))
  }

  segs <- list()
  disp_cursor <- 1 + outer_margin

  add_seg <- function(raw_start, raw_end, disp_len, type, exon_label = NA_character_) {
    if (raw_end < raw_start || disp_len <= 0L) return()
    segs[[length(segs) + 1L]] <<- data.frame(
      raw_start = as.numeric(raw_start),
      raw_end = as.numeric(raw_end),
      disp_start = as.numeric(disp_cursor),
      disp_end = as.numeric(disp_cursor + disp_len - 1L),
      type = type,
      exon_label = exon_label,
      stringsAsFactors = FALSE
    )
    disp_cursor <<- disp_cursor + disp_len
  }

  first_exon_start <- min(GenomicRanges::start(exons_plot))
  if (raw_plot_start < first_exon_start)
    add_seg(raw_plot_start, first_exon_start - 1L,
            min(flank_display, first_exon_start - raw_plot_start), "flank")

  for (i in seq_along(exons_plot)) {
    es <- GenomicRanges::start(exons_plot)[i]
    ee <- GenomicRanges::end(exons_plot)[i]
    add_seg(es, ee, ee - es + 1L, "exon", exon_labels[i])
    if (i < length(exons_plot)) {
      gs <- ee + 1L
      ge <- GenomicRanges::start(exons_plot)[i + 1L] - 1L
      if (ge >= gs)
        add_seg(gs, ge, min(intron_display, ge - gs + 1L), "intron")
    }
  }

  last_exon_end <- max(GenomicRanges::end(exons_plot))
  if (raw_plot_end > last_exon_end)
    add_seg(last_exon_end + 1L, raw_plot_end,
            min(flank_display, raw_plot_end - last_exon_end), "flank")

  segments <- do.call(rbind, segs)
  transform_pos <- function(pos) {
    pos <- as.numeric(pos)
    vapply(pos, function(pv) {
      idx <- which(pv >= segments$raw_start & pv <= segments$raw_end)[1L]
      if (is.na(idx)) {
        if (pv < min(segments$raw_start)) return(segments$disp_start[1L])
        return(segments$disp_end[nrow(segments)])
      }
      rs <- segments$raw_start[idx]; re <- segments$raw_end[idx]
      ds <- segments$disp_start[idx]; de <- segments$disp_end[idx]
      if (re <= rs) return(ds)
      ds + (pv - rs) / (re - rs) * (de - ds)
    }, numeric(1L))
  }

  exon_rows <- which(segments$type == "exon")
  exon_centers <- (segments$disp_start[exon_rows] + segments$disp_end[exon_rows]) / 2
  exon_labels2 <- segments$exon_label[exon_rows]

  list(
    segments = segments,
    transform_pos = transform_pos,
    disp_from = 0,
    disp_to = max(segments$disp_end) + outer_margin,
    exon_centers = exon_centers,
    exon_labels = exon_labels2
  )
}

.sample_compressed_coverage <- function(cov, axis_ctx, max_gap_points = 30L) {
  out <- vector("list", nrow(axis_ctx$segments))
  for (i in seq_len(nrow(axis_ctx$segments))) {
    seg <- axis_ctx$segments[i, ]
    raw_n <- as.integer(seg$raw_end - seg$raw_start + 1L)
    if (seg$type == "exon" || raw_n <= max_gap_points) {
      raw_pos <- seq.int(as.integer(seg$raw_start), as.integer(seg$raw_end))
    } else {
      raw_pos <- unique(round(seq(seg$raw_start, seg$raw_end, length.out = max_gap_points)))
      raw_pos <- as.integer(raw_pos)
    }
    disp_pos <- axis_ctx$transform_pos(raw_pos)
    out[[i]] <- data.frame(
      raw_pos = raw_pos,
      disp_pos = disp_pos,
      depth = as.numeric(cov[raw_pos]),
      type = seg$type,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

.build_itd_cov_display <- function(itd_df, dup_ends, axis_ctx) {
  rows <- vector("list", nrow(itd_df))
  idx <- 0L
  for (i in seq_len(nrow(itd_df))) {
    rle_str <- itd_df$ITDCoverageRLE[i]
    cov_vec <- .parse_itd_coverage_rle(rle_str)
    if (length(cov_vec) == 0L) next
    count <- as.numeric(cov_vec)
    x0 <- axis_ctx$transform_pos(itd_df$GenomicPosition[i])
    x1 <- axis_ctx$transform_pos(dup_ends[i])
    disp_pos <- if (length(count) == 1L) x0 else seq(x0, x1, length.out = length(count))
    idx <- idx + 1L
    rows[[idx]] <- data.frame(
      variant = i,
      disp_pos = disp_pos,
      count = count,
      raw_idx = seq_along(count),
      stringsAsFactors = FALSE
    )
  }
  if (idx == 0L) return(data.frame())
  do.call(rbind, rows[seq_len(idx)])
}

#' Generate a publication-ready PDF plot (Gviz + summary table)
#' Only the region restricted directly by targeted_exons context is shown.
#'
#' @param itd_row Data frame of TALOS results for one sample/gene.
#' @param gene_config Gene configuration list.
#' @param bam_path Path to BAM file.
#' @param sample_name Sample name.
#' @param output_pdf Output PDF file path.
#' @param width Plot dimension (inches).
#' @param height Plot dimension (inches).
#' @param show_config If TRUE, display configuration summary on a separate page.
#' @param precomputed_cov Optional precomputed coverage Rle (to avoid re-reading BAM).
#' @param include_ptd If TRUE, include zero-length (PTD) events (default FALSE).
#' @return Invisibly, output_pdf.
#' @export
plot_talos_report <- function(itd_row, gene_config, bam_path,
                               sample_name, output_pdf,
                               width = 14, height = 14,
                               show_config = FALSE,
                               precomputed_cov = NULL,
                               include_ptd = FALSE,
                               gap = 0.05) {
  itd_df <- itd_row
  if (!requireNamespace("Gviz", quietly = TRUE)) { warning("Gviz not installed – PDF plot skipped."); return(invisible(NULL)) }
  if (!requireNamespace("Rsamtools", quietly = TRUE)) stop("Rsamtools required")
  if (!requireNamespace("GenomicAlignments", quietly = TRUE)) stop("GenomicAlignments required")
  if (!requireNamespace("grid", quietly = TRUE)) stop("grid required")
  if (nrow(itd_df) == 0) { warning("No ITDs in itd_df."); return(invisible(NULL)) }

  chrom        <- gene_config$chrom
  genome_build <- ifelse(is.null(gene_config$build) || is.na(gene_config$build), "custom", gene_config$build)
  breakpoints  <- itd_df$GenomicPosition
  dup_lens     <- itd_df$Length
  dup_ends     <- breakpoints + dup_lens - 1L
  
  # Filter out zero-length (PTD) events unless include_ptd = TRUE
  keep <- if (include_ptd) rep(TRUE, nrow(itd_df)) else dup_lens > 0
  if (!any(keep)) {
    warning("No ITDs with positive length to plot (PTDs excluded). Set include_ptd=TRUE to include them.")
    return(invisible(NULL))
  }
  breakpoints <- breakpoints[keep]
  dup_lens    <- dup_lens[keep]
  dup_ends    <- dup_ends[keep]
  vafs         <- sprintf("%.1f", itd_df$AlleleFrequency[keep] * 100)
  orientations <- itd_df$Orientation[keep]
  orient_display <- ifelse(is.na(orientations) | orientations == "?", "", paste0(" | ", orientations))
  
  plot_ctx   <- .build_plot_exon_context(gene_config, flank_exons = 1L, flank_bp = 50L)
  exons_plot <- plot_ctx$exons_plot
  plot_start <- plot_ctx$plot_start
  plot_end   <- plot_ctx$plot_end

  if (is.null(precomputed_cov)) {
    param  <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges(chrom, IRanges::IRanges(plot_start, plot_end)))
    gal    <- GenomicAlignments::readGAlignments(bam_path, param = param)
    cov    <- GenomicAlignments::coverage(gal)[[chrom]]
  } else {
    cov <- precomputed_cov
  }
  axis_ctx <- .build_compressed_axis(exons_plot, plot_start, plot_end,
                                     intron_display = 40L,
                                     flank_display = 30L,
                                     outer_margin = 20L)
  cov_df <- .sample_compressed_coverage(cov, axis_ctx, max_gap_points = 30L)
  y_max  <- max(cov_df$depth, na.rm = TRUE); if (!is.finite(y_max) || y_max == 0) y_max <- 1L

  cov_track <- Gviz::DataTrack(
    start = cov_df$disp_pos, end = cov_df$disp_pos, data = cov_df$depth, chromosome = chrom,
    genome = genome_build, type = "histogram", name = "Coverage",
    col.histogram = "steelblue", fill.histogram = "#dfeefa",
    ylab = "Read depth", ylim = c(0, y_max * 1.05),
    grid = TRUE, lwd.grid = 0.4, col.grid = "#d0d0d0",
    cex.axis = 0.8, cex.title = 0.9, labelPos = "below",
    background.title = "transparent", col.border.title = NA,
    fontcolor.title = "black"
  )

  grtrack <- NULL
  if (!is.null(exons_plot) && length(exons_plot) > 0) {
    if (is.null(exons_plot$exon_idx) && !is.null(exons_plot$exon_number)) exons_plot$exon_idx <- exons_plot$exon_number
    else if (is.null(exons_plot$exon_idx)) exons_plot$exon_idx <- seq_along(exons_plot)

    exon_st <- axis_ctx$transform_pos(start(exons_plot))
    exon_en <- axis_ctx$transform_pos(end(exons_plot))
    exon_gr <- GenomicRanges::GRanges(
      seqnames = chrom, ranges = IRanges::IRanges(start = round(exon_st), end = round(exon_en)), strand = strand(exons_plot)
    )
    exon_gr$gene       <- gene_config$gene
    exon_gr$transcript <- paste0(gene_config$gene, "_tx")
    exon_gr$symbol     <- gene_config$gene
    exon_gr$exon       <- paste0("E", exons_plot$exon_idx)

    grtrack <- Gviz::GeneRegionTrack(
      exon_gr, genome = genome_build, chromosome = chrom,
      name = paste0(gene_config$gene, " exons (+ padding)"),
      fill = "#2ecc71", col = "#1a8a4a",
      transcriptAnnotation = "none", exonAnnotation = "exon",
      showExonId = TRUE, fontsize = 11, cex.feature = 1.0,
      fontcolor.feature = "black", shape = "box",
      stackHeight = 0.45,
      background.title = "transparent", col.border.title = NA,
      fontcolor.title = "black"
    )
  }

  itd_cov_df <- .build_itd_cov_display(itd_df[keep, , drop = FALSE], dup_ends, axis_ctx)
  itd_cov_track <- NULL
  if (nrow(itd_cov_df) > 0L) {
    itd_cov_track <- Gviz::DataTrack(
      start = itd_cov_df$disp_pos, end = itd_cov_df$disp_pos, data = itd_cov_df$count,
      chromosome = chrom, genome = genome_build, type = "histogram", name = "ITD cov reads",
      col.histogram = "#f39c12", fill.histogram = "#fdebd0", lwd = 2, ylim = c(0, max(itd_cov_df$count, na.rm = TRUE)), ylab = "ITD cov",
      background.title = "transparent", col.border.title = NA, fontcolor.title = "black"
    )
  }

  # Clean ITD labels: omit orientation segment when absent to avoid trailing " | "
  itd_labels <- ifelse(
    nchar(trimws(orient_display)) > 0,
    sprintf("%d bp | VAF %s%% | %s", dup_lens, vafs, trimws(orient_display)),
    sprintf("%d bp | VAF %s%%",      dup_lens, vafs)
  )

  # NB: without an explicit `group`, Gviz has no reliable way to tell the ITD
  # features apart, so with >1 ITD it was drawing every feature on the SAME
  # stacking row -- `stacking = "full"` alone did not separate them, and the
  # feature-id labels (drawn beside each box) rendered on top of one another
  # whenever two calls were close together. Giving every feature its own
  # group forces one row per ITD regardless of how close they are.
  itd_group <- factor(seq_along(breakpoints), levels = seq_along(breakpoints))

  itd_track  <- Gviz::AnnotationTrack(
    start = round(axis_ctx$transform_pos(breakpoints)), end = round(axis_ctx$transform_pos(dup_ends)), chromosome = chrom,
    genome = genome_build, name = ifelse(length(breakpoints) > 1, "ITDs", "ITD"),
    id = itd_labels,
    group = itd_group,
    stacking = "full",
    stackHeight = 0.45,
    showFeatureId = TRUE, just = "right",
    fill = "#f39c12", col = "#c0392b",
    fontcolor.feature = "#4a1c00", cex.feature = max(0.55, 0.8 - 0.04 * length(breakpoints)),
    cex.title = 0.9,
    rotation.title = 0, background.title = "transparent", col.border.title = NA, fontcolor.title = "black"
  )

  axis_track <- Gviz::GenomeAxisTrack(cex = 0.7, labelPos = "below", add35 = FALSE, add53 = FALSE, col = "black",
                                     at = axis_ctx$exon_centers, labels = axis_ctx$exon_labels)

  # Track sizes: each ITD now reliably renders on its own row (see itd_group
  # fix above) and stackHeight = 0.45 keeps each row's box thin, so the panel
  # no longer needs the old, oversized per-row allowance to avoid overlap.
  # ITD-coverage ("ITD cov reads") track size (0.75) is intentionally left
  # unchanged -- only the ITDs annotation track below it is being resized.
  itd_h <- max(1.0, 0.4 * length(breakpoints))
  if (!is.null(grtrack) && !is.null(itd_cov_track)) {
    track_list <- list(axis_track, grtrack, cov_track, itd_cov_track, itd_track)
    t_sizes <- c(0.5, 0.65, 2.2, 0.75, itd_h)
  } else if (!is.null(grtrack)) {
    track_list <- list(axis_track, grtrack, cov_track, itd_track)
    t_sizes <- c(0.5, 0.65, 2.2, itd_h)
  } else if (!is.null(itd_cov_track)) {
    track_list <- list(axis_track, cov_track, itd_cov_track, itd_track)
    t_sizes <- c(0.5, 2.2, 0.75, itd_h)
  } else {
    track_list <- list(axis_track, cov_track, itd_track)
    t_sizes <- c(0.5, 2.2, itd_h)
  }

  grDevices::pdf(output_pdf, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  main_title <- sprintf("Sample: %s  |  %s  |  %d ITD(s)  |  %s", sample_name, gene_config$gene, length(breakpoints), genome_build)
  subtitle   <- sprintf("Positions: %s  •  Lengths: %s  •  VAFs: %s%%", paste(breakpoints, collapse = ", "), paste(dup_lens, collapse = ", "), paste(vafs, collapse = ", "))

  Gviz::plotTracks(
    track_list, from = axis_ctx$disp_from, to = axis_ctx$disp_to,
    sizes = t_sizes, collapse = FALSE, main = main_title, cex.main = 1.2,
    fontcolor.title = "#2c3e50", col.axis = "black", cex.axis = 0.8,
    margins = c(18, 14, 6.5, 14), gap = gap
  )
  grid::grid.text(label = subtitle, x = grid::unit(0.5, "npc"), y = grid::unit(0.885, "npc"), just = c("center", "center"), gp = grid::gpar(fontsize = 9, col = "#2c3e50", fontface = "plain"))

  # Build summary table including new softclip metrics
  summary_df <- data.frame(
    Metric = c("Sample", "Gene", "Genome Build", "Genomic Position", "Length (bp)", "Length PE (bp)", "Length Ext (bp)", "VAF (%)", 
               "Supporting Reads", "Wildtype Reads", "Depth at Breakpoint", "HGVS cDNA", "HGVS Protein", 
               "Region", "Exon Number", "Orientation", "Hotspot", "Strand Bias", "Mean MAPQ", 
               "Support Consistency", "RefMatch (%)", "ITD Read Coverage (%)", "Sequence Source",
               "Merged From N Fragment(s)",
               "Left Softclip Count", "Right Softclip Count", "Left Softclip % (Support)", "Right Softclip % (Support)",
               "Left Softclip % (WT)", "Right Softclip % (WT)",
               "PE Softclip Support", "PE Event-size Pairs", "PE Long-span Pairs", "PE Orientation"),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(itd_df[keep, , drop = FALSE]))) {
    row <- itd_df[keep, ][i, ]
    pe_str <- if (!is.null(row$LengthPE) && !is.na(row$LengthPE))
      sprintf("%d bp (%d pairs)", as.integer(row$LengthPE),
              as.integer(if (!is.null(row$LengthPE_NSpanning)) row$LengthPE_NSpanning else 0L))
    else "N/A"
    summary_df[[paste0("ITD_", i)]] <- c(
      as.character(sample_name), as.character(row$Gene), genome_build, as.character(row$GenomicPosition), 
      as.character(row$Length), pe_str,
      ifelse(is.null(row$LengthExt) || is.na(row$LengthExt), "N/A", as.character(row$LengthExt)),
      sprintf("%.1f", row$AlleleFrequency * 100), 
      as.character(row$SupportingReads), as.character(round(row$WildtypeReads)), as.character(row$DepthAtBreakpoint), 
      ifelse(is.na(row$HGVS_cDNA), "N/A", row$HGVS_cDNA), 
      ifelse(is.na(row$HGVS_Protein), "N/A", row$HGVS_Protein), 
      ifelse(is.na(row$Region), "N/A", row$Region), 
      ifelse(is.na(row$ExonNumber), "N/A", as.character(row$ExonNumber)), 
      ifelse(is.na(row$Orientation), "N/A", row$Orientation), 
      ifelse(isTRUE(row$Hotspot), as.character(row$HotspotName), "No"), 
      ifelse(is.na(row$StrandBias), "N/A", sprintf("%.3f", row$StrandBias)), 
      ifelse(is.na(row$MeanSupportMAPQ), "N/A", sprintf("%.1f", row$MeanSupportMAPQ)), 
      ifelse(is.na(row$SupportConsistency), "N/A", paste0(round(row$SupportConsistency), "%")), 
      ifelse(is.na(row$RefMatch_Observed), "N/A", paste0(round(row$RefMatch_Observed, 1), "%")), 
      ifelse(is.na(row$ITDReadCoverage), "N/A", paste0(round(row$ITDReadCoverage, 1), "%")), 
      ifelse(is.na(row$SequenceSource), "N/A", row$SequenceSource),
      {
        n_merged <- if (is.null(row$MergedFromN) || is.na(row$MergedFromN)) 1L else as.integer(row$MergedFromN)
        flagged  <- isTRUE(row$MergeFragmentWarning)
        if (n_merged <= 1L) "1 (not merged)"
        else sprintf("%d%s", n_merged, if (flagged) " (!! unusually high !!)" else "")
      },
      # NOTE: these were previously displayed as "0" whenever NA, which
      # misrepresented "not computed for this row" (e.g. a genuinely
      # ambiguous multi-fragment merge, prior to the merge-aggregation fix)
      # as "confirmed zero soft-clip reads". Softclip counts are always a
      # real (possibly zero) integer for any non-merged call, and are now
      # summed rather than nulled for merged calls too, so NA should no
      # longer occur in practice -- but if it ever does, "N/A" is the
      # honest label, matching every sibling field in this table.
      ifelse(is.na(row$LeftSoftclipCount), "N/A", as.character(row$LeftSoftclipCount)),
      ifelse(is.na(row$RightSoftclipCount), "N/A", as.character(row$RightSoftclipCount)),
      ifelse(is.na(row$LeftSoftclipPctSupport), "N/A", sprintf("%.1f", row$LeftSoftclipPctSupport)),
      ifelse(is.na(row$RightSoftclipPctSupport), "N/A", sprintf("%.1f", row$RightSoftclipPctSupport)),
      ifelse(is.na(row$LeftSoftclipPctWT), "N/A", sprintf("%.1f", row$LeftSoftclipPctWT)),
      ifelse(is.na(row$RightSoftclipPctWT), "N/A", sprintf("%.1f", row$RightSoftclipPctWT)),
      ifelse(is.null(row$PESoftclipSupport) || is.na(row$PESoftclipSupport), "N/A", as.character(row$PESoftclipSupport)),
      ifelse(is.null(row$PESoftclipEventPairs) || is.na(row$PESoftclipEventPairs), "N/A", as.character(row$PESoftclipEventPairs)),
      ifelse(is.null(row$PESoftclipLongPairs) || is.na(row$PESoftclipLongPairs), "N/A", as.character(row$PESoftclipLongPairs)),
      ifelse(is.null(row$PEOrientationDominant) || is.na(row$PEOrientationDominant), "N/A",
             sprintf("%s (FR=%s RF=%s FF=%s RR=%s)", row$PEOrientationDominant,
                     row$PEOrientationFR, row$PEOrientationRF, row$PEOrientationFF, row$PEOrientationRR))
    )
  }

  grid::grid.newpage()
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    n_rows   <- nrow(summary_df)
    row_fill <- rep(c("white", "#f0f4f8"), length.out = n_rows)
    tbl <- gridExtra::tableGrob(summary_df, rows = NULL, theme = gridExtra::ttheme_minimal(
      core = list(fg_params = list(cex = 0.7, fontface = "plain"), 
                  bg_params = list(fill = row_fill, col = NA)), 
      colhead = list(fg_params = list(cex = 0.8, fontface = "bold", col = "white"), 
                     bg_params = list(fill = "#2c3e50", col = NA))))
    grid::grid.text(sprintf("Variant Summary – Sample: %s | %s (%d variants)", sample_name, gene_config$gene, nrow(itd_df[keep, , drop = FALSE])), 
                    x = 0.5, y = 0.91, gp = grid::gpar(fontsize = 13, fontface = "bold", col = "#2c3e50"))
    vp <- grid::viewport(x = 0.5, y = 0.49, width = 0.90, height = 0.80)
    grid::pushViewport(vp); grid::grid.draw(tbl); grid::popViewport()
  } else {
    grid::grid.text(paste(capture.output(print(summary_df)), collapse = "\n"), x = 0.05, y = 0.95, just = "left", gp = grid::gpar(cex = 0.7, fontfamily = "mono"))
  }

  if (show_config) {
    grid::grid.newpage()
    grid::grid.text("Analysis Configuration", x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 14, fontface = "bold", col = "#2c3e50"))
    ypos <- 0.85
    config_lines <- c(
      paste("Gene:", gene_config$gene), paste("Build:", genome_build), paste("Chromosome:", gene_config$chrom), 
      paste("Strand:", gene_config$strand), paste("Transcript reference:", ifelse(is.na(gene_config$transcript), "none", gene_config$transcript)), 
      paste("CDS offset:", gene_config$cds_offset), paste("Genomic window:", gene_config$genomic_start, "-", gene_config$genomic_end), 
      paste("Targeted exons:", paste(unique(gene_config$target_exons$exon_idx), collapse = ", ")), 
      paste("Settings:"), 
      paste("  min_support =", gene_config$gene_settings$min_support %||% "default"), 
      paste("  cluster_tolerance =", gene_config$gene_settings$cluster_tolerance %||% "default"), 
      paste("  min_mapq =", gene_config$gene_settings$min_mapq %||% "default")
    )
    for (line in config_lines) { grid::grid.text(line, x = 0.05, y = ypos, just = "left", gp = grid::gpar(fontsize = 10, fontfamily = "mono")); ypos <- ypos - 0.04 }
  }

  message("PDF plot saved to: ", output_pdf)
  invisible(output_pdf)
}

#' Generate an interactive plotly visualisation for TALOS results
#' Only the region restricted by targeted exons is shown.
#'
#' @param itd_df Data frame of TALOS results.
#' @param gene_config Gene configuration list.
#' @param bam_path Path to indexed BAM file.
#' @param sample_name Sample identifier.
#' @param show_config Add a config annotation string.
#' @param precomputed_cov Optional precomputed coverage Rle (to avoid re-reading BAM).
#' @param include_ptd If TRUE, include zero-length (PTD) events (default FALSE).
#' @return A plotly htmlwidget.
#' @export
plot_talos_interactive <- function(itd_df, gene_config, bam_path,
                                    sample_name = NULL,
                                    show_config = FALSE,
                                    precomputed_cov = NULL,
                                    include_ptd = FALSE) {
  if (!requireNamespace("plotly", quietly = TRUE)) { warning("plotly not installed."); return(invisible(NULL)) }
  if (!requireNamespace("Rsamtools", quietly = TRUE)) stop("Rsamtools required")
  if (!requireNamespace("GenomicAlignments", quietly = TRUE)) stop("GenomicAlignments required")
  if (nrow(itd_df) == 0) return(invisible(NULL))

  chrom        <- gene_config$chrom
  genome_build <- ifelse(is.null(gene_config$build) || is.na(gene_config$build), "custom", gene_config$build)
  breakpoints  <- itd_df$GenomicPosition
  dup_lens     <- itd_df$Length
  dup_ends     <- breakpoints + dup_lens - 1L

  # Filter out zero-length (PTD) events unless include_ptd = TRUE
  valid <- if (include_ptd) rep(TRUE, nrow(itd_df)) else dup_lens > 0
  if (!any(valid)) {
    warning("No ITDs with positive length to plot interactively (PTDs excluded). Set include_ptd=TRUE to include them.")
    return(plotly::plot_ly())
  }
  itd_df       <- itd_df[valid, , drop = FALSE]
  breakpoints  <- breakpoints[valid]
  dup_lens     <- dup_lens[valid]
  dup_ends     <- dup_ends[valid]

  plot_ctx   <- .build_plot_exon_context(gene_config, flank_exons = 1L, flank_bp = 50L)
  exons_plot <- plot_ctx$exons_plot
  plot_start <- plot_ctx$plot_start
  plot_end   <- plot_ctx$plot_end

  if (is.null(precomputed_cov)) {
    param <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges(chrom, IRanges::IRanges(plot_start, plot_end)))
    gal   <- GenomicAlignments::readGAlignments(bam_path, param = param)
    cov   <- GenomicAlignments::coverage(gal)[[chrom]]
  } else {
    cov <- precomputed_cov
  }
  axis_ctx <- .build_compressed_axis(exons_plot, plot_start, plot_end,
                                     intron_display = 40L,
                                     flank_display = 30L,
                                     outer_margin = 20L)
  cov_df <- .sample_compressed_coverage(cov, axis_ctx, max_gap_points = 30L)

  p_cov <- plotly::plot_ly(
    x = cov_df$disp_pos, y = cov_df$depth, type = "scatter", mode = "lines",
    fill = "tozeroy", name = "Coverage",
    line = list(color = "#4a90d9", width = 1.2),
    fillcolor = "rgba(74,144,217,0.22)",
    hovertemplate = "Depth: %{y}<extra></extra>"
  )
  p_cov <- plotly::layout(p_cov, yaxis = list(title = "Depth", gridcolor = "#e8e8e8"), xaxis = list(title = "", showticklabels = FALSE, range = c(axis_ctx$disp_from, axis_ctx$disp_to)), plot_bgcolor = "white", paper_bgcolor = "white")

  p_exon <- plotly::plot_ly()
  p_exon <- plotly::layout(p_exon, xaxis = list(title = "", showticklabels = FALSE, range = c(axis_ctx$disp_from, axis_ctx$disp_to), rangeslider = list(visible = FALSE)), yaxis = list(title = "Exons", showticklabels = FALSE, range = c(0, 1), fixedrange = TRUE), plot_bgcolor = "#f8fafc", paper_bgcolor = "white")

  if (!is.null(exons_plot) && length(exons_plot) > 0) {
    exon_idx_vec <- if (!is.null(exons_plot$exon_idx)) exons_plot$exon_idx else if (!is.null(exons_plot$exon_number)) exons_plot$exon_number else seq_along(exons_plot)

    for (i in seq_along(exons_plot)) {
      es <- axis_ctx$transform_pos(start(exons_plot[i])); ee <- axis_ctx$transform_pos(end(exons_plot[i]))
      p_exon <- p_exon |> plotly::add_trace(x = c(es, ee, ee, es, es), y = c(0.1, 0.1, 0.9, 0.9, 0.1), type = "scatter", mode = "lines", fill = "toself", fillcolor = "rgba(46,204,113,0.7)", line = list(color = "#1a8a4a", width = 1), name = paste0("Exon ", exon_idx_vec[i]), hovertemplate = sprintf("<b>Exon %s</b><br>%d–%d<extra></extra>", exon_idx_vec[i], start(exons_plot[i]), end(exons_plot[i])), showlegend = FALSE) |> plotly::add_annotations(x = (es + ee) / 2, y = 0.5, xref = "x", yref = "y", text = paste0("E", exon_idx_vec[i]), showarrow = FALSE, font = list(size = 10, color = "#1a5e30"))
    }
  }

  pal <- c("#e67e22", "#e74c3c", "#9b59b6", "#1abc9c", "#3498db", "#f39c12", "#c0392b", "#8e44ad", "#16a085", "#2980b9")
  # Adaptive vertical step: tighten spacing when many ITDs so they all fit in the panel
  y_step  <- max(0.10, min(0.18, 0.85 / max(1L, nrow(itd_df))))
  y_range_max <- max(1.1, nrow(itd_df) * y_step + 0.30)
  itd_cov_df <- .build_itd_cov_display(itd_df, dup_ends, axis_ctx)
  p_itd_cov <- plotly::plot_ly()
  if (nrow(itd_cov_df) > 0) {
    for (vid in unique(itd_cov_df$variant)) {
      sub <- itd_cov_df[itd_cov_df$variant == vid, , drop = FALSE]
      col_i <- pal[((vid - 1L) %% length(pal)) + 1L]
      p_itd_cov <- p_itd_cov |> plotly::add_trace(x = sub$disp_pos, y = sub$count, type = "scatter", mode = "lines", fill = "tozeroy", fillcolor = "rgba(243,156,18,0.22)", line = list(color = col_i, width = 2), name = paste0("ITD cov ", vid), hovertemplate = "ITD pos: %{text}<br>Coverage reads: %{y}<extra></extra>", text = sub$raw_idx, showlegend = FALSE)
    }
  }
  itd_cov_ymax <- if (nrow(itd_cov_df) > 0) max(itd_cov_df$count, na.rm = TRUE) else 1
  p_itd_cov <- plotly::layout(p_itd_cov, xaxis = list(title = "", showticklabels = FALSE, range = c(axis_ctx$disp_from, axis_ctx$disp_to), rangeslider = list(visible = FALSE)), yaxis = list(title = "ITD cov", range = c(0, max(1, itd_cov_ymax)), fixedrange = TRUE), plot_bgcolor = "white", paper_bgcolor = "white")

  p_itd <- plotly::plot_ly()

  for (i in seq_len(nrow(itd_df))) {
    row      <- itd_df[i, ]
    vaf_str  <- sprintf("%.1f%%", row$AlleleFrequency * 100)
    col_i    <- pal[((i - 1L) %% length(pal)) + 1L]
    y_pos    <- 0.05 + (i - 1L) * y_step
    tip_text <- paste0("<b>ITD ", i, "</b><br>",
                       "Breakpoint: ", row$GenomicPosition, "<br>",
                       "Length: ", row$Length, " bp<br>",
                       "VAF: ", vaf_str, "<br>",
                       "HGVS cDNA: ", ifelse(is.na(row$HGVS_cDNA), "N/A", row$HGVS_cDNA), "<br>",
                       "HGVS Protein: ", ifelse(is.na(row$HGVS_Protein), "N/A", row$HGVS_Protein), "<br>",
                       "Region: ", ifelse(is.na(row$Region), "N/A", row$Region), "<br>",
                       "Support: ", row$SupportingReads, " reads<br>",
                       "RefMatch: ", ifelse(is.na(row$RefMatch_Observed), "N/A", paste0(round(row$RefMatch_Observed, 1), "%")), "<br>",
                       "ITD Cov: ", ifelse(is.na(row$ITDReadCoverage), "N/A", paste0(round(row$ITDReadCoverage, 1), "%")), "<br>",
                       "Orientation: ", ifelse(is.na(row$Orientation), "?", row$Orientation), "<br>",
                       "Left Softclip: ", row$LeftSoftclipCount, " (", ifelse(is.na(row$LeftSoftclipPctSupport), "0", round(row$LeftSoftclipPctSupport,1)), "% of support)<br>",
                       "Right Softclip: ", row$RightSoftclipCount, " (", ifelse(is.na(row$RightSoftclipPctSupport), "0", round(row$RightSoftclipPctSupport,1)), "% of support)")
    label <- sprintf("%d bp  VAF %s", row$Length, vaf_str)
    if (!is.na(row$HGVS_cDNA)) label <- paste0(label, "  ", row$HGVS_cDNA)

    x0 <- axis_ctx$transform_pos(row$GenomicPosition)
    x1 <- axis_ctx$transform_pos(dup_ends[i])
    p_itd <- p_itd |> plotly::add_trace(x = c(x0, x1), y = c(y_pos, y_pos), type = "scatter", mode = "lines", line = list(color = col_i, width = 14), name = label, text = tip_text, hoverinfo = "text", showlegend = TRUE) |> plotly::add_trace(x = x0, y = y_pos, type = "scatter", mode = "markers", marker = list(symbol = "line-ns", size = 14, color = col_i, line = list(color = "white", width = 2)), hoverinfo = "skip", showlegend = FALSE)
  }

  p_itd <- plotly::layout(p_itd, xaxis = list(title = paste0(chrom, " (compressed ", genome_build, ")"), gridcolor = "#e8e8e8", range = c(axis_ctx$disp_from, axis_ctx$disp_to), rangeslider = list(visible = FALSE)), yaxis = list(title = "ITDs", showticklabels = FALSE, range = c(0, y_range_max), fixedrange = TRUE), legend = list(orientation = "h", x = 0, y = -0.18, font = list(size = 9), tracegroupgap = 0), plot_bgcolor = "white", paper_bgcolor = "white")

  title_text <- if (!is.null(sample_name)) paste0("TALOS • ", sample_name, "  |  ", gene_config$gene, "  |  ", nrow(itd_df), " ITD(s)  |  ", genome_build) else paste0("TALOS • ", gene_config$gene, "  |  ", nrow(itd_df), " ITD(s)  |  ", genome_build)

  # Heights: coverage reduced (was 0.42), ITD-cov reduced (was 0.14),
  # ITD panel expanded (was 0.30) so all calls and their labels have room.
  fig <- plotly::subplot(p_cov, p_exon, p_itd_cov, p_itd, nrows = 4, shareX = TRUE, heights = c(0.26, 0.10, 0.09, 0.55))
  fig <- plotly::layout(fig, title = list(text = title_text, x = 0.5, font = list(size = 14)), hovermode = "x unified", showlegend = TRUE, margin = list(t = 60, l = 60, r = 40, b = 75))
  
  if (show_config) {
    config_text <- paste("Config:", paste("Targeted exons:", paste(unique(gene_config$target_exons$exon_idx), collapse = ",")), paste("Transcript:", ifelse(is.na(gene_config$transcript), "none", gene_config$transcript)), paste("CDS offset:", gene_config$cds_offset), sep = " | ")
    fig <- plotly::layout(fig, annotations = list(text = config_text, x = 0, y = 1.05, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size = 9)))
  }

  fig
}


# ============================================================================
# Patch for report.R – talos_html_report() rewrite
# ============================================================================
# Code-review changes applied:
#   [CR-2] Remove `...` forwarded to rmarkdown::render(); use explicit params.
#   [CR-8] Replace fragile inline Rmd generation with external template file
#          installed at inst/rmd/TALOS_report.Rmd.
#   [CR-9] Add `mode` parameter to support ALU reports via the same template.
#   [CR-10] Validate output path before attempting to write.
# ============================================================================

#' Generate an HTML report for TALOS ITD/PTD or ALU results
#'
#' Renders the external template at \code{inst/rmd/TALOS_report.Rmd}.
#' For development (package not installed), falls back to
#' \code{./inst/rmd/TALOS_report.Rmd} relative to the working directory.
#'
#' @param result_df   Data frame from \code{\link{talos}} or \code{\link{detect_alu}}.
#' @param bam_paths   Named character vector of BAM paths (optional; used to
#'   generate interactive plotly widgets). Names must match \code{result_df$Sample}.
#' @param gene_configs Named list of resolved gene config lists.
#' @param output_file Path to the output HTML file.
#' @param title       Report title string.
#' @param show_config Include per-gene configuration block in the report.
#' @param mode        \code{"itd"} (default) or \code{"alu"}.
#' @return Invisibly, the normalised path to the written HTML file.
#' @export
talos_html_report <- function(result_df,
                               bam_paths    = NULL,
                               gene_configs = NULL,
                               output_file  = "TALOS_report.html",
                               title        = "TALOS Analysis Report",
                               show_config  = FALSE,
                               mode         = c("itd", "alu")) {

  mode <- match.arg(mode)

  # ── Dependency check ──────────────────────────────────────────────────
  needs <- c("rmarkdown", "DT")
  missing_pkg <- needs[!vapply(needs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing_pkg) > 0L)
    stop("Please install: ", paste(missing_pkg, collapse = ", "))

  # ── Resolve template path ─────────────────────────────────────────────
  template_path <- system.file("rmd", "TALOS_report.Rmd", package = "TALOS")
  if (!nzchar(template_path) || !file.exists(template_path)) {
    # Development fallback: look relative to the package root
    template_path <- file.path("inst", "rmd", "TALOS_report.Rmd")
  }
  if (!file.exists(template_path))
    stop(
      "TALOS_report.Rmd not found. Expected at:\n  ",
      system.file("rmd", "TALOS_report.Rmd", package = "TALOS"),
      "\nRe-install the package or ensure inst/rmd/TALOS_report.Rmd exists."
    )

  # ── Resolve output path ───────────────────────────────────────────────
  if (!grepl("^[A-Za-z]:|^/", output_file))
    output_file <- file.path(getwd(), output_file)
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ── Build plotly widgets ──────────────────────────────────────────────
  plot_widgets <- list()
  has_plotly   <- requireNamespace("plotly", quietly = TRUE)

  if (mode == "itd" && has_plotly &&
      !is.null(bam_paths) && !is.null(gene_configs)) {
    for (sn in names(bam_paths)) {
      bam <- bam_paths[[sn]]
      if (!file.exists(bam)) next
      for (gn in names(gene_configs)) {
        sub_df <- result_df[result_df$Sample == sn &
                              result_df$Gene   == gn, , drop = FALSE]
        if (nrow(sub_df) == 0L) next
        key <- paste(sn, gn, sep = "_")
        plot_widgets[[key]] <- tryCatch(
          plot_talos_interactive(sub_df, gene_configs[[gn]], bam, sn,
                                  show_config = show_config),
          error = function(e) {
            message(sprintf(
              "[HTML report] Could not generate plot for %s/%s: %s",
              sn, gn, conditionMessage(e)
            ))
            NULL
          }
        )
      }
    }
  }

  # ── Render ────────────────────────────────────────────────────────────
  rmarkdown::render(
    input       = template_path,
    output_file = output_file,
    params      = list(
      result_df    = result_df,
      plot_widgets = plot_widgets,
      gene_configs = gene_configs,
      has_plotly   = has_plotly && length(plot_widgets) > 0L,
      show_config  = show_config,
      title        = title,
      mode         = mode
    ),
    quiet = TRUE
  )

  message("HTML report written to: ", normalizePath(output_file))
  invisible(output_file)
}

# ============================================================================
# .write_talos_output – fixed with vapply and length sanity check
# ============================================================================
# NOTE (code review): this function used to be defined TWICE in this file.
# The earlier definition (TSV + VCF + PDF + HTML only, no plotly widget file)
# was silently shadowed by this later one, since R keeps the last definition
# in a source()'d file. That made the earlier copy dead code — any future
# edit made to it would have had zero effect, which is a maintenance trap.
# The duplicate has been removed; this is now the single source of truth.
.write_talos_output <- function(final_df, base_name, output_folder, sample_name,
                                gene_config, ref_dna, bam_path, write_vcf, plot,
                                html_report, verbose, add_config_to_report = FALSE) {
  if (is.null(base_name) || nrow(final_df) == 0) return(invisible(NULL))
  sample_folder <- file.path(output_folder, sample_name)
  if (!dir.exists(sample_folder)) dir.create(sample_folder, recursive = TRUE, showWarnings = FALSE)

  tsv_file <- file.path(sample_folder, paste0(base_name, ".tsv"))
  
  # ---- SAFE conversion of list columns to character using vapply ----
  for (col in names(final_df)) {
    if (is.list(final_df[[col]])) {
      # vapply ensures each element becomes a single string; if any result is not length 1, it errors
      converted <- vapply(final_df[[col]], function(x) {
        if (length(x) == 0 || all(is.na(x))) {
          NA_character_
        } else {
          # Collapse multi‑element vectors, handling numeric/logical/character
          paste(x, collapse = ";")
        }
      }, FUN.VALUE = character(1), USE.NAMES = FALSE)
      
      # Sanity check: length must equal number of rows
      if (length(converted) != nrow(final_df)) {
        stop(sprintf("Column '%s' conversion produced %d rows, expected %d. Check list elements.", 
                     col, length(converted), nrow(final_df)))
      }
      final_df[[col]] <- converted
    }
  }
  
  # Write TSV without row names, using NA placeholder
  write.table(final_df, file = tsv_file, sep = "\t", quote = FALSE,
              row.names = FALSE, na = ".")
  if (verbose) message(sprintf("Results written to: %s", tsv_file))

  if (write_vcf) {
    vcf_file <- file.path(sample_folder, paste0(base_name, ".vcf"))
    write_itd_vcf(itd_df = final_df, ref_dna = ref_dna, genomic_start = gene_config$genomic_start, chrom = gene_config$chrom, sample_name = sample_name, genome_build = gene_config$build, vcf_path = vcf_file, overwrite = TRUE)
  }

  pdf_path <- NULL
  if (plot && nrow(final_df) > 0) {
    pdf_file <- file.path(sample_folder, paste0(base_name, ".pdf"))
    plot_talos_report(itd_row = final_df, gene_config = gene_config, bam_path = bam_path, sample_name = sample_name, output_pdf = pdf_file, show_config = add_config_to_report)
    pdf_path <- pdf_file
  }

  # Single unified HTML report: talos_html_report() embeds the
  # plot_talos_interactive() plotly widget directly alongside the full
  # metrics table, so there is exactly one HTML deliverable per sample
  # (no separate standalone widget-only file).
  report_path <- NULL
  if (html_report && nrow(final_df) > 0) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) { warning("rmarkdown not installed. Cannot generate HTML report.") } 
    else {
      report_file <- file.path(sample_folder, paste0(base_name, "_report.html"))
      bam_vec     <- setNames(bam_path, sample_name)
      config_list <- setNames(list(gene_config), gene_config$gene)
      suppressWarnings(talos_html_report(result_df = final_df, bam_paths = bam_vec, gene_configs = config_list, output_file = report_file, title = paste("TALOS Report –", gene_config$gene, sample_name), show_config = add_config_to_report))
      if (verbose) message(sprintf("HTML report written to: %s", report_file))
      report_path <- report_file
    }
  }
  invisible(list(tsv = tsv_file, report = report_path, pdf = pdf_path))
}

.log_duration <- function(gene_name, sample_name, start_time, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  msg     <- if (elapsed < 60) sprintf("%.1f sec", elapsed) else sprintf("%.1f min", elapsed / 60)
  message(sprintf("[TALOS] Analysis complete | Gene: %s | Sample: %s | Duration: %s", gene_name, sample_name, msg))
}
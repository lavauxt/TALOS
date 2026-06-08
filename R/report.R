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
                               width = 12, height = 9,
                               show_config = FALSE,
                               precomputed_cov = NULL,
                               include_ptd = FALSE) {
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
  
  plot_start <- gene_config$genomic_start
  plot_end   <- gene_config$genomic_end

  if (is.null(precomputed_cov)) {
    param  <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges(chrom, IRanges::IRanges(plot_start, plot_end)))
    gal    <- GenomicAlignments::readGAlignments(bam_path, param = param)
    cov    <- GenomicAlignments::coverage(gal)[[chrom]]
  } else {
    cov <- precomputed_cov
  }
  pos    <- seq(plot_start, plot_end)
  depth  <- as.numeric(cov[pos])
  y_max  <- max(depth, na.rm = TRUE); if (y_max == 0) y_max <- 1L

  cov_track <- Gviz::DataTrack(
    start = pos, end = pos, data = depth, chromosome = chrom,
    genome = genome_build, type = "histogram", name = "Coverage",
    col.histogram = "steelblue", fill.histogram = "#4a90d9",
    ylab = "Read depth", ylim = c(0, y_max * 1.05),
    grid = TRUE, lwd.grid = 0.4, col.grid = "#d0d0d0",
    cex.axis = 0.8, cex.title = 0.9, labelPos = "below",
    background.title = "transparent", col.border.title = NA,
    fontcolor.title = "black"
  )

  grtrack <- NULL
  if (!is.null(gene_config$target_exons) && length(gene_config$target_exons) > 0) {
    exons_full <- gene_config$target_exons
    window_gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(plot_start, plot_end))
    ov <- GenomicRanges::findOverlaps(exons_full, window_gr, type = "any")
    if (length(ov) > 0) {
      exons_plot <- exons_full[unique(S4Vectors::queryHits(ov))]
      exons_plot <- exons_plot[order(start(exons_plot))]
      
      if (is.null(exons_plot$exon_idx) && !is.null(exons_plot$exon_number)) exons_plot$exon_idx <- exons_plot$exon_number
      else if (is.null(exons_plot$exon_idx)) exons_plot$exon_idx <- seq_along(exons_plot)
      
      exon_gr <- GenomicRanges::GRanges(
        seqnames = chrom, ranges = IRanges::IRanges(start = start(exons_plot), end = end(exons_plot)), strand = strand(exons_plot)
      )
      exon_gr$gene       <- gene_config$gene
      exon_gr$transcript <- paste0(gene_config$gene, "_tx")
      exon_gr$symbol     <- gene_config$gene
      exon_gr$exon       <- paste0("E", exons_plot$exon_idx)
      
      grtrack <- Gviz::GeneRegionTrack(
        exon_gr, genome = genome_build, chromosome = chrom,
        name = paste0(gene_config$gene, " targeted exons"),
        fill = "#2ecc71", col = "#1a8a4a",
        transcriptAnnotation = "none", exonAnnotation = "exon",
        showExonId = TRUE, fontsize = 11, cex.feature = 1.0,
        fontcolor.feature = "black", shape = "box",
        background.title = "transparent", col.border.title = NA,
        fontcolor.title = "black"
      )
    }
  }

  itd_labels <- paste0("  ", sprintf("%d bp | VAF %s%% | %s", dup_lens, vafs, orient_display))
  itd_track  <- Gviz::AnnotationTrack(
    start = breakpoints, end = dup_ends, chromosome = chrom,
    genome = genome_build, name = ifelse(length(breakpoints) > 1, "ITDs", "ITD"),
    id = itd_labels, group = itd_labels, groupAnnotation = "group", just.group = "right",
    showFeatureId = FALSE, fill = "#f39c12", col = "#c0392b",
    fontcolor.group = "#4a1c00", cex.group = 0.8, cex.title = 0.9,
    rotation.title = 0, background.title = "transparent", col.border.title = NA, fontcolor.title = "black"
  )

  axis_track <- Gviz::GenomeAxisTrack(cex = 0.7, labelPos = "below", add35 = FALSE, add53 = FALSE, col = "black")
  track_list <- if (!is.null(grtrack)) list(axis_track, grtrack, cov_track, itd_track) else list(axis_track, cov_track, itd_track)
  t_sizes <- if (!is.null(grtrack)) c(1.0, 0.6, 8, max(0.5, 0.2 * length(breakpoints))) else c(1.0, 8, max(0.5, 0.2 * length(breakpoints)))

  grDevices::pdf(output_pdf, width = width, height = height, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  main_title <- sprintf("Sample: %s  |  %s  |  %d ITD(s)  |  %s", sample_name, gene_config$gene, length(breakpoints), genome_build)
  subtitle   <- sprintf("Positions: %s  •  Lengths: %s  •  VAFs: %s%%", paste(breakpoints, collapse = ", "), paste(dup_lens, collapse = ", "), paste(vafs, collapse = ", "))

  Gviz::plotTracks(
    track_list, from = plot_start - 50L, to = plot_end + 50L,
    sizes = t_sizes, collapse = FALSE, main = main_title, cex.main = 1.2,
    fontcolor.title = "#2c3e50", col.axis = "black", cex.axis = 0.8, margins = c(5.0, 5, 5.5, 5)
  )
  grid::grid.text(label = subtitle, x = grid::unit(0.5, "npc"), y = grid::unit(0.915, "npc"), just = c("center", "center"), gp = grid::gpar(fontsize = 9, col = "#2c3e50", fontface = "plain"))

  summary_df <- data.frame(
    Metric = c("Sample", "Gene", "Genome Build", "Genomic Position", "Length (bp)", "VAF (%)", "Supporting Reads", "Wildtype Reads", "Depth at Breakpoint", "HGVS cDNA", "HGVS Protein", "Region", "Exon Number", "Orientation", "Hotspot", "Strand Bias", "Mean MAPQ", "Support Consistency", "RefMatch (%)", "ITD Read Coverage (%)", "Sequence Source"), stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(itd_df[keep, , drop = FALSE]))) {
    row <- itd_df[keep, ][i, ]
    summary_df[[paste0("ITD_", i)]] <- c(
      as.character(sample_name), as.character(row$Gene), genome_build, as.character(row$GenomicPosition), as.character(row$Length), sprintf("%.1f", row$AlleleFrequency * 100), as.character(row$SupportingReads), as.character(round(row$WildtypeReads)), as.character(row$DepthAtBreakpoint), ifelse(is.na(row$HGVS_cDNA), "N/A", row$HGVS_cDNA), ifelse(is.na(row$HGVS_Protein), "N/A", row$HGVS_Protein), ifelse(is.na(row$Region), "N/A", row$Region), ifelse(is.na(row$ExonNumber), "N/A", as.character(row$ExonNumber)), ifelse(is.na(row$Orientation), "N/A", row$Orientation), ifelse(isTRUE(row$Hotspot), as.character(row$HotspotName), "No"), ifelse(is.na(row$StrandBias), "N/A", sprintf("%.3f", row$StrandBias)), ifelse(is.na(row$MeanSupportMAPQ), "N/A", sprintf("%.1f", row$MeanSupportMAPQ)), ifelse(is.na(row$SupportConsistency), "N/A", paste0(round(row$SupportConsistency), "%")), ifelse(is.na(row$RefMatch_Observed), "N/A", paste0(round(row$RefMatch_Observed, 1), "%")), ifelse(is.na(row$ITDReadCoverage), "N/A", paste0(round(row$ITDReadCoverage, 1), "%")), ifelse(is.na(row$SequenceSource), "N/A", row$SequenceSource)
    )
  }

  grid::grid.newpage()
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    n_rows   <- nrow(summary_df)
    row_fill <- rep(c("white", "#f0f4f8"), length.out = n_rows)
    tbl <- gridExtra::tableGrob(summary_df, rows = NULL, theme = gridExtra::ttheme_minimal(core = list(fg_params = list(cex = 0.8, fontface = "plain"), bg_params = list(fill = row_fill, col = NA)), colhead = list(fg_params = list(cex = 0.9, fontface = "bold", col = "white"), bg_params = list(fill = "#2c3e50", col = NA))))
    grid::grid.text(sprintf("Variant Summary – Sample: %s | %s (%d variants)", sample_name, gene_config$gene, nrow(itd_df[keep, , drop = FALSE])), x = 0.5, y = 0.97, gp = grid::gpar(fontsize = 13, fontface = "bold", col = "#2c3e50"))
    vp <- grid::viewport(x = 0.5, y = 0.47, width = 0.90, height = 0.88)
    grid::pushViewport(vp); grid::grid.draw(tbl); grid::popViewport()
  } else {
    grid::grid.text(paste(capture.output(print(summary_df)), collapse = "\n"), x = 0.05, y = 0.95, just = "left", gp = grid::gpar(cex = 0.7, fontfamily = "mono"))
  }

  if (show_config) {
    grid::grid.newpage()
    grid::grid.text("Analysis Configuration", x = 0.5, y = 0.95, gp = grid::gpar(fontsize = 14, fontface = "bold", col = "#2c3e50"))
    ypos <- 0.85
    config_lines <- c(
      paste("Gene:", gene_config$gene), paste("Build:", genome_build), paste("Chromosome:", gene_config$chrom), paste("Strand:", gene_config$strand), paste("Transcript reference:", ifelse(is.na(gene_config$transcript), "none", gene_config$transcript)), paste("CDS offset:", gene_config$cds_offset), paste("Genomic window:", gene_config$genomic_start, "-", gene_config$genomic_end), paste("Targeted exons:", paste(unique(gene_config$target_exons$exon_idx), collapse = ", ")), paste("Settings:"), paste("  min_support =", gene_config$gene_settings$min_support %||% "default"), paste("  cluster_tolerance =", gene_config$gene_settings$cluster_tolerance %||% "default"), paste("  min_mapq =", gene_config$gene_settings$min_mapq %||% "default")
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

  plot_start <- gene_config$genomic_start
  plot_end   <- gene_config$genomic_end

  if (is.null(precomputed_cov)) {
    param <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges(chrom, IRanges::IRanges(plot_start, plot_end)))
    gal   <- GenomicAlignments::readGAlignments(bam_path, param = param)
    cov   <- GenomicAlignments::coverage(gal)[[chrom]]
  } else {
    cov <- precomputed_cov
  }
  pos   <- seq(plot_start, plot_end)
  depth <- as.numeric(cov[pos])

  p_cov <- plotly::plot_ly(
    x = pos, y = depth, type = "bar", name = "Coverage",
    marker = list(color = "#4a90d9", line = list(color = "#4a90d9", width = 0)),
    hovertemplate = "Pos: %{x}<br>Depth: %{y}<extra></extra>"
  )
  p_cov <- plotly::layout(p_cov, yaxis = list(title = "Depth", gridcolor = "#e8e8e8"), xaxis = list(title = "", showticklabels = FALSE), bargap = 0, plot_bgcolor = "white", paper_bgcolor = "white")

  p_exon <- plotly::plot_ly()
  p_exon <- plotly::layout(p_exon, xaxis = list(title = "", showticklabels = FALSE, range = c(plot_start, plot_end), rangeslider = list(visible = FALSE)), yaxis = list(title = "Exons", showticklabels = FALSE, range = c(0, 1), fixedrange = TRUE), plot_bgcolor = "#f8fafc", paper_bgcolor = "white")

  if (!is.null(gene_config$target_exons) && length(gene_config$target_exons) > 0) {
    exons_full <- gene_config$target_exons
    window_gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(plot_start, plot_end))
    ov <- GenomicRanges::findOverlaps(exons_full, window_gr, type = "any")
    if (length(ov) > 0) {
      exons_plot <- exons_full[unique(S4Vectors::queryHits(ov))]
      exons_plot <- exons_plot[order(start(exons_plot))]
      exon_idx_vec <- if (!is.null(exons_plot$exon_idx)) exons_plot$exon_idx else if (!is.null(exons_plot$exon_number)) exons_plot$exon_number else seq_along(exons_plot)
      
      for (i in seq_along(exons_plot)) {
        es <- start(exons_plot[i]); ee <- end(exons_plot[i])
        p_exon <- p_exon |> plotly::add_trace(x = c(es, ee, ee, es, es), y = c(0.1, 0.1, 0.9, 0.9, 0.1), type = "scatter", mode = "lines", fill = "toself", fillcolor = "rgba(46,204,113,0.7)", line = list(color = "#1a8a4a", width = 1), name = paste0("Exon ", exon_idx_vec[i]), hovertemplate = sprintf("<b>Exon %s</b><br>%d–%d<extra></extra>", exon_idx_vec[i], es, ee), showlegend = FALSE) |> plotly::add_annotations(x = (es + ee) / 2, y = 0.5, xref = "x", yref = "y", text = paste0("E", exon_idx_vec[i]), showarrow = FALSE, font = list(size = 10, color = "#1a5e30"))
      }
    }
  }

  pal <- c("#e67e22", "#e74c3c", "#9b59b6", "#1abc9c", "#3498db", "#f39c12", "#c0392b", "#8e44ad", "#16a085", "#2980b9")
  p_itd <- plotly::plot_ly()

  for (i in seq_len(nrow(itd_df))) {
    row      <- itd_df[i, ]
    vaf_str  <- sprintf("%.1f%%", row$AlleleFrequency * 100)
    col_i    <- pal[((i - 1L) %% length(pal)) + 1L]
    y_pos    <- 0.1 + (i - 1L) * 0.18
    tip_text <- paste0("<b>ITD ", i, "</b><br>", "Breakpoint: ", row$GenomicPosition, "<br>", "Length: ", row$Length, " bp<br>", "VAF: ", vaf_str, "<br>", "HGVS cDNA: ", ifelse(is.na(row$HGVS_cDNA), "N/A", row$HGVS_cDNA), "<br>", "HGVS Protein: ", ifelse(is.na(row$HGVS_Protein), "N/A", row$HGVS_Protein), "<br>", "Region: ", ifelse(is.na(row$Region), "N/A", row$Region), "<br>", "Support: ", row$SupportingReads, " reads<br>", "RefMatch: ", ifelse(is.na(row$RefMatch_Observed), "N/A", paste0(round(row$RefMatch_Observed, 1), "%")), "<br>", "ITD Cov: ", ifelse(is.na(row$ITDReadCoverage), "N/A", paste0(round(row$ITDReadCoverage, 1), "%")), "<br>", "Orientation: ", ifelse(is.na(row$Orientation), "?", row$Orientation))
    label <- sprintf("%d bp  VAF %s", row$Length, vaf_str)
    if (!is.na(row$HGVS_cDNA)) label <- paste0(label, "  ", row$HGVS_cDNA)

    p_itd <- p_itd |> plotly::add_trace(x = c(row$GenomicPosition, dup_ends[i]), y = c(y_pos, y_pos), type = "scatter", mode = "lines", line = list(color = col_i, width = 14), name = label, text = tip_text, hoverinfo = "text", showlegend = TRUE) |> plotly::add_trace(x = row$GenomicPosition, y = y_pos, type = "scatter", mode = "markers", marker = list(symbol = "line-ns", size = 14, color = col_i, line = list(color = "white", width = 2)), hoverinfo = "skip", showlegend = FALSE)
  }

  p_itd <- plotly::layout(p_itd, xaxis = list(title = paste0(chrom, " (", genome_build, ")"), gridcolor = "#e8e8e8", rangeslider = list(visible = FALSE)), yaxis = list(title = "ITDs", showticklabels = FALSE, range = c(0, max(1, nrow(itd_df) * 0.18 + 0.2)), fixedrange = TRUE), legend = list(orientation = "h", x = 0, y = -0.25, font = list(size = 11)), plot_bgcolor = "white", paper_bgcolor = "white")

  title_text <- if (!is.null(sample_name)) paste0("TALOS • ", sample_name, "  |  ", gene_config$gene, "  |  ", nrow(itd_df), " ITD(s)  |  ", genome_build) else paste0("TALOS • ", gene_config$gene, "  |  ", nrow(itd_df), " ITD(s)  |  ", genome_build)

  fig <- plotly::subplot(p_cov, p_exon, p_itd, nrows = 3, shareX = TRUE, heights = c(0.58, 0.14, 0.28))
  fig <- plotly::layout(fig, title = list(text = title_text, x = 0.5, font = list(size = 14)), hovermode = "x unified", showlegend = TRUE, margin = list(t = 60))
  
  if (show_config) {
    config_text <- paste("Config:", paste("Targeted exons:", paste(unique(gene_config$target_exons$exon_idx), collapse = ",")), paste("Transcript:", ifelse(is.na(gene_config$transcript), "none", gene_config$transcript)), paste("CDS offset:", gene_config$cds_offset), sep = " | ")
    fig <- plotly::layout(fig, annotations = list(text = config_text, x = 0, y = 1.05, xref = "paper", yref = "paper", showarrow = FALSE, font = list(size = 9)))
  }

  fig
}


#' Generate an HTML report for TALOS results
#' @param result_df Data frame from \code{\link{talos}}.
#' @param bam_paths Named character vector.
#' @param gene_configs Named list of gene configs.
#' @param output_file Path to the output HTML file.
#' @param title Report title.
#' @param show_config Include configuration details.
#' @param ... Unused compatibility parameters.
#' @export
talos_html_report <- function(result_df, bam_paths = NULL, gene_configs = NULL, output_file = "TALOS_report.html", title = "TALOS Analysis Report", show_config = FALSE, ...) {
  needs <- c("rmarkdown", "DT")
  missing_pkg <- needs[!vapply(needs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing_pkg) > 0) stop("Please install: ", paste(missing_pkg, collapse = ", "))

  if (!grepl("^[A-Za-z]:|^/", output_file)) output_file <- file.path(getwd(), output_file)
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  plot_widgets <- list()
  has_plotly   <- requireNamespace("plotly", quietly = TRUE)
  if (has_plotly && !is.null(bam_paths) && !is.null(gene_configs)) {
    for (sn in names(bam_paths)) {
      bam <- bam_paths[[sn]]
      if (!file.exists(bam)) next
      for (gn in names(gene_configs)) {
        sub_df <- result_df[result_df$Sample == sn & result_df$Gene == gn, , drop = FALSE]
        if (nrow(sub_df) == 0) next
        key <- paste(sn, gn, sep = "_")
        plot_widgets[[key]] <- tryCatch(plot_talos_interactive(sub_df, gene_configs[[gn]], bam, sn, show_config = show_config), error = function(e) { message(sprintf("[HTML report] Could not generate plot for %s/%s: %s", sn, gn, conditionMessage(e))); NULL })
      }
    }
  }

  rmd_template <- tempfile(fileext = ".Rmd")
  on.exit(unlink(rmd_template, force = TRUE), add = TRUE)

  rmd_lines <- c(
    "---", paste0('title: "', gsub('"', "'", title), '"'), "output:", "  html_document:", "    self_contained: true", "    theme: flatly", "    toc: true", "    toc_float: true", "params:", "  result_df: NULL", "  plot_widgets: NULL", "  gene_configs: NULL", "  has_plotly: false", "  show_config: false", "---", "", "```{r setup, include=FALSE}", "knitr::opts_chunk$set(echo=FALSE, warning=FALSE, message=FALSE)", "library(DT)", "```", "", "## Results summary", "```{r results-table}", "cols_hide <- c('ITD_Sequence','ITDCoverageRLE','SequenceImputed','SequencePartial')", "cols_show <- setdiff(names(params$result_df), cols_hide)", "df_disp   <- params$result_df[, cols_show, drop=FALSE]", "if ('AlleleFrequency' %in% colnames(df_disp)) df_disp$AlleleFrequency <- sprintf('%.1f%%', df_disp$AlleleFrequency * 100)", "DT::datatable(df_disp, filter='top', options=list(scrollX=TRUE, pageLength=25, dom='Bfrtip'), extensions='Buttons')", "```", "", "## Variant details", "```{r variant-details, results='asis'}", "if (nrow(params$result_df) == 0) { cat('*No variants detected.*\\n') } else { groups <- unique(params$result_df[, c('Sample','Gene')]); for (gi in seq_len(nrow(groups))) { sn <- groups$Sample[gi]; gn <- groups$Gene[gi]; cat('\\n\\n### Sample:', sn, '\\u2014', gn, '\\n\\n'); if (isTRUE(params$has_plotly)) { key <- paste(sn, gn, sep='_'); wgt <- params$plot_widgets[[key]]; if (!is.null(wgt)) { print(htmltools::tagList(wgt)) } else { cat('*Interactive plot not available.*\\n') } }; sub_df <- params$result_df[params$result_df$Sample==sn & params$result_df$Gene==gn, , drop=FALSE]; for (vi in seq_len(nrow(sub_df))) { row <- sub_df[vi, ]; hgvs_str <- ifelse(is.na(row$HGVS_cDNA), 'N/A', row$HGVS_cDNA); cat(sprintf('\\n\\n#### Variant %d: %s\\n', vi, hgvs_str)); cat(sprintf('* **VAF:** %.1f%%\\n', row$AlleleFrequency*100)); cat(sprintf('* **Length:** %d bp\\n', row$Length)); cat(sprintf('* **Breakpoint:** %s\\n', row$GenomicPosition)); cat(sprintf('* **Region:** %s\\n', ifelse(is.na(row$Region),'N/A',row$Region))); cat(sprintf('* **Orientation:** %s\\n', ifelse(is.na(row$Orientation),'?',row$Orientation))); cat(sprintf('* **Hotspot:** %s\\n', ifelse(isTRUE(row$Hotspot),row$HotspotName,'No'))); cat(sprintf('* **RefMatch:** %s\\n', ifelse(is.na(row$RefMatch_Observed),'N/A',paste0(round(row$RefMatch_Observed,1),'%')))); cat(sprintf('* **ITD Read Coverage:** %s\\n', ifelse(is.na(row$ITDReadCoverage),'N/A',paste0(round(row$ITDReadCoverage,1),'%')))); cat(sprintf('* **Sequence Source:** %s\\n', ifelse(is.na(row$SequenceSource),'N/A',row$SequenceSource))); if (!is.na(row$ITD_Sequence) && nchar(row$ITD_Sequence) > 0) { cat(sprintf('* **ITD Sequence:** <details><summary>Show (%d bp)</summary><pre>%s</pre></details>\\n', nchar(row$ITD_Sequence), row$ITD_Sequence)) } }; if (isTRUE(params$show_config)) { cat('\\n\\n#### Analysis configuration\\n'); cfg_str <- sprintf(paste0('* **Transcript:** %s\\n* **CDS offset:** %d\\n* **Targeted exons:** %s\\n* **Genomic window:** %d-%d\\n'), ifelse(is.na(params$gene_configs[[gn]]$transcript), 'none', params$gene_configs[[gn]]$transcript), params$gene_configs[[gn]]$cds_offset, paste(unique(params$gene_configs[[gn]]$target_exons$exon_idx), collapse=','), params$gene_configs[[gn]]$genomic_start, params$gene_configs[[gn]]$genomic_end); cat(cfg_str) } } }", "```"
  )
  writeLines(rmd_lines, rmd_template)

  rmarkdown::render(rmd_template, output_file = output_file, params = list(result_df = result_df, plot_widgets = plot_widgets, gene_configs = gene_configs, has_plotly = has_plotly, show_config = show_config), quiet = TRUE, ...)
  message("HTML report written to: ", normalizePath(output_file))
  invisible(output_file)
}

.write_talos_output <- function(final_df, output_prefix, output_folder, sample_name, gene_name, gene_config, ref_dna, bam_path, write_vcf, plot, html_report, verbose, add_config_to_report = FALSE) {
  if (is.null(output_prefix) || nrow(final_df) == 0) return(invisible(NULL))
  sample_folder <- file.path(output_folder, sample_name)
  if (!dir.exists(sample_folder)) dir.create(sample_folder, recursive = TRUE, showWarnings = FALSE)

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  base_name <- paste0(output_prefix, "_", gene_name, "_", timestamp)
  tsv_file  <- file.path(sample_folder, paste0(base_name, ".tsv"))

  for (col in names(final_df)) {
    if (is.list(final_df[[col]])) final_df[[col]] <- vapply(final_df[[col]], function(x) if (length(x) == 0) NA_character_ else paste(x, collapse = ";"), character(1L))
  }
  write.table(final_df, file = tsv_file, sep = "\t", quote = FALSE, row.names = FALSE, na = ".")
  if (verbose) message(sprintf("Results written to: %s", tsv_file))

  if (write_vcf) {
    vcf_file <- file.path(sample_folder, paste0(base_name, ".vcf"))
    write_itd_vcf(itd_df = final_df, ref_dna = ref_dna, genomic_start = gene_config$genomic_start, chrom = gene_config$chrom, sample_name = sample_name, genome_build = gene_config$build, vcf_path = vcf_file, overwrite = TRUE)
  }

  pdf_path <- NULL
  if (plot && nrow(final_df) > 0) {
    pdf_name <- file.path(sample_folder, paste0(sample_name, "_", gene_name, ".pdf"))
    plot_talos_report(itd_row = final_df, gene_config = gene_config, bam_path = bam_path, sample_name = sample_name, output_pdf = pdf_name, show_config = add_config_to_report)
    pdf_path <- pdf_name
  }

  widget_path <- NULL
  if (requireNamespace("plotly", quietly = TRUE)) {
    wgt <- tryCatch(plot_talos_interactive(final_df, gene_config, bam_path, sample_name, show_config = add_config_to_report), error = function(e) NULL)
    if (!is.null(wgt)) {
      widget_file <- file.path(sample_folder, paste0(sample_name, "_", gene_name, "_interactive.html"))
      if (requireNamespace("htmlwidgets", quietly = TRUE)) {
        htmlwidgets::saveWidget(wgt, widget_file, selfcontained = TRUE, title = paste0("TALOS – ", sample_name, " | ", gene_name))
        if (verbose) message(sprintf("Interactive plot saved to: %s", widget_file))
        widget_path <- widget_file
      }
    }
  }

  if (html_report && nrow(final_df) > 0) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) { warning("rmarkdown not installed. Cannot generate HTML report.") } 
    else {
      report_file <- file.path(sample_folder, paste0("TALOS_report_", gene_name, "_", timestamp, ".html"))
      bam_vec     <- setNames(bam_path, sample_name)
      config_list <- setNames(list(gene_config), gene_name)
      suppressWarnings(talos_html_report(result_df = final_df, bam_paths = bam_vec, gene_configs = config_list, output_file = report_file, title = paste("TALOS Report –", gene_name, sample_name), show_config = add_config_to_report))
      if (verbose) message(sprintf("HTML report written to: %s", report_file))
    }
  }
  invisible(list(tsv = tsv_file, widget = widget_path, pdf = pdf_path))
}

.log_duration <- function(gene_name, sample_name, start_time, verbose = TRUE) {
  if (!verbose) return(invisible(NULL))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  msg     <- if (elapsed < 60) sprintf("%.1f sec", elapsed) else sprintf("%.1f min", elapsed / 60)
  message(sprintf("[TALOS] Analysis complete | Gene: %s | Sample: %s | Duration: %s", gene_name, sample_name, msg))
}
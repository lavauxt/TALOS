# ============================================================================
# Main detection function (Memory Optimized)
# ============================================================================

#' Detect ITDs/PTDs using the TALOS Algorithm
#'
#' @param bam_path Path to an indexed BAM file.
#' @param gene_config Gene configuration list from \code{\link{get_gene_config}}.
#' @param genomic_ref_seq Optional override for the genomic reference sequence.
#' @param k K-mer size (default 11).
#' @param min_support Minimum support reads (default 10).
#' @param min_size Minimum duplication length (default 10).
#' @param max_missing_kmers Proportion threshold to drop reads.
#' @param cluster_tolerance Radius (bp) to cluster breakpoints.
#' @param prefilter Enable fast prefiltering of alignments.
#' @param min_ins_filter Minimum net insertion length filter.
#' @param output_prefix Base name prefix (default "TALOS").
#' @param output_folder Output folder path.
#' @param write_vcf Generate VCF report.
#' @param verbose Logging enabled.
#' @param debug Debug info.
#' @param nominal_read_len Assumed length.
#' @param max_correction Size-bias correction cap.
#' @param plot Output PDF.
#' @param sample_name Extracted/assigned sample name.
#' @param filter_intronic Auto-drop non-exonic breakpoints.
#' @param ptd_mode Run PTD softclip-based parsing.
#' @param use_cigar_bp Rely on CIGAR instead of K-mers for BP refine.
#' @param refine_bp Enable secondary refinement loop.
#' @param min_mapq Filter param.
#' @param min_wt_reads Require baseline WT.
#' @param vaf_threshold Minimum relative allele frequency.
#' @param min_strand_bias Lower bound SB.
#' @param max_strand_bias Upper bound SB.
#' @param min_mean_support_mapq Strict metric cutoff.
#' @param max_breakpoint_spread Strict metric cutoff.
#' @param min_softclip_fraction Strict metric cutoff.
#' @param min_unique_breakpoints Strict metric cutoff.
#' @param min_coverage_drop Strict metric cutoff (default 1.5, adjust based on expected coverage).
#' @param coverage_window Internal calc usage.
#' @param min_microhomology Required MH length.
#' @param min_discordant_ratio Required ratio length (for filtering).
#' @param min_entropy Sequence repeat complexity baseline.
#' @param detect_orientation Check inverse loops.
#' @param annotate_hotspots SQLite DB lookups.
#' @param hotspot_db_path Path to SQLite lookup DB.
#' @param compute_alignment_score Metric toggle.
#' @param compute_support_bases Metric toggle.
#' @param compute_consistency Metric toggle.
#' @param compute_itd_coverage Metric toggle.
#' @param compute_coverage_drop Metric toggle.
#' @param compute_microhomology Metric toggle.
#' @param compute_repeat_entropy Metric toggle.
#' @param compute_discordant_ratio Metric toggle (boolean, default TRUE).
#' @param compute_hgvs Produce HGVS format (only works if transcript is assigned).
#' @param html_report HTML overview toggle.
#' @param use_gviz Use graphical views.
#' @param min_itd_read_coverage Post-calc check ratio.
#' @param max_reads_in_region Loading safeguard.
#' @param add_config_to_report Render configs in results.
#' @param max_pairwise_alignments Max limit for computationally heavy sub-alignments.
#' @param include_sample Include sample name in output filenames (default TRUE).
#' @param include_gene Include gene name in output filenames (default TRUE).
#' @param include_timestamp Include timestamp in output filenames (default TRUE).
#' @param output_sep Separator between filename parts (default "_").
#' @param global_log Write all logs (config, compute time, errors) to a single file (default TRUE).
#' @param min_side_softclip_reads Minimum number of reads with soft‑clip on each side (default 10).
#' @param max_side_ratio Maximum ratio of left/right soft‑clip counts before filtering (default 10).
#' @param min_softclip_pct_side Minimum percentage of supporting reads that must have a soft‑clip on each side (default 1.0).
#' @param ... Other parameters passed to the underlying engine.
#' @return Data frame of ITD calls with diagnostic columns, invisibly.
#' @export
detect_itd <- function(
    bam_path, gene_config, genomic_ref_seq = NULL, k = 11, min_support = 10,
    min_size = 15, max_missing_kmers = 0.5, cluster_tolerance = 10, prefilter = TRUE,
    min_ins_filter = 3, output_prefix = "TALOS", output_folder = "./results",
    write_vcf = TRUE, verbose = TRUE, debug = FALSE, nominal_read_len = 150,
    max_correction = 2.0, plot = FALSE, sample_name = NULL, filter_intronic = FALSE,
    ptd_mode = FALSE, use_cigar_bp = TRUE, refine_bp = FALSE, min_mapq = 20,
    min_wt_reads = 0, vaf_threshold = 0.01, min_strand_bias = 0, max_strand_bias = 1,
    min_mean_support_mapq = 0, max_breakpoint_spread = Inf, min_softclip_fraction = 0,
    min_unique_breakpoints = 0, min_coverage_drop = 1.0, coverage_window = 200,
    min_microhomology = 0, min_discordant_ratio = 0, min_entropy = 0,
    detect_orientation = TRUE, do_annotate_hotspots = TRUE, hotspot_db_path = NULL,
    compute_alignment_score = TRUE, compute_support_bases = TRUE,
    compute_consistency = TRUE, compute_itd_coverage = TRUE,
    compute_coverage_drop = TRUE, compute_microhomology = TRUE,
    compute_repeat_entropy = TRUE, compute_discordant_ratio = TRUE,
    compute_hgvs = TRUE, html_report = TRUE, use_gviz = TRUE,
    min_itd_read_coverage = 0, max_reads_in_region = 200000,
    add_config_to_report = FALSE, max_pairwise_alignments = 30L,
    include_sample = TRUE, include_gene = TRUE, include_timestamp = TRUE,
    output_sep = "_", global_log = TRUE,
    min_side_softclip_reads = 10, max_side_ratio = 10,
    min_softclip_pct_side = 1.0,
    min_left_softclip_pct_wt = 1.0,
    min_right_softclip_pct_wt = 1.0,
    ...
) {
  
  start_time <- Sys.time()
  gene_name  <- if (!is.null(gene_config$gene)) gene_config$gene else "UNKNOWN"
  if (is.null(sample_name)) sample_name <- sub("\\..*", "", basename(bam_path))
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  sample_folder <- file.path(output_folder, sample_name)
  if (!dir.exists(sample_folder)) dir.create(sample_folder, recursive = TRUE, showWarnings = FALSE)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # ---- Build unified base name ----
  base_parts <- c()
  if (!is.null(output_prefix) && output_prefix != "") base_parts <- c(base_parts, output_prefix)
  if (include_sample && !is.null(sample_name) && sample_name != "") base_parts <- c(base_parts, sample_name)
  if (include_gene   && !is.null(gene_name)   && gene_name != "")   base_parts <- c(base_parts, gene_name)
  if (include_timestamp && !is.null(timestamp) && timestamp != "")  base_parts <- c(base_parts, timestamp)
  base_name <- paste(base_parts, collapse = output_sep)
  
  # ---- Global log file (single file for everything) ----
  global_log_file <- NULL
  if (global_log) {
    global_log_file <- file.path(sample_folder, paste0(base_name, ".log"))
    cat("=== TALOS Analysis Log ===\n", file = global_log_file)
    cat(paste("Sample:", sample_name, "\n"), file = global_log_file, append = TRUE)
    cat(paste("Gene:", gene_name, "\n"), file = global_log_file, append = TRUE)
    cat(paste("BAM:", bam_path, "\n"), file = global_log_file, append = TRUE)
    cat(paste("Timestamp:", timestamp, "\n"), file = global_log_file, append = TRUE)
    cat(paste("Output folder:", output_folder, "\n"), file = global_log_file, append = TRUE)
    cat(strrep("-", 40), "\n", file = global_log_file, append = TRUE)
    cat("=== Gene Configuration ===\n", file = global_log_file, append = TRUE)
    capture.output(str(gene_config, max.level = 2), file = global_log_file, append = TRUE)
    cat(strrep("-", 40), "\n", file = global_log_file, append = TRUE)
  }
  
  # Helper to write to both console and log
  log_msg <- function(..., append = TRUE) {
    msg <- paste(..., collapse = " ")
    if (verbose) message(msg)
    if (global_log && !is.null(global_log_file)) {
      cat(msg, "\n", file = global_log_file, append = append)
    }
  }
  
  log_msg("Starting TALOS detection engine...")
  
  on.exit({
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    time_msg <- if (elapsed < 60) sprintf("%.1f sec", elapsed) else sprintf("%.1f min", elapsed / 60)
    log_msg(sprintf("\n%s\nCOMPUTE TIME: %s\n%s", strrep("-", 40), time_msg, strrep("-", 40)))
    if (global_log && !is.null(global_log_file)) {
      log_msg(sprintf("Global log saved to: %s", global_log_file))
    }
  }, add = TRUE)
  
  tryCatch({
    # Use genomic_ref_seq from gene_config (always present)
    ref_dna <- gene_config$genomic_ref_seq
    if (is.null(ref_dna) || is.na(ref_dna) || nchar(ref_dna) == 0)
      stop("No genomic reference sequence available in gene_config.")
    
    ref_len <- nchar(ref_dna)
    if (ref_len < k) {
      log_msg("Reference shorter than k-mer length. Skipping.")
      return(data.frame())
    }
    
    ref_kmers <- .prepare_kmers(ref_dna, k)
    bam_data  <- .load_bam_data_streaming(bam_path, gene_config, 
                                          compute_pairs = compute_discordant_ratio,
                                          max_reads = max_reads_in_region)
    all_reads <- bam_data$reads
    
    if (length(all_reads) == 0) {
      log_msg("No reads found in target region.")
      .log_duration(gene_name, sample_name, start_time, verbose)
      return(data.frame())
    }
    
    wt_mapqs <- S4Vectors::mcols(all_reads)$mapq
    wt_mapqs[is.na(wt_mapqs)] <- 0L
    all_reads_for_wt <- all_reads[wt_mapqs >= min_mapq]
    
    if (prefilter) {
      all_reads <- .filter_reads_by_cigar(all_reads, min_mapq, min_ins_filter)
      if (length(all_reads) == 0) {
        log_msg("No reads passed CIGAR prefilter.")
        return(data.frame())
      }
    }
    
    if (ptd_mode && use_cigar_bp) {
      candidates <- .extract_candidates_ptd(all_reads, gene_config$genomic_start, ref_len)
    } else {
      candidates <- .extract_candidates_standard(all_reads, ref_kmers, ptd_mode, min_size, max_missing_kmers, refine_bp, use_cigar_bp, gene_config$genomic_start, ref_dna)
    }
    
    if (length(candidates) == 0) {
      log_msg("No candidate reads (clips or k-mer jumps) found.")
      return(data.frame())
    }
    
    bp_df    <- .candidates_to_df(candidates, gene_config$genomic_start)
    clusters <- .cluster_breakpoints(bp_df$breakpoint, cluster_tolerance)
    
    wt_info              <- .prepare_wildtype_info(all_reads_for_wt, gene_config$genomic_start, gene_config$genomic_end)
    wt_info$sample_name  <- sample_name
    
    cand_breakpoints <- bp_df$breakpoint
    cand_lengths     <- bp_df$length
    
    thresholds <- list(
        min_coverage_drop = min_coverage_drop, min_microhomology = min_microhomology,
        min_discordant_ratio = min_discordant_ratio, min_entropy = min_entropy,
        min_strand_bias = min_strand_bias, max_strand_bias = max_strand_bias,
        min_mean_support_mapq = min_mean_support_mapq, max_breakpoint_spread = max_breakpoint_spread,
        min_softclip_fraction = min_softclip_fraction, min_unique_breakpoints = min_unique_breakpoints,
        min_itd_read_coverage = min_itd_read_coverage,
        min_side_softclip_reads = min_side_softclip_reads,
        max_side_ratio          = max_side_ratio,
        min_softclip_pct_side   = min_softclip_pct_side,
        min_left_softclip_pct_wt = min_left_softclip_pct_wt,
        min_right_softclip_pct_wt = min_right_softclip_pct_wt
      )
    
    results <- list()
    for (cl in clusters) {
      cluster_candidates <- which(cand_breakpoints %in% cl)
      if (length(cluster_candidates) == 0) next
      
      cluster_lengths <- cand_lengths[cluster_candidates]
      cluster_lengths[is.na(cluster_lengths)] <- -1L
      len_groups <- split(cluster_candidates, cluster_lengths)
      
      for (len in names(len_groups)) {
        best_len <- if (as.integer(len) == -1L) NA_integer_ else as.integer(len)
        idxs <- len_groups[[len]]
        
        if (is.na(best_len)) support_rows <- bp_df[bp_df$breakpoint %in% cl & is.na(bp_df$length), ]
        else support_rows <- bp_df[bp_df$breakpoint %in% cl & !is.na(bp_df$length) & bp_df$length == best_len, ]
        
        metrics <- .compute_variant_metrics(
          cluster_bps = cl, best_len = best_len, support_rows = support_rows,
          genomic_start = gene_config$genomic_start, ref_dna = ref_dna,
          gene_config = gene_config, all_reads_cov = bam_data$cov,
          all_pairs = bam_data$pairs, wildtype_info = wt_info, ptd_mode = ptd_mode,
          min_support = min_support, min_wt_reads = min_wt_reads,
          nominal_read_len = nominal_read_len, max_correction = max_correction,
          vaf_threshold = vaf_threshold, do_alignment_score = compute_alignment_score,
          do_support_bases = compute_support_bases, do_consistency = compute_consistency,
          do_itd_coverage = compute_itd_coverage, do_coverage_drop = compute_coverage_drop,
          do_microhomology = compute_microhomology, do_repeat_entropy = compute_repeat_entropy,
          do_discordant_ratio = compute_discordant_ratio, 
          do_detect_orientation = detect_orientation,
          do_hgvs = compute_hgvs, max_pairwise_alignments = max_pairwise_alignments,
          debug = debug
        )
        if (!is.null(metrics) && .apply_filters(metrics, thresholds)) results[[length(results) + 1L]] <- metrics
      }
    }
    
    if (length(results) == 0) {
      log_msg("No variants passed all filters.")
      return(data.frame())
    }
    
    final_df <- do.call(rbind, lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
    
    if (do_annotate_hotspots) final_df <- annotate_hotspots(final_df, db_path = hotspot_db_path, genome_build = gene_config$build)
    else { final_df$Hotspot <- FALSE; final_df$HotspotName <- NA_character_ }
    
    final_df <- .annotate_exonic_region(final_df, gene_config$target_exons)
    if (filter_intronic && nrow(final_df) > 0) final_df <- final_df[is.na(final_df$Region) | final_df$Region == "exonic", ]
    
    final_df$TranscriptRef <- gene_config$transcript
    if (is.na(gene_config$transcript)) final_df$TranscriptRef <- "none"
    
    .write_talos_output(
      final_df, base_name = base_name, output_folder = output_folder,
      sample_name = sample_name, gene_config = gene_config,
      ref_dna = ref_dna, bam_path = bam_path,
      write_vcf = write_vcf, plot = plot, html_report = html_report,
      verbose = verbose, add_config_to_report = add_config_to_report
    )
    
    log_msg(sprintf("Analysis completed successfully. %d variant(s) found.", nrow(final_df)))
    .log_duration(gene_name, sample_name, start_time, verbose)
    return(final_df)
    
  }, error = function(e) {
    err_msg <- sprintf("\n%s\nERROR: %s\n%s", strrep("-", 40), conditionMessage(e), strrep("-", 40))
    log_msg(err_msg)
    stop(e)
  })
}

#' TALOS wrapper with automatic exon fetching
#' @param bam_path BAM file context path.
#' @param gene Gene symbol (e.g., "FLT3")
#' @param build "hg19" or "hg38"
#' @param yaml_path Path to optional YAML config overrides
#' @param padding Base pairs mapped on bounds (default 500L).
#' @param compute_discordant_ratio Boolean: whether to compute discordant pair metric (default TRUE).
#' @param add_config_to_report Include config in plots/HTML report
#' @param max_pairwise_alignments Limit maximum sequences for robust alignment matching
#' @param output_prefix Base name prefix (default "TALOS")
#' @param include_sample Include sample name in filenames
#' @param include_gene Include gene name in filenames
#' @param include_timestamp Include timestamp in filenames
#' @param output_sep Separator for filename parts
#' @param global_log Write all logs to a single file
#' @param min_side_softclip_reads Minimum number of reads with soft‑clip on each side (default 10).
#' @param max_side_ratio Maximum ratio of left/right soft‑clip counts before filtering (default 10).
#' @param min_softclip_pct_side Minimum percentage of supporting reads with soft‑clip on each side (default 1.0).
#' @param ... Extra settings.
#' @export
talos <- function(
    bam_path, gene, build = "hg19", padding = 500L,
    min_support = NULL, min_size = NULL, max_correction = NULL,
    plot = TRUE, sample_name = NULL, filter_intronic = NULL,
    ptd_mode = NULL, use_cigar_bp = NULL, refine_bp = NULL,
    min_mapq = NULL, min_wt_reads = NULL, cluster_tolerance = NULL,
    vaf_threshold = NULL, min_strand_bias = NULL, max_strand_bias = NULL,
    min_mean_support_mapq = NULL, max_breakpoint_spread = NULL,
    min_softclip_fraction = NULL, min_unique_breakpoints = NULL,
    min_coverage_drop = NULL, coverage_window = NULL, min_microhomology = NULL,
    min_discordant_ratio = NULL, min_entropy = NULL, detect_orientation = NULL,
    do_annotate_hotspots = NULL, hotspot_db_path = NULL,
    compute_alignment_score = TRUE, compute_support_bases = TRUE,
    compute_consistency = TRUE, compute_itd_coverage = TRUE,
    compute_coverage_drop = TRUE,
    compute_microhomology = TRUE, compute_repeat_entropy = TRUE,
    compute_discordant_ratio = TRUE,
    compute_hgvs = TRUE,
    min_itd_read_coverage = NULL, html_report = NULL, use_gviz = NULL,
    bsgenome = NULL, yaml_path = system.file("extdata", "gene_config.yaml", package = "TALOS"),
    add_config_to_report = FALSE, max_pairwise_alignments = NULL,
    output_prefix = "TALOS", include_sample = TRUE, include_gene = TRUE,
    include_timestamp = TRUE, output_sep = "_", global_log = TRUE,
    min_side_softclip_reads = 10, max_side_ratio = 10,
    min_softclip_pct_side = 1.0,
    min_left_softclip_pct_wt = 1.0,
    min_right_softclip_pct_wt = 1.0,
    ...

) {
  
  config <- get_gene_config(gene = gene, build = build, padding = padding, config_path = yaml_path, bsgenome = bsgenome)
  config$build <- build
  
  defaults <- list(
    min_support = 10, min_size = 15, max_correction = 2.0, filter_intronic = FALSE,
    ptd_mode = FALSE, use_cigar_bp = TRUE, refine_bp = FALSE, min_mapq = 20,
    min_wt_reads = 0, cluster_tolerance = 10, vaf_threshold = 0.01,
    min_strand_bias = 0, max_strand_bias = 1, min_mean_support_mapq = 0,
    max_breakpoint_spread = Inf, min_softclip_fraction = 0, min_unique_breakpoints = 0,
    min_coverage_drop = 1.0, coverage_window = 200, min_microhomology = 0,
    min_discordant_ratio = 0, min_entropy = 0, detect_orientation = TRUE,
    do_annotate_hotspots = TRUE, compute_alignment_score = TRUE, compute_support_bases = TRUE,
    compute_consistency = TRUE, compute_itd_coverage = TRUE,
    compute_coverage_drop = TRUE, compute_microhomology = TRUE,
    compute_repeat_entropy = TRUE, compute_discordant_ratio = TRUE,
    compute_hgvs = TRUE,
    html_report = TRUE, use_gviz = TRUE, min_itd_read_coverage = 0,
    max_pairwise_alignments = 30L,
    min_side_softclip_reads = 10, max_side_ratio = 10,
    min_softclip_pct_side = 1.0,
    min_left_softclip_pct_wt = 1.0,
    min_right_softclip_pct_wt = 1.0
  )

  
  yaml_vals <- config$gene_settings
  
  resolve <- function(user_val, yaml_val, default_val) {
    if (!is.null(user_val)) return(user_val)
    if (!is.null(yaml_val)) return(yaml_val)
    return(default_val)
  }
  
  p <- lapply(names(defaults), function(nm) {
    resolve(get(nm), yaml_vals[[nm]], defaults[[nm]])
  })
  names(p) <- names(defaults)
  
  if (is.null(config$exons)) config$exons <- config$target_exons
  
  detect_itd(
    bam_path = bam_path, gene_config = config,
    min_support = p$min_support, min_size = p$min_size,
    max_correction = p$max_correction, plot = plot, sample_name = sample_name,
    filter_intronic = p$filter_intronic, ptd_mode = p$ptd_mode,
    use_cigar_bp = p$use_cigar_bp, refine_bp = p$refine_bp,
    min_mapq = p$min_mapq, min_wt_reads = p$min_wt_reads,
    cluster_tolerance = p$cluster_tolerance, vaf_threshold = p$vaf_threshold,
    min_strand_bias = p$min_strand_bias, max_strand_bias = p$max_strand_bias,
    min_mean_support_mapq = p$min_mean_support_mapq, max_breakpoint_spread = p$max_breakpoint_spread,
    min_softclip_fraction = p$min_softclip_fraction, min_unique_breakpoints = p$min_unique_breakpoints,
    min_coverage_drop = p$min_coverage_drop, coverage_window = p$coverage_window,
    min_microhomology = p$min_microhomology, min_discordant_ratio = p$min_discordant_ratio,
    min_entropy = p$min_entropy, detect_orientation = p$detect_orientation,
    do_annotate_hotspots = p$do_annotate_hotspots, hotspot_db_path = hotspot_db_path,
    compute_alignment_score = p$compute_alignment_score, compute_support_bases = p$compute_support_bases,
    compute_consistency = p$compute_consistency, compute_itd_coverage = p$compute_itd_coverage,
    compute_coverage_drop = p$compute_coverage_drop, compute_microhomology = p$compute_microhomology,
    compute_repeat_entropy = p$compute_repeat_entropy, compute_discordant_ratio = p$compute_discordant_ratio,
    compute_hgvs = p$compute_hgvs, html_report = p$html_report, use_gviz = p$use_gviz,
    min_itd_read_coverage = p$min_itd_read_coverage, add_config_to_report = add_config_to_report,
    max_pairwise_alignments = p$max_pairwise_alignments,
    output_prefix = output_prefix, include_sample = include_sample, include_gene = include_gene,
    include_timestamp = include_timestamp, output_sep = output_sep, global_log = global_log,
    min_side_softclip_reads = p$min_side_softclip_reads,
    max_side_ratio          = p$max_side_ratio,
    min_softclip_pct_side   = p$min_softclip_pct_side,
    min_left_softclip_pct_wt = p$min_left_softclip_pct_wt,
    min_right_softclip_pct_wt = p$min_right_softclip_pct_wt,
    ...
  )
}
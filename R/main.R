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
#' @param search_window Max distance (bp) the split-read anchor search
#'   (.estimate_itd_from_anchors()) looks for a re-anchor match. This is a
#'   hard ceiling on any anchor-confirmed length estimate -- a true
#'   duplication longer than this can only be sized via the noisier k-mer
#'   walk fallback. Default 5000 (unchanged for FLT3-scale ITDs). KMT2A
#'   widens this since its PTDs can genuinely exceed 5kb; safe to widen
#'   because the anchor search now requires an unambiguous (single-hit)
#'   match per prefix length rather than pooling every repeat-driven hit.
#' @param prefilter Enable fast prefiltering of alignments.
#' @param min_ins_filter Minimum net insertion length filter.
#' @param output_prefix Base name prefix (default "TALOS").
#' @param output_folder Output folder path.
#' @param write_vcf Generate VCF report.
#' @param verbose Logging enabled.
#' @param debug If TRUE, writes one line per candidate cluster.
#' @param nominal_read_len Assumed read length.
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
#' @param coverage_window Flank (bp) used as the background window for SpanDepthFoldChange.
#' @param min_span_coverage Minimum depth required across the *entire* candidate span
#'   (breakpoint to breakpoint+Length), not just at the breakpoint -- see
#'   compute_span_coverage(). Default 0 (off). This is the filter meant to be
#'   tuned per-gene for large PTDs (KMT2A); short ITDs (FLT3) have spans small
#'   enough that LocalCoverage already covers this case.
#' @param min_span_depth_fold_change Minimum ratio of mean span depth to mean
#'   flanking depth -- see compute_span_depth_fold_change(). Default 0 (off);
#'   soft/informational signal, only enable after checking real numbers since
#'   a low-VAF true event can also fail this.
#' @param max_span_depth_fold_change Maximum ratio of mean span depth to mean
#'   flanking depth. Default Inf (off). A real heterozygous-ish duplication
#'   should show roughly (1+VAF)-fold elevation, not many-fold -- a huge
#'   ratio is Manta's "collapsed reference region" signature (multi-mapping
#'   reads from elsewhere pooling into one locus), not a duplication.
#' @param min_microhomology Required MH length.
#' @param min_discordant_ratio Required ratio length (for filtering).
#' @param min_entropy Sequence repeat complexity baseline.
#' @param min_alignment_score Minimum sequence identity (0–1) between read-derived ITD sequence and the corresponding reference segment. Bypassed automatically for imputed sequences (default 0.2).
#' @param detect_orientation Check inverse loops.
#' @param annotate_hotspots SQLite DB lookups.
#' @param hotspot_db_path Path to SQLite lookup DB.
#' @param compute_alignment_score Metric toggle.
#' @param compute_support_bases Metric toggle.
#' @param compute_consistency Metric toggle. Default FALSE: SupportConsistency
#'   never gates .apply_filters() (purely informational) and is circular for
#'   imputed/ref-derived PTD sequence (checking reads against a sequence
#'   assembled from a subset of those same reads, or filled from reference).
#'   AlignmentScore/RefMatch_Observed/ITDReadCoverage already cover this
#'   ground without the circularity. Set TRUE to restore it for debugging.
#' @param compute_itd_coverage Metric toggle.
#' @param compute_coverage_drop Metric toggle.
#' @param compute_microhomology Metric toggle.
#' @param compute_repeat_entropy Metric toggle.
#' @param compute_discordant_ratio Metric toggle (boolean, default TRUE).
#' @param compute_hgvs Produce HGVS format (only works if transcript is assigned).
#' @param html_report HTML overview toggle.
#' @param use_gviz Use graphical views.
#' @param min_itd_read_coverage Minimum percentage of the ITD length that must be covered by supporting reads (default 50). Low values indicate the sequence was largely imputed rather than directly observed. KMT2A is automatically relaxed to 10 to accommodate large PTDs where only the junctions are read-covered.
#' @param max_reads_in_region Loading safeguard.
#' @param add_config_to_report Render configs in results.
#' @param max_pairwise_alignments Max limit for computationally heavy sub-alignments.
#' @param include_sample Include sample name in output filenames (default TRUE).
#' @param include_gene Include gene name in output filenames (default TRUE).
#' @param include_timestamp Include timestamp in output filenames (default TRUE).
#' @param output_sep Separator between filename parts (default "_").
#' @param global_log Write all logs (config, compute time, errors) to a single file (default TRUE).
#' @param min_side_softclip_reads Minimum number of reads with soft‑clip on each side (default 0, i.e. disabled).
#' @param max_side_ratio Maximum ratio of left/right soft‑clip counts before filtering (default 10).
#' @param min_softclip_pct_side Minimum percentage of supporting reads that must have a soft‑clip on each side (default 0, i.e. disabled).
#' @param min_left_softclip_pct_wt Minimum percentage of wildtype reads with left soft‑clip (default 0, i.e. disabled).
#' @param min_right_softclip_pct_wt Minimum percentage of wildtype reads with right soft‑clip (default 0, i.e. disabled).
#' @param min_abs_side_softclip Minimum absolute soft‑clip count on each side (default 0, i.e. disabled).
#' @param max_itd_length Maximum duplication length to consider as ITD; longer duplications are automatically converted to PTD (zero length) if \code{convert_long_to_ptd} is TRUE (default 1000).
#' @param convert_long_to_ptd If TRUE, duplications longer than \code{max_itd_length} are reported as PTDs (length 0). If FALSE, they are skipped (default TRUE).
#' @param min_length Minimum duplication length to report (NULL = no lower bound).  
#' @param max_length Maximum duplication length to report (NULL = no upper bound).  
#' @param use_kmers Enable k-mer based analysis (default TRUE).
#' @param exon_padding Number of flanking exons to add to the genomic window (default 0).
#' @param prefer_extended_length If TRUE and duplication length > nominal_read_len,
#'        use extension-based ITD length (mode of itdsize_ext) and refined breakpoint
#'        as the primary call. Default TRUE.
#' @param ptd_allow_asymmetric Allow PTD detection with soft-clips on only one side (default TRUE).
#' @param ptd_use_local_assembly Attempt to reconstruct duplication from soft-clips (default TRUE).
#' @param merge_ptd_intervals If TRUE, merge overlapping/adjacent large duplication intervals into a single PTD call (default FALSE).
#' @param merge_gap Maximum gap (bp) between intervals to merge (default 500).
#' @param min_ptd_length Minimum length (bp) for intervals to be considered for merging (default 200).
#' @param max_plausible_fragments If a merged interval pools more than this many
#'        eligible fragments, it is flagged (\code{MergeFragmentWarning = TRUE})
#'        rather than silently reported at full confidence: a single true PTD
#'        has historically fragmented into at most ~2-3 rows pre-merge, so a
#'        larger pool is more consistent with independent scattered signal
#'        (e.g. repeat-driven misalignment) getting merged by genomic
#'        proximity alone. Informational only -- does not remove rows from
#'        the result. Default 4.
#' @param use_exon_graph If TRUE, additionally run the directed exon-graph PTD
#'        reconstruction (requires exon annotations in gene_config$all_exons).
#'        Graph-derived candidates supplement rather than replace the legacy
#'        candidate extraction below; the two are combined and reconciled by
#'        the merge_ptd_intervals step. No longer requires ptd_mode=TRUE.
#'        Default FALSE.
#' @param ... Other parameters passed to the underlying engine.
#' @return Data frame of ITD calls with diagnostic columns, invisibly.
#' @export
detect_itd <- function(
    bam_path, gene_config, genomic_ref_seq = NULL, k = 11, min_support = 10,
    min_size = 15, max_missing_kmers = 0.5, cluster_tolerance = 10, prefilter = TRUE,
    search_window = 5000L,
    min_ins_filter = 3, output_prefix = "TALOS", output_folder = "./results",
    write_vcf = TRUE, verbose = TRUE, debug = FALSE, nominal_read_len = 150,
    max_correction = 2.0, plot = FALSE, sample_name = NULL, filter_intronic = FALSE,
    ptd_mode = FALSE, use_cigar_bp = TRUE, refine_bp = FALSE, min_mapq = 20,
    min_wt_reads = 0, vaf_threshold = 0.01, min_strand_bias = 0, max_strand_bias = 1,
    min_mean_support_mapq = 0, max_breakpoint_spread = Inf, min_softclip_fraction = 0,
    min_unique_breakpoints = 0, min_coverage_drop = 1.0, min_local_coverage = 0, coverage_window = 200,
    min_span_coverage = 0, min_span_depth_fold_change = 0, max_span_depth_fold_change = Inf,
    min_microhomology = 0, min_discordant_ratio = 0, min_entropy = 0,
    min_alignment_score = 0.2,
    detect_orientation = TRUE, do_annotate_hotspots = TRUE, hotspot_db_path = NULL,
    compute_alignment_score = TRUE, compute_support_bases = TRUE,
    compute_consistency = FALSE, compute_itd_coverage = TRUE,
    compute_coverage_drop = TRUE, compute_microhomology = TRUE,
    compute_repeat_entropy = TRUE, compute_discordant_ratio = TRUE,
    compute_hgvs = TRUE, html_report = TRUE, use_gviz = TRUE,
    min_itd_read_coverage = 50, max_reads_in_region = 200000,
    add_config_to_report = FALSE, max_pairwise_alignments = 30L,
    include_sample = TRUE, include_gene = TRUE, include_timestamp = TRUE,
    output_sep = "_", global_log = TRUE,
    min_side_softclip_reads = 0, max_side_ratio = 10,
    min_softclip_pct_side = 0,
    min_left_softclip_pct_wt = 0,
    min_right_softclip_pct_wt = 0,
    min_abs_side_softclip = 0,
    max_itd_length = 1000,
    convert_long_to_ptd = TRUE,
    min_length = NULL,     
    max_length = NULL,
    use_kmers = TRUE,
    exon_padding = 0L,
    prefer_extended_length = TRUE,
    ptd_allow_asymmetric = TRUE,
    ptd_use_local_assembly = TRUE,
    merge_ptd_intervals = FALSE,      
    merge_gap = 500L,                 
    min_ptd_length = 200L,
    max_plausible_fragments = 4L,
    use_exon_graph = FALSE,   
    ...
) {
  
  start_time <- Sys.time()
  gene_name  <- if (!is.null(gene_config$gene)) gene_config$gene else "UNKNOWN"
  if (is.null(sample_name)) sample_name <- sub("\\..*", "", basename(bam_path))
  
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  sample_folder <- file.path(output_folder, sample_name)
  if (!dir.exists(sample_folder)) dir.create(sample_folder, recursive = TRUE, showWarnings = FALSE)
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  base_parts <- c()
  if (!is.null(output_prefix) && output_prefix != "") base_parts <- c(base_parts, output_prefix)
  if (include_sample && !is.null(sample_name) && sample_name != "") base_parts <- c(base_parts, sample_name)
  if (include_gene   && !is.null(gene_name)   && gene_name != "")   base_parts <- c(base_parts, gene_name)
  if (include_timestamp && !is.null(timestamp) && timestamp != "")  base_parts <- c(base_parts, timestamp)
  base_name <- paste(base_parts, collapse = output_sep)
  
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
                                          max_reads = max_reads_in_region,
                                          verbose = verbose)
    all_reads <- bam_data$reads
    
    if (length(all_reads) == 0) {
      log_msg("No reads found in target region.")
      return(data.frame())
    }
    
    wt_mapqs <- S4Vectors::mcols(all_reads)$mapq
    wt_mapqs[is.na(wt_mapqs)] <- 0L
    all_reads_for_wt <- all_reads[wt_mapqs >= min_mapq]
    
    if (prefilter) {
      all_reads <- .filter_reads_by_cigar(all_reads, min_mapq, min_ins_filter, verbose = verbose)
      if (length(all_reads) == 0) {
        log_msg("No reads passed CIGAR prefilter.")
        return(data.frame())
      }
    }
    
    candidates <- list()
    bp_df <- data.frame()
    
    if (use_exon_graph && !is.null(gene_config$all_exons)) {
      if (verbose) message("[TALOS] Running Exon-Graph PTD reconstruction...")
      
      if (!exists("extract_graph_candidates", mode = "function")) {
        stop("extract_graph_candidates() not found. Please source engine_graph_ptd.R")
      }
      
      graph_cycles <- extract_graph_candidates(bam_path, gene_config)
      
      if (length(graph_cycles) > 0) {
        graph_candidates <- list()
        for (cand in graph_cycles) {
          support_qnames <- unlist(strsplit(cand$support_qnames, ";"))
          junction_reads <- all_reads[.safe_qnames(all_reads) %in% support_qnames]
          
          if (length(junction_reads) == 0) next
          
          junction_bps <- c(BiocGenerics::start(junction_reads), BiocGenerics::end(junction_reads))
          refined_bp <- .kde_breakpoint(junction_bps)
          local_bp   <- refined_bp - gene_config$genomic_start + 1L
          
          sc_seqs <- sapply(seq_len(length(junction_reads)), function(i) {
            sc <- .get_softclips(GenomicAlignments::cigar(junction_reads)[i], 
                                 as.character(S4Vectors::mcols(junction_reads)$seq[i]))
            paste(na.omit(c(sc$lead, sc$trail)), collapse="")
          })
          mapqs <- S4Vectors::mcols(junction_reads)$mapq
          
          junction_seq <- .build_local_debruijn(sc_seqs)
          if (is.na(junction_seq)) {
            junction_seq <- .assemble_weighted_consensus(sc_seqs, mapqs)
          }
          
          read_qnames  <- .safe_qnames(junction_reads)
          read_mapqs   <- S4Vectors::mcols(junction_reads)$mapq; read_mapqs[is.na(read_mapqs)] <- 0L
          read_flags   <- S4Vectors::mcols(junction_reads)$flag
          read_reverse <- bitwAnd(read_flags, 0x10L) != 0L
          read_cigars  <- GenomicAlignments::cigar(junction_reads)
          
          for (r in seq_len(length(junction_reads))) {
            graph_candidates[[length(graph_candidates) + 1]] <- list(
              read_name         = read_qnames[r],
              local_breakpoint  = local_bp,
              length            = 0L, 
              type              = "ptd_graph",
              mapq              = read_mapqs[r],
              is_reverse        = read_reverse[r],
              cigar             = read_cigars[r],
              read_seq          = junction_seq,
              itd_seq           = junction_seq
            )
          }
        }
        if (length(graph_candidates) > 0) {
          candidates <- c(candidates, graph_candidates)
          if (verbose) message("[TALOS] Graph extraction contributed ", length(graph_candidates),
                               " read-level PTD candidate(s) from ", length(graph_cycles), " cycle(s).")
        }
      } else {
        log_msg("Graph extraction found no PTD cycles.")
      }
    }
    
    legacy_candidates <- list()
    if (ptd_mode && use_cigar_bp) {
      legacy_candidates <- .extract_candidates_ptd(all_reads, gene_config$genomic_start, ref_len, verbose = verbose)
      if (use_kmers) {
        jump_candidates <- .extract_candidates_standard(
          reads = all_reads, ref_kmers = ref_kmers, ptd_mode = TRUE, min_size = min_size,
          max_missing_kmers = max_missing_kmers, refine_bp = refine_bp, use_cigar_bp = use_cigar_bp,
          genomic_start = gene_config$genomic_start, ref_dna = ref_dna,
          max_itd_length = max_itd_length, convert_long_to_ptd = convert_long_to_ptd,
          search_window = search_window, verbose = verbose
        )
        legacy_candidates <- c(legacy_candidates, jump_candidates)
      }
    } else if (!use_kmers) {
      legacy_candidates <- .extract_candidates_cigar(all_reads, gene_config$genomic_start, ref_len,
                                              min_size = min_size, verbose = verbose)
    } else {
      legacy_candidates <- .extract_candidates_standard(
        reads = all_reads, ref_kmers = ref_kmers, ptd_mode = ptd_mode, min_size = min_size,
        max_missing_kmers = max_missing_kmers, refine_bp = refine_bp, use_cigar_bp = use_cigar_bp,
        genomic_start = gene_config$genomic_start, ref_dna = ref_dna,
        max_itd_length = max_itd_length, convert_long_to_ptd = convert_long_to_ptd,
        search_window = search_window, verbose = verbose
      )
    }
    candidates <- c(candidates, legacy_candidates)
    
    if (length(candidates) == 0) {
      log_msg("No candidate reads (clips or k-mer jumps) found.")
      return(data.frame())
    }
    bp_df <- .candidates_to_df(candidates, gene_config$genomic_start)
    
    clusters <- .cluster_breakpoints(bp_df$breakpoint, cluster_tolerance)
    
    ptd_assembled_seqs <- list()
    if (ptd_mode && ptd_use_local_assembly && nrow(bp_df) > 0) {
      for (cl in clusters) {
        cl_df <- bp_df[bp_df$breakpoint %in% cl, ]
        assembled <- .assemble_ptd_consensus(cl_df, min_reads = 3L, min_len = 15L)
        if (!is.na(assembled) && nchar(assembled) > 0) {
          med_cl <- as.character(round(median(cl)))
          ptd_assembled_seqs[[med_cl]] <- assembled
          log_msg(sprintf("[TALOS] Assembled PTD consensus of length %d bp for cluster at %s",
                          nchar(assembled), med_cl))
        }
      }
    }
    
    .lookup_assembled_seq <- function(bp, store, tol = 50L) {
      if (length(store) == 0L) return(NA_character_)
      keys <- suppressWarnings(as.integer(names(store)))
      dists <- abs(keys - as.integer(bp))
      best  <- which.min(dists)
      if (dists[best] <= tol) store[[best]] else NA_character_
    }
    
    ext_df <- .extend_candidates(
      bp_df[!is.na(bp_df$length) & bp_df$length > 0L, ],
      ref_dna       = ref_dna,
      genomic_start = gene_config$genomic_start,
      min_size      = min_size,
      search_window = search_window,
      verbose       = verbose
    )
    
    wt_info <- .prepare_wildtype_info(all_reads_for_wt, gene_config$genomic_start, gene_config$genomic_end)
    wt_info$sample_name <- sample_name
    
    thresholds <- list(
      min_coverage_drop = min_coverage_drop, min_local_coverage = min_local_coverage,
      min_span_coverage = min_span_coverage,
      min_span_depth_fold_change = min_span_depth_fold_change,
      max_span_depth_fold_change = max_span_depth_fold_change,
      min_microhomology = min_microhomology,
      min_discordant_ratio = min_discordant_ratio, min_entropy = min_entropy,
      min_strand_bias = min_strand_bias, max_strand_bias = max_strand_bias,
      min_mean_support_mapq = min_mean_support_mapq, max_breakpoint_spread = max_breakpoint_spread,
      min_softclip_fraction = min_softclip_fraction, min_unique_breakpoints = min_unique_breakpoints,
      min_itd_read_coverage = min_itd_read_coverage,
      min_alignment_score   = min_alignment_score,
      min_side_softclip_reads = min_side_softclip_reads,
      max_side_ratio          = max_side_ratio,
      min_softclip_pct_side   = min_softclip_pct_side,
      min_left_softclip_pct_wt = min_left_softclip_pct_wt,
      min_right_softclip_pct_wt = min_right_softclip_pct_wt,
      min_abs_side_softclip    = min_abs_side_softclip  
    )
    
    results <- list()
    
    if (prefer_extended_length && !ptd_mode && nrow(ext_df) >= 3) {
      ext_valid <- ext_df[!is.na(ext_df$itdsize) & ext_df$itdsize >= min_size, ]
      if (nrow(ext_valid) >= 3) {
        refined_clusters <- .cluster_breakpoints(ext_valid$breakpoint_refined, cluster_tolerance)
        for (rcl in refined_clusters) {
          rcl_reads <- ext_valid[ext_valid$breakpoint_refined %in% rcl, ]
          if (nrow(rcl_reads) < 3) {
            if (debug) log_msg(sprintf(
              "[CANDIDATE] GenomicPosition=%s SupportingReads=%d -> SKIPPED (refined cluster has <3 supporting reads; metrics not computed, min_support/other filters never applied)",
              as.character(round(median(rcl, na.rm = TRUE))), nrow(rcl_reads)))
            next
          }
          
          anchor_vals  <- rcl_reads$itdsize[rcl_reads$itdsize_source == "anchor"]
          best_len     <- if (length(anchor_vals) > 0L) .mode(anchor_vals) else .mode(rcl_reads$itdsize)
          median_bp <- as.integer(median(rcl_reads$breakpoint_refined, na.rm = TRUE))
          support_read_names <- rcl_reads$read_name
          support_rows <- bp_df[bp_df$read_name %in% support_read_names &
                                abs(bp_df$breakpoint - median_bp) <= cluster_tolerance, ]

          additional <- bp_df[abs(bp_df$breakpoint - median_bp) <= cluster_tolerance & 
                              !(bp_df$read_name %in% support_read_names), ]
          support_rows <- rbind(support_rows, additional)
          support_rows <- unique(support_rows)
          
          if (nrow(support_rows) == 0) {
            if (debug) log_msg(sprintf(
              "[CANDIDATE] GenomicPosition=%d SupportingReads=0 -> SKIPPED (no bp_df rows matched refined cluster; metrics not computed)",
              median_bp))
            next
          }
          
          length_ext <- best_len
          if (verbose) {
            log_msg(sprintf("Extension-based cluster at %d: %d bp (mode of %d anchor-confirmed / %d total values, %d supporting reads)",
                            median_bp, best_len, length(anchor_vals), nrow(rcl_reads), nrow(support_rows)))
          }
          metrics <- .compute_variant_metrics(
            cluster_bps = rep(median_bp, nrow(support_rows)), best_len = best_len, support_rows = support_rows,
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
            span_flank = coverage_window,
            length_ext = length_ext,
            pre_assembled_seq = .lookup_assembled_seq(median_bp, ptd_assembled_seqs),
            debug = debug, verbose = verbose
          )
          if (!is.null(metrics)) {
            passed <- .apply_filters(metrics, thresholds, min_length, max_length,
                                      ptd_mode = ptd_mode, ptd_allow_asymmetric = ptd_allow_asymmetric)
            if (debug) log_msg(.format_candidate_debug(metrics, passed))
            if (passed) results[[length(results) + 1L]] <- metrics
          }
        }
      }
    }
    
    if (length(results) == 0) {
      ext_by_read <- if (nrow(ext_df) > 0) {
        tmp <- tapply(ext_df$itdsize, ext_df$read_name, function(x) x[1L])
        tmp[!is.na(tmp)]
      } else {
        c()
      }
      bp_refined_by_read <- if (nrow(ext_df) > 0) {
        tmp <- tapply(ext_df$breakpoint_refined, ext_df$read_name, function(x) x[1L])
        tmp[!is.na(tmp)]
      } else {
        c()
      }
      
      cand_breakpoints <- bp_df$breakpoint
      cand_lengths     <- bp_df$length
      clusters <- .cluster_breakpoints(bp_df$breakpoint, cluster_tolerance)
      
      for (cl in clusters) {
        cluster_candidates <- which(cand_breakpoints %in% cl)
        if (length(cluster_candidates) == 0) {
          if (debug) log_msg(sprintf(
            "[CANDIDATE] GenomicPosition=%s SupportingReads=0 -> SKIPPED (no candidates matched cluster; metrics not computed)",
            as.character(round(median(cl, na.rm = TRUE)))))
          next
        }
        
        ext_sizes <- c()
        refined_bps <- c()
        support_read_names <- bp_df$read_name[cluster_candidates]
        for (rn in support_read_names) {
          sz <- ext_by_read[rn]
          if (length(sz) == 0) sz <- NA_integer_
          else if (length(sz) > 1) sz <- sz[1L]
          bp_ref <- bp_refined_by_read[rn]
          if (length(bp_ref) == 0) bp_ref <- NA_integer_
          else if (length(bp_ref) > 1) bp_ref <- bp_ref[1L]
          
          if (!is.na(sz) && !is.na(bp_ref)) {
            ext_sizes <- c(ext_sizes, sz)
            refined_bps <- c(refined_bps, bp_ref)
          }
        }
        
        use_ext <- prefer_extended_length && !ptd_mode && length(ext_sizes) >= 3
        if (use_ext) {
          best_len <- .mode(ext_sizes)
          median_bp <- as.integer(median(refined_bps, na.rm = TRUE))
          cl <- rep(median_bp, length(cl))
          support_rows <- bp_df[bp_df$read_name %in% names(ext_by_read), ]
          support_rows <- support_rows[abs(support_rows$breakpoint - median_bp) <= cluster_tolerance, ]
          if (nrow(support_rows) == 0) support_rows <- bp_df[bp_df$breakpoint %in% cl, ]
          length_ext <- best_len
          metrics <- .compute_variant_metrics(
            cluster_bps = rep(median_bp, nrow(support_rows)), best_len = best_len, support_rows = support_rows,
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
            span_flank = coverage_window,
            length_ext = length_ext,
            pre_assembled_seq = .lookup_assembled_seq(median_bp, ptd_assembled_seqs),
            debug = debug, verbose = verbose
          )
          if (!is.null(metrics)) {
            passed <- .apply_filters(metrics, thresholds, min_length, max_length,
                                      ptd_mode = ptd_mode, ptd_allow_asymmetric = ptd_allow_asymmetric)
            if (debug) log_msg(.format_candidate_debug(metrics, passed))
            if (passed) results[[length(results) + 1L]] <- metrics
          }
          next
        }
        
        # Fallback to k‑mer based length grouping
        cluster_lengths <- cand_lengths[cluster_candidates]
        cluster_lengths[is.na(cluster_lengths)] <- -1L
        len_groups <- split(cluster_candidates, cluster_lengths)
        for (len in names(len_groups)) {
          best_len <- if (as.integer(len) == -1L) NA_integer_ else as.integer(len)
          idxs <- len_groups[[len]]
          if (is.na(best_len)) support_rows <- bp_df[bp_df$breakpoint %in% cl & is.na(bp_df$length), ]
          else support_rows <- bp_df[bp_df$breakpoint %in% cl & !is.na(bp_df$length) & bp_df$length == best_len, ]
          
          ext_sizes_cl <- unlist(ext_by_read[support_rows$read_name], use.names = FALSE)
          ext_sizes_cl <- ext_sizes_cl[!is.na(ext_sizes_cl) & ext_sizes_cl > 0L]
          length_ext <- if (length(ext_sizes_cl) >= 3L) {
            as.integer(names(sort(table(ext_sizes_cl), decreasing = TRUE))[1L])
          } else NA_integer_
          
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
            span_flank = coverage_window,
            length_ext = length_ext,
            pre_assembled_seq = .lookup_assembled_seq(round(median(cl)), ptd_assembled_seqs),
            debug = debug, verbose = verbose
          )
          if (!is.null(metrics)) {
            passed <- .apply_filters(metrics, thresholds, min_length, max_length,
                                      ptd_mode = ptd_mode, ptd_allow_asymmetric = ptd_allow_asymmetric)
            if (debug) log_msg(.format_candidate_debug(metrics, passed))
            if (passed) results[[length(results) + 1L]] <- metrics
          }
        }
      }
    }
    
    if (length(results) == 0) {
      log_msg("No variants passed all filters.")
      return(data.frame())
    }
    
    final_df <- do.call(rbind, lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
    final_df$MergedFromN <- 1L   
    final_df$MergeFragmentWarning <- FALSE  
    final_df$Hotspot     <- FALSE
    final_df$HotspotName <- NA_character_
    final_df$Region      <- NA_character_
    final_df$ExonNumber  <- NA_integer_
    if (merge_ptd_intervals && nrow(final_df) > 0) {
      keep <- which(!is.na(final_df$Length) & final_df$Length >= min_ptd_length)
      skip <- setdiff(seq_len(nrow(final_df)), keep)  
      if (length(keep) >= 2) {
        log_msg(sprintf("Merging %d large duplication intervals (gap \u2264 %d bp)...",
                        length(keep), merge_gap))
        gr <- GenomicRanges::GRanges(
          seqnames = final_df$Genome[keep],
          ranges = IRanges::IRanges(
            start = final_df$GenomicPosition[keep],
            end   = final_df$GenomicPosition[keep] + final_df$Length[keep] - 1L
          )
        )
        merged_gr <- GenomicRanges::reduce(gr, min.gapwidth = merge_gap)
        if (length(merged_gr) > 0) {
          hits <- GenomicRanges::findOverlaps(gr, merged_gr)
          group_of <- integer(length(gr))
          group_of[S4Vectors::queryHits(hits)] <- S4Vectors::subjectHits(hits)
          
          merged_rows <- lapply(seq_along(merged_gr), function(i) {
            grp_keep <- keep[group_of == i]  
            n_members <- length(grp_keep)
            mg_start <- GenomicRanges::start(merged_gr[i])
            mg_len   <- GenomicRanges::width(merged_gr[i])
            tmpl_idx <- grp_keep[which.min(final_df$GenomicPosition[grp_keep])]
            template <- final_df[tmpl_idx, ]
            support_w <- final_df$SupportingReads[grp_keep]
            support_w[is.na(support_w) | support_w < 0] <- 0
            wmean <- function(vals, digits = NA_integer_) {
              ok <- !is.na(vals) & support_w > 0
              if (!any(ok)) return(NA_real_)
              out <- stats::weighted.mean(vals[ok], support_w[ok])
              if (!is.na(digits)) out <- round(out, digits)
              out
            }
            
            template$GenomicPosition  <- mg_start
            template$Length           <- mg_len
            if (!is.null(bam_data$cov)) {
              mg_span_cov <- compute_span_coverage(bam_data$cov, gene_config$chrom,
                                                   mg_start, mg_len, verbose = debug)
              template$SpanMinCoverage  <- mg_span_cov$min_cov
              template$SpanMeanCoverage <- mg_span_cov$mean_cov
              template$SpanDepthFoldChange <- compute_span_depth_fold_change(
                bam_data$cov, gene_config$chrom, mg_start, mg_len,
                flank = coverage_window, verbose = debug
              )
            }
            template$SupportingReads  <- sum(final_df$SupportingReads[grp_keep], na.rm = TRUE)
            template$WildtypeReads    <- sum(final_df$WildtypeReads[grp_keep],   na.rm = TRUE)
            template$DepthAtBreakpoint <- template$SupportingReads + template$WildtypeReads
            template$AlleleFrequency  <- if (template$DepthAtBreakpoint > 0L)
              template$SupportingReads / template$DepthAtBreakpoint else NA_real_

            template$TotalSupportBases  <- as.integer(sum(final_df$TotalSupportBases[grp_keep],  na.rm = TRUE))
            template$LeftSoftclipCount  <- as.integer(sum(final_df$LeftSoftclipCount[grp_keep],   na.rm = TRUE))
            template$RightSoftclipCount <- as.integer(sum(final_df$RightSoftclipCount[grp_keep],  na.rm = TRUE))
            template$BothSoftclipCount  <- as.integer(sum(final_df$BothSoftclipCount[grp_keep],   na.rm = TRUE))

            template$LeftSoftclipPctSupport  <- if (template$SupportingReads > 0L)
              template$LeftSoftclipCount  / template$SupportingReads * 100 else NA_real_
            template$RightSoftclipPctSupport <- if (template$SupportingReads > 0L)
              template$RightSoftclipCount / template$SupportingReads * 100 else NA_real_
            template$LeftSoftclipPctWT  <- if (template$WildtypeReads > 0L)
              template$LeftSoftclipCount  / template$WildtypeReads * 100 else NA_real_
            template$RightSoftclipPctWT <- if (template$WildtypeReads > 0L)
              template$RightSoftclipCount / template$WildtypeReads * 100 else NA_real_
            

            template$RefMatch_Observed  <- wmean(final_df$RefMatch_Observed[grp_keep], 1L)
            template$RefMatch_Total     <- wmean(final_df$RefMatch_Total[grp_keep],    1L)
            template$ITDReadCoverage    <- wmean(final_df$ITDReadCoverage[grp_keep])
            template$SupportConsistency <- wmean(final_df$SupportConsistency[grp_keep])
            template$AlignmentScore     <- wmean(final_df$AlignmentScore[grp_keep])
            template$StrandBias         <- wmean(final_df$StrandBias[grp_keep], 4L)
            template$MeanSupportMAPQ    <- wmean(final_df$MeanSupportMAPQ[grp_keep], 1L)

            if (n_members > 1L) {
              template$ITD_Sequence       <- NA_character_
              template$HGVS_cDNA          <- NA_character_
              template$HGVS_Protein       <- NA_character_
              template$ITDCoverageRLE     <- NA_character_
              template$LengthPE           <- NA_real_
              template$LengthPE_NSpanning <- NA_integer_
              template$LengthExt          <- NA_integer_
              template$SequenceImputed    <- TRUE
              template$SequencePartial    <- TRUE
              template$SequenceSource     <- "merged_intervals"
              template$Hotspot     <- FALSE
              template$HotspotName <- NA_character_
              template$Region      <- NA_character_
              template$ExonNumber  <- NA_integer_
            }
            template$PESoftclipSupport   <- sum(final_df$PESoftclipSupport[grp_keep], na.rm = TRUE)
            template$PESoftclipEventPairs <- sum(final_df$PESoftclipEventPairs[grp_keep], na.rm = TRUE)
            template$PESoftclipLongPairs <- sum(final_df$PESoftclipLongPairs[grp_keep], na.rm = TRUE)
            template$PEOrientationFR     <- sum(final_df$PEOrientationFR[grp_keep], na.rm = TRUE)
            template$PEOrientationRF     <- sum(final_df$PEOrientationRF[grp_keep], na.rm = TRUE)
            template$PEOrientationFF     <- sum(final_df$PEOrientationFF[grp_keep], na.rm = TRUE)
            template$PEOrientationRR     <- sum(final_df$PEOrientationRR[grp_keep], na.rm = TRUE)
            template$PEOrientationOther  <- sum(final_df$PEOrientationOther[grp_keep], na.rm = TRUE)
            pe_counts <- c(FR = template$PEOrientationFR, RF = template$PEOrientationRF,
                          FF = template$PEOrientationFF, RR = template$PEOrientationRR,
                          Other = template$PEOrientationOther)
            template$PEOrientationDominant <- if (all(pe_counts == 0)) NA_character_
                                              else names(pe_counts)[which.max(pe_counts)]
            template$MergedFromN <- n_members

            template$MergeFragmentWarning <- n_members > max_plausible_fragments
            if (isTRUE(template$MergeFragmentWarning)) {
              log_msg(sprintf(
                paste("WARNING: merged interval at %d (length %d bp) pools %d fragments,",
                      "exceeding max_plausible_fragments=%d. A single true PTD has",
                      "historically fragmented into at most ~2-3 rows pre-merge, so this",
                      "many disjoint clusters merging by proximity alone is more consistent",
                      "with scattered background signal than one real event -- review before",
                      "reporting. Flagged via MergeFragmentWarning."),
                mg_start, mg_len, n_members, max_plausible_fragments))
            }

            if (!is.na(template$SpanMinCoverage) &&
                template$SpanMinCoverage < thresholds$min_span_coverage) {
              log_msg(sprintf(
                paste("REJECTED merged interval at %d (length %d bp, %d fragment(s)):",
                      "SpanMinCoverage=%.1f < min_span_coverage=%.1f -- part of this span",
                      "has essentially no real coverage, so the merge is not trusted as",
                      "one true event."),
                mg_start, mg_len, n_members, template$SpanMinCoverage, thresholds$min_span_coverage))
              return(NULL)
            }
            
            template
          })
          merged_rows <- Filter(Negate(is.null), merged_rows)
          if (length(merged_rows) > 0) {
            merged_df <- do.call(rbind, merged_rows)
            log_msg(sprintf("Merged into %d interval(s).", nrow(merged_df)))
            final_df <- if (length(skip) > 0L)
              rbind(merged_df, final_df[skip, , drop = FALSE])
            else
              merged_df
          } else {
            log_msg("All merged intervals rejected by min_span_coverage; none survive.")
            final_df <- final_df[skip, , drop = FALSE]
          }
        } else {
          log_msg("No intervals could be merged.")
        }
      } else if (length(keep) == 1L) {
        log_msg("Only one large duplication found; no merging performed.")
      } else {
        log_msg(sprintf("No duplications >= %d bp found; skipping merge.", min_ptd_length))
      }
    }

    
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
    return(final_df)
    
  }, error = function(e) {
    err_msg <- sprintf("\n%s\nERROR: %s\n%s", strrep("-", 40), conditionMessage(e), strrep("-", 40))
    log_msg(err_msg)
    stop(e)
  })
}

==========================================================================

#' @rdname detect_itd
#' @export
talos <- function(
    bam_path, gene, build = "hg19", padding = 500L,
    min_support = NULL, min_size = NULL, max_correction = NULL,
    max_missing_kmers = NULL,
    plot = TRUE, sample_name = NULL, filter_intronic = NULL,
    ptd_mode = NULL, use_cigar_bp = NULL, refine_bp = NULL,
    min_mapq = NULL, min_wt_reads = NULL, cluster_tolerance = NULL,
    search_window = NULL,
    vaf_threshold = NULL, min_strand_bias = NULL, max_strand_bias = NULL,
    min_mean_support_mapq = NULL, max_breakpoint_spread = NULL,
    min_softclip_fraction = NULL, min_unique_breakpoints = NULL,
    min_coverage_drop = NULL, min_local_coverage = NULL, coverage_window = NULL, min_microhomology = NULL,
    min_span_coverage = NULL, min_span_depth_fold_change = NULL, max_span_depth_fold_change = NULL,
    min_discordant_ratio = NULL, min_entropy = NULL, min_alignment_score = NULL,
    detect_orientation = NULL,
    do_annotate_hotspots = NULL, hotspot_db_path = NULL,
    compute_alignment_score = TRUE, compute_support_bases = TRUE,
    compute_consistency = FALSE, compute_itd_coverage = TRUE,
    compute_coverage_drop = TRUE,
    compute_microhomology = TRUE, compute_repeat_entropy = TRUE,
    compute_discordant_ratio = TRUE,
    compute_hgvs = TRUE,
    min_itd_read_coverage = NULL, html_report = NULL, use_gviz = NULL,
    bsgenome = NULL, yaml_path = system.file("extdata", "gene_config.yaml", package = "TALOS"),
    add_config_to_report = FALSE, max_pairwise_alignments = NULL,
    output_prefix = "TALOS", include_sample = TRUE, include_gene = TRUE,
    include_timestamp = TRUE, output_sep = "_", global_log = TRUE,
    min_side_softclip_reads = NULL, max_side_ratio = NULL,
    min_softclip_pct_side = NULL,
    min_left_softclip_pct_wt = NULL,
    min_right_softclip_pct_wt = NULL,
    min_abs_side_softclip = NULL,
    max_itd_length = NULL,
    convert_long_to_ptd = NULL,
    min_length = NULL,       
    max_length = NULL,
    use_kmers = TRUE,
    exon_padding = 0L,
    prefer_extended_length = NULL,
    ptd_allow_asymmetric = NULL,
    ptd_use_local_assembly = NULL,
    merge_ptd_intervals = NULL,      
    merge_gap = NULL,                
    min_ptd_length = NULL,           
    max_plausible_fragments = NULL,
    use_exon_graph = NULL,   
    verbose = TRUE,
    debug = FALSE,
    ...
) {
  
  config <- get_gene_config(gene = gene, build = build, padding = padding, 
                            config_path = yaml_path, bsgenome = bsgenome,
                            exon_padding = exon_padding)
  config$build <- build
  

  defaults <- list(
    min_support = 10, min_size = 15, max_correction = 2.0, filter_intronic = FALSE,
    max_missing_kmers = 0.5, ptd_mode = FALSE, use_cigar_bp = TRUE, refine_bp = FALSE, min_mapq = 20,
    min_wt_reads = 0, cluster_tolerance = 10, vaf_threshold = 0.01,
    search_window = 5000L,
    min_strand_bias = 0, max_strand_bias = 1, min_mean_support_mapq = 0,
    max_breakpoint_spread = Inf, min_softclip_fraction = 0, min_unique_breakpoints = 0,
    min_coverage_drop = 1.0, min_local_coverage = 0, coverage_window = 200, min_microhomology = 0,
    min_span_coverage = 0, min_span_depth_fold_change = 0, max_span_depth_fold_change = Inf,
    min_discordant_ratio = 0, min_entropy = 0, min_alignment_score = 0.2,
    detect_orientation = TRUE,
    do_annotate_hotspots = TRUE, compute_alignment_score = TRUE, compute_support_bases = TRUE,
    compute_consistency = FALSE, compute_itd_coverage = TRUE,
    compute_coverage_drop = TRUE, compute_microhomology = TRUE,
    compute_repeat_entropy = TRUE, compute_discordant_ratio = TRUE,
    compute_hgvs = TRUE,
    html_report = TRUE, use_gviz = TRUE, min_itd_read_coverage = 50,
    max_pairwise_alignments = 30L,
    min_side_softclip_reads = 0, max_side_ratio = 10,
    min_softclip_pct_side = 0,
    min_left_softclip_pct_wt = 0,
    min_right_softclip_pct_wt = 0,
    min_abs_side_softclip = 0,
    max_itd_length = 1000,
    convert_long_to_ptd = TRUE,
    min_length = NULL,       
    max_length = NULL,
    prefer_extended_length = TRUE,
    ptd_allow_asymmetric = TRUE,
    ptd_use_local_assembly = TRUE,
    merge_ptd_intervals = FALSE,
    merge_gap = 500L,
    min_ptd_length = 200L,
    max_plausible_fragments = 4L,
    use_exon_graph = FALSE
  )
  
  yaml_vals <- config$gene_settings
  

  resolve <- function(user_val, yaml_val, default_val) {
    if (!is.null(user_val)) return(user_val)
    if (!is.null(yaml_val)) return(yaml_val)
    return(default_val)
  }
  
  p <- lapply(names(defaults), function(nm) {
    user_val <- get(nm)  
    yaml_val <- yaml_vals[[nm]]
    resolve(user_val, yaml_val, defaults[[nm]])
  })
  names(p) <- names(defaults)
  
  if (gene == "KMT2A") {
    if (verbose) message("[TALOS] Applying relaxed thresholds for KMT2A PTD detection")
    p$min_support             <- min(p$min_support, 5L)
    p$min_side_softclip_reads <- min(p$min_side_softclip_reads, 10L)
    p$min_abs_side_softclip   <- min(p$min_abs_side_softclip,   10L)
    p$min_coverage_drop       <- min(p$min_coverage_drop, 1.0)
    p$ptd_allow_asymmetric    <- TRUE
    p$min_itd_read_coverage   <- min(p$min_itd_read_coverage, 10)
    if (is.null(use_exon_graph)) p$use_exon_graph <- TRUE
  }
  
  if (is.null(config$exons)) config$exons <- config$target_exons
  
  detect_itd(
    bam_path = bam_path, gene_config = config,
    min_support = p$min_support, min_size = p$min_size,
    max_correction = p$max_correction, plot = plot, sample_name = sample_name,
    filter_intronic = p$filter_intronic, ptd_mode = p$ptd_mode,
    max_missing_kmers = p$max_missing_kmers,
    use_cigar_bp = p$use_cigar_bp, refine_bp = p$refine_bp,
    min_mapq = p$min_mapq, min_wt_reads = p$min_wt_reads,
    cluster_tolerance = p$cluster_tolerance, vaf_threshold = p$vaf_threshold,
    search_window = p$search_window,
    min_strand_bias = p$min_strand_bias, max_strand_bias = p$max_strand_bias,
    min_mean_support_mapq = p$min_mean_support_mapq, max_breakpoint_spread = p$max_breakpoint_spread,
    min_softclip_fraction = p$min_softclip_fraction, min_unique_breakpoints = p$min_unique_breakpoints,
    min_coverage_drop = p$min_coverage_drop, min_local_coverage = p$min_local_coverage, coverage_window = p$coverage_window,
    min_span_coverage = p$min_span_coverage, min_span_depth_fold_change = p$min_span_depth_fold_change,
    max_span_depth_fold_change = p$max_span_depth_fold_change,
    min_microhomology = p$min_microhomology, min_discordant_ratio = p$min_discordant_ratio,
    min_entropy = p$min_entropy, min_alignment_score = p$min_alignment_score,
    detect_orientation = p$detect_orientation,
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
    min_abs_side_softclip = p$min_abs_side_softclip,
    max_itd_length = p$max_itd_length,
    convert_long_to_ptd = p$convert_long_to_ptd,
    min_length = p$min_length,   
    max_length = p$max_length,
    use_kmers = use_kmers,
    exon_padding = exon_padding,
    prefer_extended_length = p$prefer_extended_length,
    ptd_allow_asymmetric = p$ptd_allow_asymmetric,
    ptd_use_local_assembly = p$ptd_use_local_assembly,
    merge_ptd_intervals = p$merge_ptd_intervals,
    merge_gap = p$merge_gap,
    min_ptd_length = p$min_ptd_length,
    max_plausible_fragments = p$max_plausible_fragments,
    use_exon_graph = p$use_exon_graph,
    verbose = verbose,
    debug = debug,
    ...
  )
}
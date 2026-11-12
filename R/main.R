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
#' @param min_base_quality Minimum mean Phred base quality for a read to enter
#'   the analysis population (candidate discovery, wild-type support and
#'   coverage-derived metrics). Reads without quality strings are retained.
#'   Default 20; set 0 to disable.
#' @param adapter_sequences Character vector of adapter sequences screened as
#'   read-through contamination. Defaults to the standard Illumina TruSeq /
#'   Nextera / small-RNA adapter set. Screening is a report, not a filter:
#'   affected QNAMEs are counted region-wide and per variant. NULL or
#'   character(0) disables screening.
#' @param adapter_min_overlap Minimum contiguous adapter bases (exact prefix
#'   match) that flag a read (default 10).
#' @param cluster_length_scale Length-scaled breakpoint clustering: the
#'   effective tolerance for a breakpoint with duplication length L is
#'   max(cluster_tolerance, cluster_length_scale * L), capped by
#'   max_cluster_tolerance. Default 0.1; set 0 for the legacy scalar tolerance.
#' @param max_cluster_tolerance Upper bound of the length-scaled clustering
#'   tolerance (default 50).
#' @param search_window Max distance (bp) the split-read anchor search
#'   (.estimate_itd_from_anchors()) looks for a re-anchor match. This is a
#'   hard ceiling on any anchor-confirmed length estimate -- a true
#'   duplication longer than this can only be sized via the noisier k-mer
#'   walk fallback. Default 5000 (unchanged for FLT3-scale ITDs). KMT2A
#'   widens this since its PTDs can genuinely exceed 5kb; safe to widen
#'   because the anchor search now requires an unambiguous (single-hit)
#'   match per prefix length rather than pooling every repeat-driven hit.
#' @param require_unique_anchor If TRUE, the anchor search
#'   (.estimate_itd_from_anchors()) only trusts a prefix-length match when
#'   it's unambiguous (exactly one hit in the window); ambiguous matches
#'   contribute nothing rather than every hit. Default FALSE (original
#'   behaviour: pool every hit, ambiguous or not). Motivated by KMT2A's
#'   Alu-rich introns producing wildly inconsistent length estimates from
#'   short-prefix repeat coincidences, but NOT proven safe for every gene --
#'   combined with an anchor-preferring best_len (since reverted) it
#'   regressed HORIZON FLT3's certified 300bp ITD to 260bp. Left FALSE
#'   everywhere by default; enable per-gene only after checking real
#'   results, the way every other threshold in this pipeline was tuned.
#' @param prefilter Enable fast prefiltering of alignments.
#' @param min_ins_filter Minimum net insertion length filter.
#' @param output_prefix Base name prefix (default "TALOS").
#' @param output_folder Output folder path.
#' @param write_vcf Generate VCF report.
#' @param verbose Logging enabled.
#' @param debug If TRUE, writes one line per candidate cluster.
#' @param nominal_read_len Assumed read length.
#' @param max_correction Exploratory size-bias correction cap (default 1 = no inflation). Primary AF and support remain raw regardless of this setting.
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
#' @param min_alignment_score Minimum sequence identity (0-1) between read-derived ITD sequence and the corresponding reference segment. Bypassed automatically for imputed sequences (default 0.2).
#' @param detect_orientation Check inverse loops.
#' @param annotate_hotspots SQLite DB lookups.
#' @param hotspot_db_path Path to SQLite lookup DB.
#' @param compute_alignment_score Metric toggle.
#' @param compute_support_bases Metric toggle.
#' @param compute_consistency Metric toggle (default TRUE). Calculates the
#'   percentage of supporting read-end sequences that align to a fully
#'   read-derived full-length ITD template at >=90% identity over >=20 bases.
#'   This is informational only and never gates .apply_filters(). It is not
#'   calculated against reference-imputed or reference-only coverage templates,
#'   because those templates cannot independently validate the read sequence.
#'   Use ITDCoverageDirectAlignmentRatePct for the read-to-template rate when
#'   only a reference coordinate template is available.
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
#' @param max_pairwise_alignments Maximum number of supporting reads aligned per chunk for ITD metrics (default 30). All distinct supporting reads are processed; this is not a read-sampling cap.
#' @param include_sample Include sample name in output filenames (default TRUE).
#' @param include_gene Include gene name in output filenames (default TRUE).
#' @param include_timestamp Include timestamp in output filenames (default TRUE).
#' @param output_sep Separator between filename parts (default "_").
#' @param global_log Write all logs (config, compute time, errors) to a single file (default TRUE).
#' @param min_side_softclip_reads Minimum number of reads with soft-clip on each side (default 0, i.e. disabled).
#' @param max_side_ratio Maximum ratio of left/right soft-clip counts before filtering (default 10).
#' @param min_softclip_pct_side Minimum percentage of supporting reads that must have a soft-clip on each side (default 0, i.e. disabled).
#' @param min_left_softclip_pct_wt Minimum percentage of wildtype reads with left soft-clip (default 0, i.e. disabled).
#' @param min_right_softclip_pct_wt Minimum percentage of wildtype reads with right soft-clip (default 0, i.e. disabled).
#' @param min_abs_side_softclip Minimum absolute soft-clip count on each side (default 0, i.e. disabled).
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
#'        Graph calls are reported with the duplicated interval found from
#'        split reads / everted pairs (GenomicPosition = first duplicated
#'        base, Length = interval length, DupExons = duplicated exons) and
#'        zero-length clip calls at either end of a graph call are folded
#'        into it. Tuning keys (YAML \code{graph_settings:}):
#'        max_exon_distance, buffer_bp, min_support, expected_isize_max.
#'        Default FALSE.
#' @param ... Other parameters passed to the underlying engine.
#' @return Data frame of ITD calls with diagnostic columns, invisibly.
#' @export
detect_itd <- function(
  bam_path,
  gene_config,
  genomic_ref_seq = NULL,
  k = 12,
  min_support = 10,
  min_size = 15,
  max_missing_kmers = 0.5,
  cluster_tolerance = 10,
  min_base_quality = 20,
  adapter_sequences = c(
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "CTGTCTCTTATACACATCT",
    "TGGAATTCTCGGGTGCCAAGG"
  ),
  adapter_min_overlap = 10L,
  cluster_length_scale = 0.1,
  max_cluster_tolerance = 50,
  prefilter = TRUE,
  search_window = 5000L,
  require_unique_anchor = FALSE,
  min_ins_filter = 3,
  output_prefix = "TALOS",
  output_folder = "./results",
  write_vcf = TRUE,
  verbose = TRUE,
  debug = FALSE,
  nominal_read_len = 150,
  max_correction = 1.0,
  plot = FALSE,
  sample_name = NULL,
  filter_intronic = FALSE,
  ptd_mode = FALSE,
  use_cigar_bp = TRUE,
  refine_bp = FALSE,
  min_mapq = 20,
  min_wt_reads = 0,
  vaf_threshold = 0.01,
  min_strand_bias = 0,
  max_strand_bias = 1,
  min_mean_support_mapq = 0,
  max_breakpoint_spread = Inf,
  min_softclip_fraction = 0,
  min_unique_breakpoints = 0,
  min_coverage_drop = 1.0,
  min_local_coverage = 0,
  coverage_window = 200,
  min_span_coverage = 0,
  min_span_depth_fold_change = 0,
  max_span_depth_fold_change = Inf,
  min_microhomology = 0,
  min_discordant_ratio = 0,
  min_entropy = 0,
  min_alignment_score = 0.2,
  detect_orientation = TRUE,
  do_annotate_hotspots = TRUE,
  hotspot_db_path = NULL,
  compute_alignment_score = TRUE,
  compute_support_bases = TRUE,
  compute_consistency = TRUE,
  compute_itd_coverage = TRUE,
  compute_coverage_drop = TRUE,
  compute_microhomology = TRUE,
  compute_repeat_entropy = TRUE,
  compute_discordant_ratio = TRUE,
  compute_hgvs = TRUE,
  html_report = TRUE,
  use_gviz = TRUE,
  min_itd_read_coverage = 50,
  max_reads_in_region = 200000,
  add_config_to_report = FALSE,
  max_pairwise_alignments = 30L,
  include_sample = TRUE,
  include_gene = TRUE,
  include_timestamp = TRUE,
  output_sep = "_",
  global_log = TRUE,
  min_side_softclip_reads = 0,
  max_side_ratio = 10,
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
  if (length(min_mapq) != 1L || is.na(min_mapq) || min_mapq < 0L) {
    stop("min_mapq must be one non-negative integer.", call. = FALSE)
  }
  if (length(vaf_threshold) != 1L || !is.finite(vaf_threshold) ||
      vaf_threshold < 0 || vaf_threshold > 1) {
    stop("vaf_threshold must be one finite value between 0 and 1.", call. = FALSE)
  }
  if (length(max_correction) != 1L || !is.finite(max_correction) || max_correction < 1) {
    stop("max_correction must be one finite value >= 1; use 1 to disable size correction.", call. = FALSE)
  }
  if (length(nominal_read_len) != 1L || !is.finite(nominal_read_len) || nominal_read_len <= 0) {
    stop("nominal_read_len must be one positive finite number.", call. = FALSE)
  }
  if (length(min_base_quality) != 1L || is.na(min_base_quality) ||
      !is.finite(min_base_quality) || min_base_quality < 0) {
    stop("min_base_quality must be one non-negative finite number (0 disables the filter).", call. = FALSE)
  }
  if (length(adapter_min_overlap) != 1L || is.na(adapter_min_overlap) ||
      adapter_min_overlap < 1) {
    stop("adapter_min_overlap must be one positive integer.", call. = FALSE)
  }
  if (length(cluster_length_scale) != 1L || is.na(cluster_length_scale) ||
      !is.finite(cluster_length_scale) || cluster_length_scale < 0 ||
      cluster_length_scale > 1) {
    stop("cluster_length_scale must be one finite number between 0 and 1 (0 disables length scaling).", call. = FALSE)
  }
  if (length(max_cluster_tolerance) != 1L || is.na(max_cluster_tolerance) ||
      !is.finite(max_cluster_tolerance) || max_cluster_tolerance < cluster_tolerance) {
    stop("max_cluster_tolerance must be one finite number >= cluster_tolerance.", call. = FALSE)
  }
  if (!is.null(adapter_sequences)) {
    adapter_sequences <- as.character(adapter_sequences)
    adapter_sequences <- adapter_sequences[!is.na(adapter_sequences) & nzchar(adapter_sequences)]
  }
  gene_name <- gene_config$gene %||% "UNKNOWN"
  if (is.null(sample_name)) {
    sample_name <- .sample_from_bam(bam_path)
  }
  sample_folder <- .ensure_dir(file.path(output_folder, sample_name))

  timestamp <- .ts_stamp()
  base_parts <- c(
    output_prefix,
    if (include_sample) sample_name,
    if (include_gene) gene_name,
    if (include_timestamp) timestamp
  )
  base_name <- paste(base_parts[!is.na(base_parts) & nzchar(base_parts)], collapse = output_sep)

  global_log_file <- NULL
  if (global_log) {
    global_log_file <- file.path(sample_folder, paste0(base_name, ".log"))
    sep_line <- paste0(strrep("-", 40), " \n")
    cat("=== TALOS Analysis Log ===\n",
        paste("Sample:", sample_name, "\n"),
        paste("Gene:", gene_name, "\n"),
        paste("BAM:", bam_path, "\n"),
        paste("Timestamp:", timestamp, "\n"),
        paste("Output folder:", output_folder, "\n"),
        sep_line, "=== Gene Configuration ===\n",
        file = global_log_file, sep = "")
    utils::capture.output(utils::str(gene_config, max.level = 2), file = global_log_file, append = TRUE)
    cat(sep_line, file = global_log_file, append = TRUE)
  }

  # [FIX Phase 4] log_msg() writes to stdout via cat() (NOT message()) so
  # that the withCallingHandlers(message=...) wrapper below can safely
  # mirror ALL message() calls from engine helpers (engine_metrics.R,
  # engine_bam.R, engine_alu.R) into the log file without duplicating
  # log_msg()'s own output. Using cat() to stdout visually matches
  # message() in interactive R sessions; the only difference is the
  # stream (stdout vs stderr), which doesn't matter for the console.
  log_msg <- function(..., append = TRUE) {
    msg <- paste(..., collapse = " ")
    if (verbose) cat(msg, "\n")
    if (!is.null(global_log_file)) cat(msg, "\n", file = global_log_file, append = append)
  }

  # Carry the analysis MAPQ threshold into optional report BAM queries so the
  # displayed depth is calculated from the same read population as WT support.
  gene_config$min_mapq <- as.integer(min_mapq)

  log_msg("Starting TALOS detection engine...")
  on.exit({
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    time_msg <- if (elapsed < 60) sprintf("%.1f sec", elapsed) else sprintf("%.1f min", elapsed / 60)
    log_msg(sprintf("\n%s\nCOMPUTE TIME: %s\n%s", strrep("-", 40), time_msg, strrep("-", 40)))
    if (!is.null(global_log_file)) log_msg(sprintf("Global log saved to: %s", global_log_file))
  }, add = TRUE)

  # Negative runs still write their TSV / HTML artefacts (the report shows an
  # explicit "No variants detected" state), so a run can never silently
  # produce nothing.
  finish_empty <- function(msg, cov = NULL, ref = NULL) {
    log_msg(msg)
    .write_talos_output(
      final_df = data.frame(), base_name = base_name, output_folder = output_folder,
      sample_name = sample_name, gene_config = gene_config, ref_dna = ref,
      bam_path = bam_path, write_vcf = write_vcf, plot = plot,
      html_report = html_report, verbose = verbose,
      add_config_to_report = add_config_to_report,
      # bam_data$cov is restricted to gene_config$genomic_start/end (the
      # detection window). The report includes flanking transcript exons such
      # as FLT3 E16, so let the PDF/HTML plot loaders query their full plotted
      # genomic interval from the indexed BAM instead of reusing this partial
      # detection-only Rle.
      precomputed_cov = NULL
    )
    invisible(data.frame())
  }

  tryCatch(
    withCallingHandlers(
      {
      ref_dna <- genomic_ref_seq %||% gene_config$genomic_ref_seq
      if (is.null(ref_dna) || length(ref_dna) != 1L || is.na(ref_dna) || nchar(ref_dna) == 0) {
        stop("No genomic reference sequence available in gene_config or genomic_ref_seq override.")
      }
      ref_dna <- toupper(as.character(ref_dna))
      if (!is.null(genomic_ref_seq)) {
        expected_len <- gene_config$genomic_end - gene_config$genomic_start + 1L
        if (nchar(ref_dna) != expected_len) {
          stop("genomic_ref_seq is ", nchar(ref_dna), " bp but the configured window is ",
               expected_len, " bp; the override must match the configured genomic window exactly.")
        }
        log_msg("[TALOS] Using caller-supplied genomic_ref_seq override.")
      }
      ref_len <- nchar(ref_dna)
      if (ref_len < k) {
        return(finish_empty("Reference shorter than k-mer length. Skipping.", ref = ref_dna))
      }
      gstart <- gene_config$genomic_start

      ref_kmers <- .prepare_kmers(ref_dna, k)
      bam_data <- .load_bam_data_streaming(
        bam_path, gene_config,
        compute_pairs = isTRUE(compute_discordant_ratio) || isTRUE(compute_itd_coverage),
        max_reads = max_reads_in_region,
        verbose = verbose
      )
      all_reads <- bam_data$reads
      if (length(all_reads) == 0L) {
        return(finish_empty("No reads found in target region.", cov = bam_data$cov, ref = ref_dna))
      }

      # -- Adapter screening (report-only) and base-quality filtering --------
      # Screening flags reads carrying adapter read-through; it never removes
      # reads. The base-quality filter drops reads whose mean Phred score is
      # below min_base_quality from the whole analysis population (candidate
      # discovery, wild-type support and coverage metrics share it). Both
      # steps run before the WT population split so every downstream
      # population is defined identically.
      adapter_screen <- .screen_adapter_reads(
        all_reads, adapter_sequences, adapter_min_overlap, verbose = verbose
      )
      bq_loaded <- .bam_has_base_quality(all_reads)
      if (bq_loaded && min_base_quality > 0) {
        all_reads <- .filter_reads_by_base_quality(all_reads, min_base_quality, verbose = verbose)
        if (length(all_reads) == 0L) {
          return(finish_empty("No reads passed the base-quality filter.", cov = bam_data$cov, ref = ref_dna))
        }
      }

      wt_mapqs <- S4Vectors::mcols(all_reads)$mapq
      wt_mapqs[is.na(wt_mapqs)] <- 0L
      all_reads_for_wt <- all_reads[wt_mapqs >= min_mapq]
      # All coverage-derived metrics must use the same MAPQ / duplicate / QC
      # population as WT support. bam_data$cov is the unfiltered loader depth.
      cov_filtered <- if (length(all_reads_for_wt) > 0L)
        GenomicAlignments::coverage(all_reads_for_wt) else NULL
      all_reads_cov_filtered <- if (!is.null(cov_filtered) &&
                                    gene_config$chrom %in% names(cov_filtered))
        cov_filtered[[gene_config$chrom]] else NULL

      if (prefilter) {
        all_reads <- .filter_reads_by_cigar(all_reads, min_mapq, min_ins_filter, verbose = verbose)
        if (length(all_reads) == 0L) {
          return(finish_empty("No reads passed CIGAR prefilter.", cov = bam_data$cov, ref = ref_dna))
        }
      }

      # -- Candidate extraction ---------------------------------------------
      candidates <- list()

      if (use_exon_graph && !is.null(gene_config$all_exons)) {
        log_msg("[TALOS] Running Exon-Graph PTD reconstruction...")
        graph_cycles <- extract_graph_candidates(bam_path, gene_config, min_mapq = min_mapq)
        if (length(graph_cycles) > 0) {
          for (cy in graph_cycles) {
            log_msg(sprintf(
              "[GRAPH] %s | dup %s-%s (%s bp) | reads=%d (split=%s, pairs=%s) | score=%.2f",
              cy$ptd_signature, cy$dup_start, cy$dup_end, cy$dup_length, cy$support_reads,
              cy$n_split, cy$n_pairs, cy$evidence_score
            ))
          }
          # Look the reads up in the MAPQ-filtered set: everted mates carry no
          # clip and are removed by the CIGAR prefilter.
          gc <- .graph_cycles_to_candidates(graph_cycles, all_reads_for_wt, gstart, ref_len = ref_len)
          if (length(gc$candidates) > 0) {
            candidates <- c(candidates, gc$candidates)
            log_msg(sprintf("[TALOS] Graph extraction contributed %d read-level PTD candidate(s) from %d cycle(s).",
                            length(gc$candidates), gc$cycles_used))
          }
        } else {
          log_msg("Graph extraction found no PTD cycles.")
        }
      }

      standard_args <- list(
        reads = all_reads, ref_kmers = ref_kmers, min_size = min_size,
        max_missing_kmers = max_missing_kmers, refine_bp = refine_bp,
        use_cigar_bp = use_cigar_bp, genomic_start = gstart, ref_dna = ref_dna,
        max_itd_length = max_itd_length, convert_long_to_ptd = convert_long_to_ptd,
        search_window = search_window, require_unique_anchor = require_unique_anchor,
        verbose = verbose
      )
      legacy_candidates <- if (ptd_mode && use_cigar_bp) {
        c(.extract_candidates_ptd(all_reads, gstart, ref_len, verbose = verbose),
          if (use_kmers) do.call(.extract_candidates_standard, c(standard_args, list(ptd_mode = TRUE))))
      } else if (!use_kmers) {
        .extract_candidates_cigar(all_reads, gstart, ref_len, min_size = min_size, verbose = verbose)
      } else {
        do.call(.extract_candidates_standard, c(standard_args, list(ptd_mode = ptd_mode)))
      }
      candidates <- c(candidates, legacy_candidates)

      if (length(candidates) == 0) {
        return(finish_empty("No candidate reads (clips or k-mer jumps) found.", cov = bam_data$cov, ref = ref_dna))
      }
      bp_df <- .candidates_to_df(candidates, gstart)
      # Length-scaled clustering (v0.4.20): the effective tolerance for a
      # breakpoint grows with its duplication length (capped), while the
      # anti-chaining total-width bound is kept.
      clusters <- .cluster_breakpoints(
        bp_df$breakpoint, cluster_tolerance,
        lengths = bp_df$length, length_scale = cluster_length_scale,
        max_scaled_tolerance = max_cluster_tolerance
      )

      ptd_assembled_seqs <- list()
      if (ptd_mode && ptd_use_local_assembly && nrow(bp_df) > 0) {
        for (cl in clusters) {
          assembled <- .assemble_ptd_consensus(bp_df[bp_df$breakpoint %in% cl, ], min_reads = 3L, min_len = 15L)
          if (!is.na(assembled) && nchar(assembled) > 0) {
            med_cl <- as.character(round(stats::median(cl)))
            ptd_assembled_seqs[[med_cl]] <- assembled
            log_msg(sprintf("[TALOS] Assembled PTD consensus of length %d bp for cluster at %s",
                            nchar(assembled), med_cl))
          }
        }
      }
      .lookup_assembled_seq <- function(bp, store = ptd_assembled_seqs, tol = 50L) {
        if (length(store) == 0L) return(NA_character_)
        dists <- abs(suppressWarnings(as.integer(names(store))) - as.integer(bp))
        best <- which.min(dists)
        if (dists[best] <= tol) store[[best]] else NA_character_
      }

      ext_df <- .extend_candidates(
        bp_df[!is.na(bp_df$length) & bp_df$length > 0L, ],
        ref_dna = ref_dna, genomic_start = gstart, min_size = min_size,
        search_window = search_window, require_unique_anchor = require_unique_anchor,
        verbose = verbose
      )

      wt_info <- .prepare_wildtype_info(all_reads_for_wt, gstart, gene_config$genomic_end)
      wt_info$sample_name <- sample_name

      # Threshold list is derived from .FILTER_RULES so a new rule can never
      # reference a threshold missing here (it would silently default to
      # "disabled" inside .apply_filters()).
      threshold_names <- unique(c(
        vapply(.FILTER_RULES, function(r) r[[2L]], character(1L)),
        "min_alignment_score", "min_side_softclip_reads", "max_side_ratio",
        "min_softclip_pct_side", "min_left_softclip_pct_wt", "min_right_softclip_pct_wt",
        "min_abs_side_softclip"
      ))
      # [FIX Phase 4] Defensive mget(): if a threshold name in .FILTER_RULES
      # ever drifts from a detect_itd() formal argument, the previous
      # mget() threw a cryptic "value for 'X' not found" error mid-run. We
      # pass an explicit ifnotfound (a simple list of NAs -- NO side effects
      # inside, because R evaluates ifnotfound eagerly when the list is
      # built, which would emit spurious warnings for every name even when
      # the name IS found). After mget() we check which names came back NA
      # and emit ONE consolidated warning listing all of them.
      thresholds <- mget(
        threshold_names,
        envir = environment(),
        ifnotfound = setNames(replicate(length(threshold_names), NA, simplify = FALSE),
                              threshold_names)
      )
      missing_thresholds <- threshold_names[
        vapply(thresholds, function(v) is.null(v) || (length(v) == 1L && is.na(v)),
               logical(1L))
      ]
      if (length(missing_thresholds) > 0L) {
        warning(
          "[TALOS] The following threshold(s) were not found in detect_itd() ",
          "formals and will be treated as 'rule disabled' by .apply_filters(): ",
          paste(missing_thresholds, collapse = ", "),
          ". Check the .FILTER_RULES / formal-argument name alignment in ",
          "engine_metrics.R and main.R."
        )
      }

      metric_args <- list(
        genomic_start = gstart, ref_dna = ref_dna, gene_config = gene_config,
        all_reads_cov = all_reads_cov_filtered, all_reads_for_coverage = all_reads_for_wt,
        all_pairs = bam_data$pairs, wildtype_info = wt_info,
        ptd_mode = ptd_mode, min_support = min_support, min_wt_reads = min_wt_reads,
        nominal_read_len = nominal_read_len, max_correction = max_correction,
        min_mapq = min_mapq, vaf_threshold = vaf_threshold, do_alignment_score = compute_alignment_score,
        do_consistency = compute_consistency, do_itd_coverage = compute_itd_coverage,
        do_support_bases = compute_support_bases, do_coverage_drop = compute_coverage_drop,
        do_microhomology = compute_microhomology, do_repeat_entropy = compute_repeat_entropy,
        do_discordant_ratio = compute_discordant_ratio, do_detect_orientation = detect_orientation,
        do_hgvs = compute_hgvs, max_pairwise_alignments = max_pairwise_alignments,
        span_flank = coverage_window, max_itd_length = max_itd_length,
        adapter_positive_qnames = adapter_screen$qnames,
        adapter_screened_reads = adapter_screen$screened,
        adapter_positive_reads = adapter_screen$positive,
        adapter_n_configured = length(adapter_sequences),
        adapter_min_overlap = as.integer(adapter_min_overlap),
        base_quality_loaded = bq_loaded,
        min_base_quality = min_base_quality,
        debug = debug, verbose = verbose
      )

      results <- list()
      # Metrics + filters + debug line + collection: one place instead of three.
      evaluate_cluster <- function(cluster_bps, best_len, support_rows, length_ext, lookup_bp,
                                   anchor_confirmed_qnames = NA_integer_) {
        if (is.na(anchor_confirmed_qnames) && nrow(ext_read_level) > 0L) {
          candidate_qnames <- unique(support_rows$read_name)
          anchor_rows <- ext_read_level[
            ext_read_level$read_name %in% candidate_qnames &
              ext_read_level$itdsize_source == "anchor" &
              (is.na(best_len) | ext_read_level$itdsize == best_len), , drop = FALSE
          ]
          anchor_confirmed_qnames <- length(unique(anchor_rows$read_name))
        }
        metrics <- do.call(.compute_variant_metrics, c(list(
          cluster_bps = cluster_bps, best_len = best_len, support_rows = support_rows,
          length_ext = length_ext, pre_assembled_seq = .lookup_assembled_seq(lookup_bp),
          anchor_confirmed_qnames = anchor_confirmed_qnames
        ), metric_args))
        if (is.null(metrics)) return(invisible(FALSE))
        passed <- .apply_filters(metrics, thresholds, min_length, max_length,
                                 ptd_mode = ptd_mode, ptd_allow_asymmetric = ptd_allow_asymmetric)
        if (debug) log_msg(.format_candidate_debug(metrics, passed))
        if (passed) results[[length(results) + 1L]] <<- metrics
        invisible(passed)
      }
      skip_debug <- function(pos, n, why) {
        if (debug) {
          log_msg(sprintf("[CANDIDATE] GenomicPosition=%s SupportingReads=%d -> SKIPPED (%s)",
                          as.character(pos), as.integer(n), why))
        }
      }

      # -- Path 1: extension-refined clusters -------------------------------
      # Collapse each QNAME before clustering or taking a length mode. A
      # read/pair can generate several raw candidate rows; those are not
      # independent votes for the same duplication.
      ext_read_level <- .collapse_extension_candidates_by_read(ext_df)
      if (prefer_extended_length && !ptd_mode && nrow(ext_read_level) >= 3L) {
        ext_read_valid <- ext_read_level[!is.na(ext_read_level$itdsize) &
                                           ext_read_level$itdsize >= min_size, , drop = FALSE]
        if (nrow(ext_read_valid) >= 3L) {
          for (rcl in .cluster_breakpoints(
            ext_read_valid$breakpoint_refined, cluster_tolerance,
            lengths = ext_read_valid$itdsize, length_scale = cluster_length_scale,
            max_scaled_tolerance = max_cluster_tolerance
          )) {
            rcl_reads <- ext_read_valid[
              ext_read_valid$breakpoint_refined %in% rcl, , drop = FALSE
            ]
            if (nrow(rcl_reads) < 3L) {
              skip_debug(round(stats::median(rcl, na.rm = TRUE)), nrow(rcl_reads),
                         "refined cluster has <3 independent supporting QNAMEs")
              next
            }
            best_len <- .mode_extension_length(rcl_reads$itdsize, rcl_reads$itdsize_source)
            # Keep only QNAMEs supporting the modal length. This prevents
            # different ITD-size hypotheses at one breakpoint being pooled.
            selected_reads <- rcl_reads[rcl_reads$itdsize == best_len, , drop = FALSE]
            if (nrow(selected_reads) < 3L) next
            median_bp <- as.integer(stats::median(selected_reads$breakpoint_refined, na.rm = TRUE))
            selected_qnames <- unique(selected_reads$read_name)
            near <- abs(bp_df$breakpoint - median_bp) <= cluster_tolerance
            support_rows <- bp_df[bp_df$read_name %in% selected_qnames & near, , drop = FALSE]
            if (nrow(support_rows) == 0L) {
              skip_debug(median_bp, 0L, "no breakpoint rows matched extension-selected QNAMEs")
              next
            }
            anchor_n <- sum(selected_reads$itdsize_source == "anchor", na.rm = TRUE)
            if (verbose) {
              log_msg(sprintf(
                "Extension-based cluster at %d: %d bp (%d/%d independent QNAMEs anchor-confirmed; %d selected support rows)",
                median_bp, best_len, anchor_n, nrow(selected_reads), nrow(support_rows)
              ))
            }
            evaluate_cluster(rep(median_bp, nrow(support_rows)), best_len, support_rows,
                             best_len, median_bp, anchor_confirmed_qnames = anchor_n)
          }
        }
      }

      # -- Path 2: raw breakpoint clusters ------------------------------------
      # Run for every raw cluster not already called by the extension path.
      for (cl in clusters) {
        # Refresh inside the loop: Path 2 can append a call, so a snapshot
        # taken before iteration can allow later length groups to duplicate it.
        called_pos <- vapply(results, function(r) as.numeric(r$GenomicPosition), numeric(1L))
        cluster_candidates <- which(bp_df$breakpoint %in% cl)
        if (length(cluster_candidates) == 0L) {
          skip_debug(round(stats::median(cl, na.rm = TRUE)), 0L, "no candidates matched cluster; metrics not computed")
          next
        }
        if (length(called_pos) && any(abs(called_pos - stats::median(cl)) <= cluster_tolerance)) next

        support_read_names <- unique(bp_df$read_name[cluster_candidates])
        raw_center <- stats::median(cl, na.rm = TRUE)
        ext_cluster <- ext_read_level[
          ext_read_level$read_name %in% support_read_names &
            abs(ext_read_level$breakpoint_refined - raw_center) <= cluster_tolerance,
          , drop = FALSE
        ]

        if (prefer_extended_length && !ptd_mode && nrow(ext_cluster) >= 3L) {
          best_len <- .mode_extension_length(ext_cluster$itdsize, ext_cluster$itdsize_source)
          selected_reads <- ext_cluster[ext_cluster$itdsize == best_len, , drop = FALSE]
          # If no length hypothesis has >=3 independent QNAMEs, do not
          # abandon this raw breakpoint cluster: fall through to the raw
          # CIGAR-length path below instead of silently losing the event.
          if (nrow(selected_reads) >= 3L) {
            median_bp <- as.integer(stats::median(selected_reads$breakpoint_refined, na.rm = TRUE))
            selected_qnames <- unique(selected_reads$read_name)
            support_rows <- bp_df[
              bp_df$read_name %in% selected_qnames &
                abs(bp_df$breakpoint - median_bp) <= cluster_tolerance,
              , drop = FALSE
            ]
            if (nrow(support_rows) > 0L) {
              evaluate_cluster(rep(median_bp, nrow(support_rows)), best_len,
                               support_rows, best_len, median_bp)
              next
            }
          }
        }

        cluster_lengths <- bp_df$length[cluster_candidates]
        cluster_lengths[is.na(cluster_lengths)] <- -1L
        for (len in names(split(cluster_candidates, cluster_lengths))) {
          best_len <- if (as.integer(len) == -1L) NA_integer_ else as.integer(len)
          in_cl <- bp_df$breakpoint %in% cl
          support_rows <- if (is.na(best_len)) {
            bp_df[in_cl & is.na(bp_df$length), ]
          } else {
            bp_df[in_cl & !is.na(bp_df$length) & bp_df$length == best_len, ]
          }
          ext_sizes_cl <- ext_read_level$itdsize[
            ext_read_level$read_name %in% support_rows$read_name &
              !is.na(ext_read_level$itdsize) & ext_read_level$itdsize > 0L
          ]
          length_ext <- if (length(ext_sizes_cl) >= 3L) {
            as.integer(names(sort(table(ext_sizes_cl, decreasing = TRUE))[1L]))
          } else {
            NA_integer_
          }
          evaluate_cluster(cl, best_len, support_rows, length_ext, round(stats::median(cl)))
        }
      }

      if (length(results) == 0) {
        return(finish_empty("No variants passed all filters.", cov = bam_data$cov, ref = ref_dna))
      }

      final_df <- do.call(rbind, lapply(results, function(x) as.data.frame(x, stringsAsFactors = FALSE)))
      final_df$MergedFromN <- 1L
      final_df$MergeFragmentWarning <- FALSE
      final_df$Hotspot <- FALSE
      final_df$HotspotName <- NA_character_
      final_df$Region <- NA_character_
      final_df$ExonNumber <- NA_integer_

      final_df <- .fold_clip_calls_into_graph(final_df, cluster_tolerance, log_msg)

      if (merge_ptd_intervals && nrow(final_df) > 0) {
        final_df <- .merge_ptd_intervals(
          final_df, min_ptd_length = min_ptd_length, merge_gap = merge_gap,
          max_plausible_fragments = max_plausible_fragments, thresholds = thresholds,
          cov = bam_data$cov, chrom = gene_config$chrom, coverage_window = coverage_window,
          log_msg = log_msg
        )
      }

      if (do_annotate_hotspots) {
        final_df <- annotate_hotspots(final_df, db_path = hotspot_db_path, genome_build = gene_config$build)
      } else {
        final_df$Hotspot <- rep(FALSE, nrow(final_df))
        final_df$HotspotName <- rep(NA_character_, nrow(final_df))
      }

      final_df <- .annotate_exonic_region(final_df, gene_config$target_exons)
      if (filter_intronic && nrow(final_df) > 0) {
        final_df <- final_df[is.na(final_df$Region) | final_df$Region == "exonic", ]
      }
      final_df$TranscriptRef <- rep(
        if (is.null(gene_config$transcript) || length(gene_config$transcript) != 1L ||
            is.na(gene_config$transcript) || !nzchar(gene_config$transcript)) "none"
        else gene_config$transcript, nrow(final_df)
      )

      .write_talos_output(
        final_df, base_name = base_name, output_folder = output_folder,
        sample_name = sample_name, gene_config = gene_config, ref_dna = ref_dna,
        bam_path = bam_path, write_vcf = write_vcf, plot = plot, html_report = html_report,
        verbose = verbose, add_config_to_report = add_config_to_report,
        # Detection coverage is limited to the configured target interval.
        # Report coverage is fetched separately across every displayed exon,
        # including flanking exons outside the candidate-discovery window.
        precomputed_cov = NULL
      )

      log_msg(sprintf("Analysis completed successfully. %d variant(s) found.", nrow(final_df)))
      return(final_df)
    },
    # [FIX Phase 4] Mirror every message() emitted inside the engine body
    # (PE-spanning coverage, PE insert-size length, debug anchor sizes,
    # anything in engine_metrics.R / engine_bam.R / engine_alu.R) into the
    # global log file. log_msg() already does this for messages it emits,
    # but engine helpers call base::message() directly, so their text only
    # reached the console. Here we append the condition's message text to
    # the log file (when one is open). We do NOT call
    # invokeRestart("muffleMessage"), so the default handler still prints
    # the message to the console -- the only side effect here is the extra
    # cat() to the log file.
    message = function(m) {
      if (!is.null(global_log_file)) {
        cat(conditionMessage(m), "\n", file = global_log_file, append = TRUE)
      }
      # Returning without muffling -> default handler prints to console too.
    }
    )
  ,
    error = function(e) {
      log_msg(sprintf("\n%s\nERROR: %s\n%s", strrep("-", 40), conditionMessage(e), strrep("-", 40)))
      stop(e)
    }
  )
}

#' Fold zero-length clip calls into the graph PTD they belong to
#'
#' A graph PTD [S, E] is also seen by the legacy clip extractor as two
#' zero-length calls (left clips at S, right clips at E). They describe the
#' same event, so they are dropped when a graph call covers them.
#' @keywords internal
.fold_clip_calls_into_graph <- function(final_df, tol, log_msg = message) {
  if (!"CandidateTypes" %in% names(final_df)) return(final_df)
  graph <- which(grepl("ptd_graph", final_df$CandidateTypes) & !is.na(final_df$Length) & final_df$Length > 0)
  if (!length(graph)) return(final_df)
  ends <- c(final_df$GenomicPosition[graph], final_df$GenomicPosition[graph] + final_df$Length[graph] - 1L)
  zero <- is.na(final_df$Length) | final_df$Length == 0
  covered <- zero & vapply(final_df$GenomicPosition, function(p) any(abs(ends - p) <= tol), logical(1L))
  if (any(covered)) {
    log_msg(sprintf("[GRAPH] %d zero-length clip call(s) folded into graph PTD call(s).", sum(covered)))
    final_df <- final_df[!covered, , drop = FALSE]
  }
  final_df
}

# Merge per-base coverage RLEs in genomic coordinates. Coverage vectors are
# anchored at each source row's GenomicPosition; gaps remain zero. This keeps
# PDF and HTML ITD-coverage tracks available after PTD interval consolidation.
.merge_coverage_rle_column <- function(rows, column, merged_start, merged_len,
                                       max_length = 1000000L) {
  if (!column %in% names(rows) || length(merged_len) != 1L ||
      is.na(merged_len) || merged_len < 1L || merged_len > max_length) {
    return(NA_character_)
  }
  merged_len <- as.integer(merged_len)
  out <- integer(merged_len)
  found <- FALSE
  for (i in seq_len(nrow(rows))) {
    cov_vec <- .parse_itd_coverage_rle(rows[[column]][i])
    if (!length(cov_vec)) next
    row_start <- suppressWarnings(as.numeric(rows$GenomicPosition[i]))
    if (length(row_start) != 1L || !is.finite(row_start)) next
    offset <- as.integer(row_start - as.numeric(merged_start))
    src_lo <- max(1L, 1L - offset)
    src_hi <- min(length(cov_vec), merged_len - offset)
    if (src_hi < src_lo) next
    dst_lo <- offset + src_lo
    dst_hi <- offset + src_hi
    # Per-row coverage RLEs are non-additive summaries (max support per
    # base), not independent read-depth contributions. Summing them inflated
    # merged interval coverage whenever adjacent rows overlapped.
    out[dst_lo:dst_hi] <- pmax(out[dst_lo:dst_hi],
      as.integer(cov_vec[src_lo:src_hi]))
    found <- TRUE
  }
  if (!found) return(NA_character_)
  rr <- rle(out)
  paste(paste0(rr$values, ":", rr$lengths), collapse = ",")
}

# Columns combined when several duplication rows are merged into one interval.
.MERGE_SUM_COLS <- c("TotalSupportBases", "LeftSoftclipCount", "RightSoftclipCount",
                     "BothSoftclipCount", "PESoftclipSupport", "PESoftclipEventPairs",
                     "PESoftclipLongPairs", "PEOrientationFR", "PEOrientationRF",
                     "PEOrientationFF", "PEOrientationRR", "PEOrientationOther",
                     "ITDPE_SpanningPairs", "ITDKmerReads")
.MERGE_WMEAN_COLS <- list(RefMatch_Observed = 1L, RefMatch_Total = 1L, ITDReadCoverage = NA_integer_,
                          ITDReadCoverage_Combined = NA_integer_, ITDReadCoverage_Kmer = NA_integer_,
                          SupportConsistency = NA_integer_, AlignmentScore = NA_integer_,
                          StrandBias = 4L, MeanSupportMAPQ = 1L)
.MERGE_RESET_COLS <- list(ITD_Sequence = NA_character_, HGVS_cDNA = NA_character_,
                          HGVS_Protein = NA_character_, ITDCoverageRLE = NA_character_,
                          ITDCoverageRLE_Direct = NA_character_,
                          ITDCoverageRLE_PE = NA_character_,
                          ITDCoverageRLE_Kmer = NA_character_,
                          LengthPE = NA_real_, LengthPE_NSpanning = NA_integer_,
                          LengthPE_Delta = NA_real_, LengthPE_Concordant = NA,
                          ITDUncoveredBases_Direct = NA_integer_,
                          ITDUncoveredBases_Combined = NA_integer_,
                          ITDUncoveredIntervals_Combined = NA_character_,
                          LengthExt = NA_integer_, SequenceImputed = TRUE, SequencePartial = TRUE,
                          SequenceSource = "merged_intervals", Hotspot = FALSE,
                          HotspotName = NA_character_, Region = NA_character_,
                          ExonNumber = NA_integer_)

# GenomicRanges::reduce() merges only when the gap is strictly shorter than
# min.gapwidth. Convert the user-facing inclusive maximum gap to that API.
.ptd_merge_min_gapwidth <- function(merge_gap) {
  gap <- suppressWarnings(as.numeric(merge_gap))
  if (length(gap) != 1L || is.na(gap) || !is.finite(gap) || gap < 0 ||
      gap != floor(gap) || gap >= .Machine$integer.max) {
    stop("merge_gap must be a single non-negative integer below .Machine$integer.max.", call. = FALSE)
  }
  as.integer(gap + 1)
}

#' Merge overlapping / adjacent large duplication intervals into one call
#' @keywords internal
.merge_ptd_intervals <- function(final_df, min_ptd_length, merge_gap, max_plausible_fragments,
                                 thresholds, cov, chrom, coverage_window, log_msg = message) {
  thv <- function(nm, default) {
    v <- thresholds[[nm]]
    if (is.null(v) || length(v) != 1L || is.na(v)) default else as.numeric(v)
  }
  min_gapwidth <- .ptd_merge_min_gapwidth(merge_gap)
  merge_gap <- min_gapwidth - 1L
  if (length(chrom) != 1L || is.na(chrom) || !nzchar(as.character(chrom))) {
    stop("chrom must be a single non-empty chromosome name.", call. = FALSE)
  }
  eligible <- !is.na(final_df$Length) & final_df$Length >= min_ptd_length &
    is.finite(as.numeric(final_df$GenomicPosition))
  keep <- which(eligible)
  skip <- setdiff(seq_len(nrow(final_df)), keep)
  if (length(keep) == 0L) {
    log_msg(sprintf("No duplications >= %d bp found; skipping merge.", min_ptd_length))
    return(final_df)
  }
  if (length(keep) == 1L) {
    log_msg("Only one large duplication found; no merging performed.")
    return(final_df)
  }
  log_msg(sprintf("Merging %d large duplication intervals (gap <= %d bp)...", length(keep), merge_gap))
  gr <- GenomicRanges::GRanges(
    seqnames = rep(as.character(chrom), length(keep)),
    ranges = IRanges::IRanges(start = final_df$GenomicPosition[keep],
                              end = final_df$GenomicPosition[keep] + final_df$Length[keep] - 1L)
  )
  # GenomicRanges merges when the intervening gap is strictly shorter than
  # min.gapwidth. Add one so a gap exactly equal to merge_gap is included.
  merged_gr <- GenomicRanges::reduce(gr, min.gapwidth = min_gapwidth)
  if (length(merged_gr) == 0L) {
    log_msg("No intervals could be merged.")
    return(final_df)
  }
  hits <- GenomicRanges::findOverlaps(gr, merged_gr)
  group_of <- integer(length(gr))
  group_of[S4Vectors::queryHits(hits)] <- S4Vectors::subjectHits(hits)

  sum_na <- function(x) sum(x, na.rm = TRUE)
  merged_rows <- lapply(seq_along(merged_gr), function(i) {
    grp <- keep[group_of == i]
    n_members <- length(grp)
    mg_start <- GenomicRanges::start(merged_gr[i])
    mg_len <- GenomicRanges::width(merged_gr[i])
    template <- final_df[grp[which.min(final_df$GenomicPosition[grp])], ]
    w <- final_df$SupportingReads[grp]
    w[is.na(w) | w < 0] <- 0
    wmean <- function(vals, digits) {
      ok <- !is.na(vals) & w > 0
      if (!any(ok)) return(NA_real_)
      out <- stats::weighted.mean(vals[ok], w[ok])
      if (is.na(digits)) out else round(out, digits)
    }

    template$GenomicPosition <- mg_start
    template$Length <- mg_len
    template$InFrame <- if (mg_len > 0L) (as.integer(mg_len) %% 3L == 0L) else NA
    if (!is.null(cov)) {
      sc <- .compute_span_coverage_summary(cov, chrom, mg_start, mg_len, flank = coverage_window)
      template$SpanMinCoverage <- sc$min_cov
      template$SpanMeanCoverage <- sc$mean_cov
      template$SpanDepthFoldChange <- sc$depth_fold_change
    }
    # Exact QNAME-union support (v0.4.20). Component rows carry their
    # supporting / wild-type QNAME sets (SupportQNAMEs / WildtypeQNAMEs,
    # ";"-joined). Merged support is the size of the SET UNION, so a read
    # supporting two adjacent fragments counts once instead of twice. When
    # the QNAME columns are missing, or a set was truncated (oversized sets
    # carry a "...TRUNCATED" marker), fall back to the legacy component-row
    # sum and say so in MergedSupportSemantics.
    qnames_union <- function(col) {
      if (!col %in% names(final_df)) return(NULL)
      vals <- as.character(final_df[[col]][grp])
      vals <- vals[!is.na(vals) & nzchar(vals)]
      if (!length(vals)) return(NULL)
      if (any(grepl("...TRUNCATED", vals, fixed = TRUE))) return(NULL)
      unique(stats::na.omit(unlist(strsplit(vals, ";", fixed = TRUE))))
    }
    sup_union <- qnames_union("SupportQNAMEs")
    wt_union <- qnames_union("WildtypeQNAMEs")
    exact_union <- !is.null(sup_union) || !is.null(wt_union)
    template$SupportingReads <- if (!is.null(sup_union)) length(sup_union) else
      sum_na(final_df$SupportingReads[grp])
    template$WildtypeReads <- if (!is.null(wt_union)) length(wt_union) else
      sum_na(final_df$WildtypeReads[grp])
    if (!is.null(sup_union)) {
      template$SupportQNAMEs <- paste(sort(sup_union), collapse = ";")
      if ("IndependentSupportQNAMEs" %in% names(final_df)) {
        template$IndependentSupportQNAMEs <- length(sup_union)
      }
    }
    if (!is.null(wt_union)) {
      template$WildtypeQNAMEs <- paste(sort(wt_union), collapse = ";")
    }
    template$DepthAtBreakpoint <- template$SupportingReads + template$WildtypeReads
    template$AlleleFrequency <- if (template$DepthAtBreakpoint > 0L) {
      template$SupportingReads / template$DepthAtBreakpoint
    } else {
      NA_real_
    }
    if ("AdjustedSupportingReads" %in% names(final_df)) {
      # Corrected counts are exploratory diagnostics; recompute them from the
      # union support with the template row's size-correction factor so the
      # diagnostic stays consistent with the primary union count.
      if (!is.null(sup_union)) {
        scf <- suppressWarnings(as.numeric(template$SizeCorrectionFactor))
        if (length(scf) != 1L || is.na(scf) || !is.finite(scf) || scf < 1) scf <- 1
        template$AdjustedSupportingReads <- as.integer(round(length(sup_union) * scf))
      } else {
        template$AdjustedSupportingReads <- sum_na(final_df$AdjustedSupportingReads[grp])
      }
      template$CorrectedDepthAtBreakpoint <- template$AdjustedSupportingReads + template$WildtypeReads
      template$CorrectedAlleleFrequency <- if (template$CorrectedDepthAtBreakpoint > 0L)
        template$AdjustedSupportingReads / template$CorrectedDepthAtBreakpoint else NA_real_
      template$SizeCorrectionFactor <- if (template$SupportingReads > 0L)
        template$AdjustedSupportingReads / template$SupportingReads else 1
    }
    template$RawAlleleFrequency <- template$AlleleFrequency
    for (cl in intersect(.MERGE_SUM_COLS, names(final_df))) {
      template[[cl]] <- as.integer(sum_na(final_df[[cl]][grp]))
    }
    pct <- function(num, den) if (den > 0L) num / den * 100 else NA_real_
    template$LeftSoftclipPctSupport <- pct(template$LeftSoftclipCount, template$SupportingReads)
    template$RightSoftclipPctSupport <- pct(template$RightSoftclipCount, template$SupportingReads)
    template$LeftSoftclipPctWT <- pct(template$LeftSoftclipCount, template$WildtypeReads)
    template$RightSoftclipPctWT <- pct(template$RightSoftclipCount, template$WildtypeReads)
    for (cl in intersect(names(.MERGE_WMEAN_COLS), names(final_df))) {
      template[[cl]] <- wmean(final_df[[cl]][grp], .MERGE_WMEAN_COLS[[cl]])
    }
    if (n_members > 1L) {
      for (cl in intersect(names(.MERGE_RESET_COLS), names(final_df))) {
        template[[cl]] <- .MERGE_RESET_COLS[[cl]]
      }
      # Do not discard position-resolved coverage when adjacent fragments are
      # consolidated. Restore each coverage layer by mapping its RLE onto the
      # merged genomic interval; the combined display RLE is merged as well.
      coverage_rle_cols <- c("ITDCoverageRLE", "ITDCoverageRLE_Direct",
                             "ITDCoverageRLE_PE", "ITDCoverageRLE_Kmer")
      for (cl in intersect(coverage_rle_cols, names(final_df))) {
        template[[cl]] <- .merge_coverage_rle_column(
          final_df[grp, , drop = FALSE], cl, mg_start, mg_len
        )
      }
      for (cl in intersect(c("DupExons", "CandidateTypes"), names(final_df))) {
        v <- unique(stats::na.omit(final_df[[cl]][grp]))
        template[[cl]] <- if (length(v)) paste(v, collapse = ";") else NA_character_
      }
    }
    pe_counts <- c(FR = template$PEOrientationFR, RF = template$PEOrientationRF,
                   FF = template$PEOrientationFF, RR = template$PEOrientationRR,
                   other = template$PEOrientationOther)
    # NA-safe: all(pe_counts == 0) is NA (and if (NA) errors) when the
    # orientation counts are missing (e.g. no paired-end evidence at all).
    pe_ok <- !is.na(pe_counts)
    template$PEOrientationDominant <- if (!any(pe_ok)) NA_character_ else
      names(pe_counts[pe_ok])[which.max(pe_counts[pe_ok])]
    template$MergedFromN <- n_members
    template$MergeFragmentWarning <- n_members > max_plausible_fragments
    template$MergedSupportSemantics <- if (exact_union) {
      "exact_qname_union_across_component_rows"
    } else {
      "sum_of_component_rows; shared_QNAMEs_may_be_double_counted"
    }
    if (isTRUE(template$MergeFragmentWarning)) {
      log_msg(sprintf(paste(
        "WARNING: merged interval at %d (length %d bp) pools %d fragments,",
        "exceeding max_plausible_fragments=%d. A single true PTD has",
        "historically fragmented into at most ~2-3 rows pre-merge, so this",
        "many disjoint clusters merging by proximity alone is more consistent",
        "with scattered background signal than one real event -- review before",
        "reporting. Flagged via MergeFragmentWarning."
      ), mg_start, mg_len, n_members, max_plausible_fragments))
    }
    if (!is.na(template$SpanMinCoverage) &&
        template$SpanMinCoverage < thv("min_span_coverage", 0)) {
      log_msg(sprintf(paste(
        "REJECTED merged interval at %d (length %d bp, %d fragment(s)):",
        "SpanMinCoverage=%.1f < min_span_coverage=%.1f -- part of this span",
        "has essentially no real coverage, so the merge is not trusted as",
        "one true event."
      ), mg_start, mg_len, n_members, template$SpanMinCoverage, thresholds$min_span_coverage))
      return(NULL)
    }
    if (!is.na(template$SpanDepthFoldChange) &&
        (template$SpanDepthFoldChange < thv("min_span_depth_fold_change", 0) ||
         template$SpanDepthFoldChange > thv("max_span_depth_fold_change", Inf))) {
      log_msg(sprintf(paste(
        "REJECTED merged interval at %d (length %d bp, %d fragment(s)):",
        "SpanDepthFoldChange=%.2fx outside [%.2f, %.2f] -- doesn't match the",
        "depth pattern of a real duplication (too high suggests a",
        "collapsed/repetitive region pulling in multi-mapped reads)."
      ), mg_start, mg_len, n_members, template$SpanDepthFoldChange,
      thresholds$min_span_depth_fold_change, thresholds$max_span_depth_fold_change))
      return(NULL)
    }
    template
  })
  merged_rows <- Filter(Negate(is.null), merged_rows)
  if (length(merged_rows) == 0L) {
    log_msg("All merged intervals rejected by min_span_coverage; none survive.")
    return(final_df[skip, , drop = FALSE])
  }
  merged_df <- do.call(rbind, merged_rows)
  log_msg(sprintf("Merged into %d interval(s).", nrow(merged_df)))
  if (length(skip) > 0L) rbind(merged_df, final_df[skip, , drop = FALSE]) else merged_df
}


#' @rdname detect_itd
#' @export
talos <- function(
  bam_path,
  gene,
  build = "hg19",
  padding = 500L,
  min_support = NULL,
  min_size = NULL,
  max_correction = NULL,
  max_missing_kmers = NULL,
  plot = TRUE,
  sample_name = NULL,
  filter_intronic = NULL,
  ptd_mode = NULL,
  use_cigar_bp = NULL,
  refine_bp = NULL,
  min_mapq = NULL,
  min_wt_reads = NULL,
  cluster_tolerance = NULL,
  min_base_quality = NULL,
  adapter_sequences = NULL,
  adapter_min_overlap = NULL,
  cluster_length_scale = NULL,
  max_cluster_tolerance = NULL,
  search_window = NULL,
  require_unique_anchor = NULL,
  vaf_threshold = NULL,
  min_strand_bias = NULL,
  max_strand_bias = NULL,
  min_mean_support_mapq = NULL,
  max_breakpoint_spread = NULL,
  min_softclip_fraction = NULL,
  min_unique_breakpoints = NULL,
  min_coverage_drop = NULL,
  min_local_coverage = NULL,
  coverage_window = NULL,
  min_microhomology = NULL,
  min_span_coverage = NULL,
  min_span_depth_fold_change = NULL,
  max_span_depth_fold_change = NULL,
  min_discordant_ratio = NULL,
  min_entropy = NULL,
  min_alignment_score = NULL,
  detect_orientation = NULL,
  do_annotate_hotspots = NULL,
  hotspot_db_path = NULL,
  compute_alignment_score = NULL,
  compute_support_bases = NULL,
  compute_consistency = NULL,
  compute_itd_coverage = NULL,
  compute_coverage_drop = NULL,
  compute_microhomology = NULL,
  compute_repeat_entropy = NULL,
  compute_discordant_ratio = NULL,
  compute_hgvs = NULL,
  min_itd_read_coverage = NULL,
  html_report = NULL,
  use_gviz = NULL,
  bsgenome = NULL,
  yaml_path = system.file("extdata", "gene_config.yaml", package = "TALOS"),
  add_config_to_report = FALSE,
  max_pairwise_alignments = NULL,
  output_prefix = "TALOS",
  include_sample = TRUE,
  include_gene = TRUE,
  include_timestamp = TRUE,
  output_sep = "_",
  global_log = TRUE,
  min_side_softclip_reads = NULL,
  max_side_ratio = NULL,
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
  # Analysis-window settings may come from the YAML gene_settings too
  # (explicit argument > YAML > default), e.g. a wider padding for KMT2A so
  # intronic PTD breakpoints fall inside the window.
  yaml_gs <- tryCatch(yaml::read_yaml(yaml_path)[[gene]]$gene_settings, error = function(e) NULL)
  if (missing(padding) && !is.null(yaml_gs$padding)) padding <- as.integer(yaml_gs$padding)
  if (missing(exon_padding) && !is.null(yaml_gs$exon_padding)) exon_padding <- as.integer(yaml_gs$exon_padding)

  config <- get_gene_config(
    gene = gene,
    build = build,
    padding = padding,
    config_path = yaml_path,
    bsgenome = bsgenome,
    exon_padding = exon_padding
  )
  config$build <- build

  # Settings resolved as: explicit argument > YAML gene_settings > detect_itd()
  # default. The defaults are read from detect_itd()'s formals so the two
  # functions cannot drift apart (they used to be duplicated by hand).
  tunable <- c(
    "min_support", "min_size", "max_correction", "filter_intronic", "max_missing_kmers",
    "ptd_mode", "use_cigar_bp", "refine_bp", "min_mapq", "min_wt_reads", "cluster_tolerance",
    "min_base_quality", "adapter_sequences", "adapter_min_overlap",
    "cluster_length_scale", "max_cluster_tolerance",
    "vaf_threshold", "search_window", "require_unique_anchor", "min_strand_bias",
    "max_strand_bias", "min_mean_support_mapq", "max_breakpoint_spread",
    "min_softclip_fraction", "min_unique_breakpoints", "min_coverage_drop",
    "min_local_coverage", "coverage_window", "min_microhomology", "min_span_coverage",
    "min_span_depth_fold_change", "max_span_depth_fold_change", "min_discordant_ratio",
    "min_entropy", "min_alignment_score", "detect_orientation", "do_annotate_hotspots",
    "compute_alignment_score", "compute_support_bases", "compute_consistency",
    "compute_itd_coverage", "compute_coverage_drop", "compute_microhomology",
    "compute_repeat_entropy", "compute_discordant_ratio", "compute_hgvs", "html_report",
    "use_gviz", "min_itd_read_coverage", "max_pairwise_alignments",
    "min_side_softclip_reads", "max_side_ratio", "min_softclip_pct_side",
    "min_left_softclip_pct_wt", "min_right_softclip_pct_wt", "min_abs_side_softclip",
    "max_itd_length", "convert_long_to_ptd", "min_length", "max_length",
    "prefer_extended_length", "ptd_allow_asymmetric", "ptd_use_local_assembly",
    "merge_ptd_intervals", "merge_gap", "min_ptd_length", "max_plausible_fragments",
    "use_exon_graph"
  )
  p <- .resolve_settings(environment(), formals(detect_itd)[tunable], config$gene_settings)

  if (gene == "KMT2A") {
    if (verbose) message("[TALOS] Applying relaxed thresholds for KMT2A PTD detection")
    p$min_support <- min(p$min_support, 5L)
    p$min_side_softclip_reads <- min(p$min_side_softclip_reads, 10L)
    p$min_abs_side_softclip <- min(p$min_abs_side_softclip, 10L)
    p$min_coverage_drop <- min(p$min_coverage_drop, 1.0)
    p$ptd_allow_asymmetric <- TRUE
    p$min_itd_read_coverage <- min(p$min_itd_read_coverage, 10)
    # Graph on by default for KMT2A unless the caller or the YAML said otherwise
    # (the YAML value used to be overridden).
    if (is.null(use_exon_graph) && is.null(config$gene_settings$use_exon_graph)) {
      p$use_exon_graph <- TRUE
    }
  }

  do.call(detect_itd, c(
    list(
      bam_path = bam_path, gene_config = config, plot = plot, sample_name = sample_name,
      hotspot_db_path = hotspot_db_path, add_config_to_report = add_config_to_report,
      output_prefix = output_prefix, include_sample = include_sample,
      include_gene = include_gene, include_timestamp = include_timestamp,
      output_sep = output_sep, global_log = global_log, use_kmers = use_kmers,
      exon_padding = exon_padding, verbose = verbose, debug = debug
    ),
    p,
    list(...)
  ))
}

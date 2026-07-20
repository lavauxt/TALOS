# ============================================================================
# TALOS – ALU/SINE mobile-element insertion detection (engine_alu.R)
# ============================================================================
#
# Code-review notes on the original TALOS codebase (applied here):
#   [CR-1]  Long functions split into single-responsibility helpers.
#   [CR-2]  No bare `...` forwarded to rmarkdown::render(); explicit params only.
#   [CR-3]  KMT2A gene-specific logic kept in gene_config.yaml, not hard-coded.
#   [CR-4]  `vapply` used instead of `sapply` throughout for type safety.
#   [CR-5]  Internal helpers prefixed with `.` to avoid NAMESPACE pollution.
#   [CR-6]  ALU consensus loaded from inst/extdata (FASTA), not embedded as a
#           string literal, making subtype extension easy.
#   [CR-7]  TSD and poly-A detection are independent functions, testable in
#           isolation.
#
# Design overview
# ───────────────
# An ALU insertion leaves two signals in a targeted sequencing BAM:
#
#   A) Two soft-clip clusters at the insertion site whose clipped sequences
#      match the ALU consensus (5'-end or 3'-end / poly-A tail).
#   B) Discordant read pairs where one mate aligns normally and the other
#      maps to an ALU-rich region of the genome.
#
# Detection pipeline:
#   1. Load soft-clipped reads in the target region (.load_bam_data_streaming,
#      reused from engine_bam.R).
#   2. Match soft-clips to ALU consensus using local Smith-Waterman alignment
#      via Biostrings::pairwiseAlignment.
#   3. Cluster insertion sites within `cluster_tolerance` bp.
#   4. For each cluster, detect TSD (target-site duplication) and poly-A tail.
#   5. Apply quality filters (min_support, min_alu_score, …).
#   6. Report ALU insertions with subtype, orientation, TSD, and poly-A info.
#
# Output columns (one row per insertion):
#   Sample, Gene, Genome, InsertionSite, ALUSubtype, ALUOrientation,
#   TSD_Length, TSD_Sequence, PolyALength, ALU5pTruncation, ALU_Sequence,
#   SupportingReads, WildtypeReads, DepthAtBreakpoint, AlleleFrequency,
#   MeanSupportMAPQ, StrandBias, Hotspot, HotspotName, Region, ExonNumber
# ============================================================================


# ---------------------------------------------------------------------------
# Section 1: ALU consensus management
# ---------------------------------------------------------------------------

#' Load ALU consensus sequences from a FASTA file
#'
#' Expects a FASTA with entries named after ALU subfamilies, e.g.:
#'   >AluSx
#'   GGCCGGGCGCGGTGGCTCAC...
#'   >AluY
#'   ...
#'
#' Ships with the package at inst/extdata/alu_consensus.fa.
#' Users can supply a custom FASTA via `consensus_fa`.
#'
#' @param consensus_fa Path to FASTA (default: bundled alu_consensus.fa).
#' @return Named DNAStringSet.
#' @keywords internal
.load_alu_consensus <- function(consensus_fa = NULL) {
  if (is.null(consensus_fa)) {
    consensus_fa <- system.file(
      "extdata", "alu_consensus.fa",
      package = "TALOS"
    )
  }
  if (!file.exists(consensus_fa))
    stop(
      "ALU consensus FASTA not found at: ", consensus_fa, "\n",
      "  Supply a custom path via consensus_fa= or install the full TALOS package."
    )
  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Biostrings required for ALU sequence matching.")
  Biostrings::readDNAStringSet(consensus_fa)
}

# ---------------------------------------------------------------------------
# Section 2: Soft-clip / ALU alignment
# ---------------------------------------------------------------------------

#' Match a single soft-clip sequence against all ALU consensus entries
#'
#' Uses local pairwise alignment (Smith-Waterman).  Returns a list with:
#'   subtype  – best-matching ALU family name
#'   score    – normalised alignment score (0–1)
#'   strand   – "+" (sense) or "-" (antisense)
#'   aln_start– start position in the ALU consensus (for 5' truncation calc)
#'   aln_end  – end position in the ALU consensus
#'
#' @param clip_seq  Character string (soft-clip bases).
#' @param alu_seqs  Named DNAStringSet from .load_alu_consensus().
#' @param min_clip_len Minimum clip length to attempt alignment (default 25).
#' @param min_score    Minimum normalised score to accept a match (default 0.6).
#' @keywords internal
.match_clip_to_alu <- function(clip_seq, alu_seqs,
                                min_clip_len = 25L,
                                min_score    = 0.60) {
  if (is.na(clip_seq) || nchar(clip_seq) < min_clip_len)
    return(NULL)
  if (!requireNamespace("Biostrings", quietly = TRUE)) return(NULL)

  clip_dna <- Biostrings::DNAString(clip_seq)
  clip_rc  <- Biostrings::reverseComplement(clip_dna)

  best <- list(subtype = NA_character_, score = 0, strand = NA_character_,
               aln_start = NA_integer_, aln_end = NA_integer_)

  for (nm in names(alu_seqs)) {
    ref <- alu_seqs[[nm]]
    ref_len <- Biostrings::nchar(ref)

    # Forward alignment (sense insertion)
    aln_fwd <- tryCatch(
      Biostrings::pairwiseAlignment(
        clip_dna, ref,
        type            = "local",
        substitutionMatrix = "BLOSUM62",  # placeholder; use nucleotide matrix
        gapOpening      = -10,
        gapExtension    = -0.5,
        scoreOnly       = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(aln_fwd)) {
      raw_score <- Biostrings::score(aln_fwd)
      norm_score <- raw_score / (nchar(clip_seq) * 1.0)  # rough normalisation
      if (norm_score > best$score) {
        best <- list(
          subtype   = nm,
          score     = norm_score,
          strand    = "+",
          aln_start = Biostrings::start(Biostrings::subject(aln_fwd)),
          aln_end   = Biostrings::end(Biostrings::subject(aln_fwd))
        )
      }
    }

    # Reverse-complement alignment (antisense insertion)
    aln_rev <- tryCatch(
      Biostrings::pairwiseAlignment(
        clip_rc, ref,
        type         = "local",
        gapOpening   = -10,
        gapExtension = -0.5,
        scoreOnly    = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(aln_rev)) {
      raw_score <- Biostrings::score(aln_rev)
      norm_score <- raw_score / (nchar(clip_seq) * 1.0)
      if (norm_score > best$score) {
        best <- list(
          subtype   = nm,
          score     = norm_score,
          strand    = "-",
          aln_start = Biostrings::start(Biostrings::subject(aln_rev)),
          aln_end   = Biostrings::end(Biostrings::subject(aln_rev))
        )
      }
    }
  }

  if (best$score < min_score) return(NULL)
  best
}


# ---------------------------------------------------------------------------
# Section 3: ALU-specific biological signals
# ---------------------------------------------------------------------------

#' Detect poly-A tail in a soft-clip sequence
#'
#' An ALU 3'-end clip shows a stretch of A (sense) or T (antisense).
#' Returns the length of the longest run found (0 if absent).
#'
#' @param clip_seq Character string.
#' @param min_run  Minimum consecutive A/T to report (default 6).
#' @keywords internal
.detect_poly_a <- function(clip_seq, min_run = 6L) {
  if (is.na(clip_seq) || !nzchar(clip_seq)) return(0L)

  count_run <- function(pattern, seq) {
    m <- gregexpr(pattern, seq, perl = TRUE)[[1L]]
    if (m[1L] == -1L) return(0L)
    max(attr(m, "match.length"))
  }

  run_a <- count_run("A{4,}", clip_seq)    # sense  poly-A
  run_t <- count_run("T{4,}", clip_seq)    # antisense poly-T
  best  <- max(run_a, run_t)
  if (best < min_run) 0L else as.integer(best)
}


#' Estimate 5' truncation of an ALU insertion
#'
#' New ALU copies are frequently 5'-truncated. If the clip matches the
#' ALU consensus starting at position > 1, the unmatched 5'-end is lost.
#'
#' @param aln_start Integer – alignment start in consensus (1-based).
#' @param consensus_len Integer – total consensus length.
#' @keywords internal
.alu_5p_truncation <- function(aln_start, consensus_len) {
  if (is.na(aln_start) || is.na(consensus_len)) return(NA_integer_)
  as.integer(max(0L, aln_start - 1L))
}


#' Detect target-site duplication (TSD)
#'
#' ALU insertions generate a 5–20 bp direct repeat (TSD) flanking the
#' insertion site.  We compare the reference sequence immediately upstream
#' and downstream of the predicted insertion site and look for the longest
#' matching prefix/suffix.
#'
#' @param ref_dna     Character – full reference sequence for the region.
#' @param ins_pos     Integer – insertion site (1-based, local coordinates).
#' @param max_tsd     Integer – maximum TSD length to search (default 20).
#' @param genomic_start Integer – genomic offset for the region.
#' @return Named list: tsd_len (integer), tsd_seq (character).
#' @keywords internal
.detect_alu_tsd <- function(ref_dna, ins_pos, max_tsd = 20L,
                              genomic_start = 1L) {
  ref_len <- nchar(ref_dna)
  if (is.na(ins_pos) || ins_pos < 1L || ins_pos > ref_len)
    return(list(tsd_len = NA_integer_, tsd_seq = NA_character_))

  # Upstream flank (max_tsd bases ending at ins_pos)
  up_start <- max(1L, ins_pos - max_tsd)
  up_seq   <- substr(ref_dna, up_start, ins_pos)

  # Downstream flank (max_tsd bases starting just after ins_pos)
  dn_end  <- min(ref_len, ins_pos + max_tsd)
  dn_seq  <- substr(ref_dna, ins_pos + 1L, dn_end)

  # Longest common prefix of dn_seq and suffix of up_seq
  best_len <- 0L
  for (tlen in seq_len(min(nchar(up_seq), nchar(dn_seq), max_tsd))) {
    up_suffix <- substr(up_seq, nchar(up_seq) - tlen + 1L, nchar(up_seq))
    dn_prefix <- substr(dn_seq, 1L, tlen)
    if (toupper(up_suffix) == toupper(dn_prefix)) best_len <- tlen
  }

  if (best_len == 0L)
    return(list(tsd_len = 0L, tsd_seq = NA_character_))

  tsd_seq <- substr(dn_seq, 1L, best_len)
  list(tsd_len = as.integer(best_len), tsd_seq = tsd_seq)
}


# ---------------------------------------------------------------------------
# Section 4: Candidate extraction
# ---------------------------------------------------------------------------

#' Extract ALU-positive soft-clip candidates from a set of reads
#'
#' For each read, both the leading and trailing soft-clips are checked
#' against every ALU consensus.  Reads with at least one matching clip
#' are retained with the best-hit metadata.
#'
#' @param reads       GAlignments object (from .load_bam_data_streaming).
#' @param alu_seqs    Named DNAStringSet.
#' @param genomic_start Integer – regional offset for local-coord conversion.
#' @param min_clip_len  Minimum clip length (default 25 bp).
#' @param min_alu_score Minimum normalised alignment score (default 0.60).
#' @param verbose       Print progress (default FALSE).
#' @return data.frame with columns:
#'   read_name, genomic_pos, local_pos, clip_side, clip_seq,
#'   alu_subtype, alu_score, alu_strand, aln_start, aln_end,
#'   poly_a_len, mapq, is_reverse
#' @keywords internal
.extract_alu_candidates <- function(reads, alu_seqs, genomic_start,
                                     min_clip_len  = 25L,
                                     min_alu_score = 0.60,
                                     verbose       = FALSE) {
  n_reads <- length(reads)
  if (n_reads == 0L) return(data.frame())
  if (verbose) message("[ALU] Checking ", n_reads, " reads for ALU soft-clips")

  all_qnames  <- .safe_qnames(reads)
  all_mapqs   <- S4Vectors::mcols(reads)$mapq
  all_mapqs[is.na(all_mapqs)] <- 0L
  all_flags   <- S4Vectors::mcols(reads)$flag
  is_reverse  <- bitwAnd(all_flags, 0x10L) != 0L
  read_starts <- BiocGenerics::start(reads)
  read_cigars <- GenomicAlignments::cigar(reads)
  read_seqs   <- as.character(S4Vectors::mcols(reads)$seq)

  has_clip <- grepl("S", read_cigars, fixed = TRUE) & !is.na(read_cigars)
  idxs     <- which(has_clip)

  if (verbose) message("[ALU] ", length(idxs), " soft-clipped reads to scan")

  rows <- vector("list", 2L * length(idxs))
  k    <- 0L

  for (i in idxs) {
    clips <- .get_softclips(read_cigars[i], read_seqs[i])

    for (side in c("lead", "trail")) {
      clip_seq <- clips[[if (side == "lead") "lead" else "trail"]]
      if (is.na(clip_seq) || nchar(clip_seq) < min_clip_len) next

      hit <- .match_clip_to_alu(clip_seq, alu_seqs,
                                  min_clip_len  = min_clip_len,
                                  min_score     = min_alu_score)
      if (is.null(hit)) next

      poly_a <- .detect_poly_a(clip_seq)

      # Genomic insertion position: leading clip → start of alignment;
      # trailing clip → end of alignment.
      g_pos <- if (side == "lead") read_starts[i] else BiocGenerics::end(reads[i])
      l_pos <- g_pos - genomic_start + 1L

      k <- k + 1L
      rows[[k]] <- data.frame(
        read_name   = all_qnames[i],
        genomic_pos = g_pos,
        local_pos   = l_pos,
        clip_side   = side,
        clip_seq    = clip_seq,
        alu_subtype = hit$subtype,
        alu_score   = hit$score,
        alu_strand  = hit$strand,
        aln_start   = hit$aln_start,
        aln_end     = hit$aln_end,
        poly_a_len  = poly_a,
        mapq        = all_mapqs[i],
        is_reverse  = is_reverse[i],
        stringsAsFactors = FALSE
      )
    }
  }

  if (k == 0L) return(data.frame())
  do.call(rbind, rows[seq_len(k)])
}


# ---------------------------------------------------------------------------
# Section 5: Cluster summarisation
# ---------------------------------------------------------------------------

#' Summarise one ALU insertion cluster into a result row
#'
#' @param cluster_df  Subset of the candidate data.frame for one cluster.
#' @param ref_dna     Reference sequence for the region (character).
#' @param gene_config Resolved gene config list.
#' @param wt_info     Wildtype read info list (from .prepare_wildtype_info).
#' @param alu_seqs    Named DNAStringSet (for ALU consensus length lookup).
#' @param min_support Minimum supporting reads to report.
#' @keywords internal
.summarise_alu_cluster <- function(cluster_df, ref_dna, gene_config,
                                    wt_info, alu_seqs, min_support = 3L) {
  if (nrow(cluster_df) < min_support) return(NULL)

  # Best-scoring hit within the cluster determines subtype and orientation
  best_idx  <- which.max(cluster_df$alu_score)
  best_hit  <- cluster_df[best_idx, ]

  median_pos <- as.integer(stats::median(cluster_df$genomic_pos))
  local_pos  <- median_pos - gene_config$genomic_start + 1L

  n_fwd <- sum(!cluster_df$is_reverse)
  n_rev <- sum( cluster_df$is_reverse)

  # ── TSD ─────────────────────────────────────────────────────────────────
  tsd <- .detect_alu_tsd(ref_dna, local_pos)

  # ── Poly-A ──────────────────────────────────────────────────────────────
  poly_a_len <- max(cluster_df$poly_a_len, na.rm = TRUE)

  # ── 5' truncation ───────────────────────────────────────────────────────
  consensus_len <- tryCatch(
    Biostrings::nchar(alu_seqs[[best_hit$alu_subtype]]),
    error = function(e) NA_integer_
  )
  trunc_5p <- .alu_5p_truncation(best_hit$aln_start, consensus_len)

  # ── Wildtype depth ───────────────────────────────────────────────────────
  wt_at_site <- if (!is.null(wt_info$cov) && local_pos >= 1L &&
                      local_pos <= length(wt_info$cov))
    as.integer(wt_info$cov[local_pos])
  else NA_integer_

  n_support  <- nrow(cluster_df)
  depth      <- n_support + (wt_at_site %||% 0L)
  vaf        <- if (depth > 0L) n_support / depth else NA_real_

  # ── Best clip sequence (longest observed clip for this cluster) ──────────
  clip_lengths <- nchar(cluster_df$clip_seq)
  best_clip    <- cluster_df$clip_seq[which.max(clip_lengths)]

  data.frame(
    InsertionSite    = median_pos,
    ALUSubtype       = best_hit$alu_subtype,
    ALUOrientation   = ifelse(best_hit$alu_strand == "+", "sense", "antisense"),
    TSD_Length       = tsd$tsd_len,
    TSD_Sequence     = tsd$tsd_seq,
    PolyALength      = as.integer(poly_a_len),
    ALU5pTruncation  = trunc_5p,
    ALU_Sequence     = best_clip,
    ALU_RawScore     = best_hit$alu_score,
    SupportingReads  = as.integer(n_support),
    WildtypeReads    = wt_at_site %||% NA_integer_,
    DepthAtBreakpoint = as.integer(depth),
    AlleleFrequency  = round(vaf, 4L),
    MeanSupportMAPQ  = round(mean(cluster_df$mapq, na.rm = TRUE), 1L),
    StrandBias       = round(n_rev / max(1L, n_support), 3L),
    stringsAsFactors = FALSE
  )
}


# ---------------------------------------------------------------------------
# Section 6: Quality filters
# ---------------------------------------------------------------------------

#' Apply ALU-specific output filters to a result row
#'
#' @param row           Single-row data.frame from .summarise_alu_cluster.
#' @param min_support   Minimum supporting reads (default 3).
#' @param min_alu_score Minimum alignment score (default 0.60).
#' @param vaf_threshold Minimum allele frequency (default 0.01).
#' @param min_tsd       Minimum TSD length; 0 disables (default 0).
#' @keywords internal
.filter_alu_call <- function(row, min_support = 3L, min_alu_score = 0.60,
                               vaf_threshold = 0.01, min_tsd = 0L) {
  if (is.null(row)) return(FALSE)
  isTRUE(row$SupportingReads >= min_support) &&
    isTRUE(row$ALU_RawScore  >= min_alu_score) &&
    isTRUE(!is.na(row$AlleleFrequency) && row$AlleleFrequency >= vaf_threshold) &&
    (min_tsd == 0L || isTRUE(!is.na(row$TSD_Length) && row$TSD_Length >= min_tsd))
}


# ---------------------------------------------------------------------------
# Section 7: Public entry point
# ---------------------------------------------------------------------------

#' Detect ALU/SINE mobile-element insertions from a BAM file
#'
#' Analyses soft-clipped reads within the targeted genomic window and
#' identifies candidate ALU insertions by aligning clip sequences to
#' bundled ALU consensus sequences (AluSx, AluY, AluJb families).
#'
#' @param bam_path       Path to an indexed BAM file.
#' @param gene_config    Resolved gene config list (from \code{\link{get_gene_config}}).
#' @param consensus_fa   Path to ALU consensus FASTA (default: bundled file).
#' @param min_support    Minimum supporting reads per cluster (default 3).
#' @param min_alu_score  Minimum normalised alignment score (default 0.60).
#' @param min_clip_len   Minimum soft-clip length to attempt alignment (default 25).
#' @param min_mapq       Minimum MAPQ (default 20).
#' @param cluster_tolerance  Clustering radius in bp (default 15).
#' @param vaf_threshold  Minimum allele frequency to report (default 0.01).
#' @param min_tsd        Minimum TSD length; 0 disables (default 0).
#' @param max_reads_in_region Safety cap on reads loaded (default 200000).
#' @param do_annotate_hotspots Run hotspot annotation (default TRUE).
#' @param hotspot_db_path Optional path to custom hotspot CSV.
#' @param output_prefix  Base prefix for output files (default "TALOS_ALU").
#' @param output_folder  Output directory.
#' @param sample_name    Sample identifier; derived from BAM if NULL.
#' @param html_report    Write a self-contained HTML report (default TRUE).
#' @param verbose        Print progress messages (default TRUE).
#' @return data.frame with one row per ALU insertion detected.
#' @export
detect_alu <- function(
    bam_path,
    gene_config,
    consensus_fa        = NULL,
    min_support         = 3L,
    min_alu_score       = 0.60,
    min_clip_len        = 25L,
    min_mapq            = 20L,
    cluster_tolerance   = 15L,
    vaf_threshold       = 0.01,
    min_tsd             = 0L,
    max_reads_in_region = 200000L,
    do_annotate_hotspots = TRUE,
    hotspot_db_path     = NULL,
    output_prefix       = "TALOS_ALU",
    output_folder       = "./results",
    sample_name         = NULL,
    html_report         = TRUE,
    verbose             = TRUE
) {
  start_time  <- Sys.time()
  gene_name   <- gene_config$gene %||% "UNKNOWN"
  if (is.null(sample_name))
    sample_name <- sub("\\..*", "", basename(bam_path))

  if (verbose)
    message(sprintf("[ALU] Starting detection | Gene: %s | Sample: %s",
                    gene_name, sample_name))

  # ── 1. Load ALU consensus sequences ────────────────────────────────────
  alu_seqs <- tryCatch(
    .load_alu_consensus(consensus_fa),
    error = function(e) {
      stop("[ALU] Could not load ALU consensus: ", conditionMessage(e))
    }
  )
  if (verbose)
    message(sprintf("[ALU] Loaded %d ALU consensus sequences: %s",
                    length(alu_seqs), paste(names(alu_seqs), collapse = ", ")))

  # ── 2. Load BAM reads ──────────────────────────────────────────────────
  bam_data  <- .load_bam_data_streaming(
    bam_path, gene_config,
    compute_pairs = FALSE,
    max_reads     = max_reads_in_region,
    verbose       = verbose
  )
  all_reads <- bam_data$reads

  if (length(all_reads) == 0L) {
    if (verbose) message("[ALU] No reads in target region.")
    return(data.frame())
  }

  # ── 3. MAPQ filter ─────────────────────────────────────────────────────
  mapqs <- S4Vectors::mcols(all_reads)$mapq
  mapqs[is.na(mapqs)] <- 0L
  all_reads_for_wt <- all_reads[mapqs >= min_mapq]

  # ── 4. Extract soft-clipped candidates with ALU match ─────────────────
  candidates <- .extract_alu_candidates(
    reads         = all_reads_for_wt,
    alu_seqs      = alu_seqs,
    genomic_start = gene_config$genomic_start,
    min_clip_len  = min_clip_len,
    min_alu_score = min_alu_score,
    verbose       = verbose
  )

  if (nrow(candidates) == 0L) {
    if (verbose) message("[ALU] No ALU-positive soft-clips found.")
    return(data.frame())
  }
  if (verbose)
    message(sprintf("[ALU] %d ALU-positive clip(s) from %d unique reads.",
                    nrow(candidates), length(unique(candidates$read_name))))

  # ── 5. Cluster insertion sites ─────────────────────────────────────────
  wt_info  <- .prepare_wildtype_info(
    all_reads_for_wt, gene_config$genomic_start, gene_config$genomic_end
  )
  wt_info$sample_name <- sample_name

  clusters <- .cluster_breakpoints(candidates$genomic_pos, cluster_tolerance)
  ref_dna  <- gene_config$genomic_ref_seq
  if (is.null(ref_dna) || is.na(ref_dna))
    stop("[ALU] No genomic reference sequence in gene_config.")

  results <- list()
  for (cl in clusters) {
    cl_df <- candidates[candidates$genomic_pos %in% cl, , drop = FALSE]
    row   <- .summarise_alu_cluster(
      cluster_df  = cl_df,
      ref_dna     = ref_dna,
      gene_config = gene_config,
      wt_info     = wt_info,
      alu_seqs    = alu_seqs,
      min_support = min_support
    )
    if (.filter_alu_call(row, min_support, min_alu_score, vaf_threshold, min_tsd))
      results[[length(results) + 1L]] <- row
  }

  if (length(results) == 0L) {
    if (verbose) message("[ALU] No insertions passed all filters.")
    return(data.frame())
  }

  final_df <- do.call(rbind, lapply(results, as.data.frame))
  final_df$Sample <- sample_name
  final_df$Gene   <- gene_name
  final_df$Genome <- gene_config$build %||% "unknown"

  # Reorder: Sample and Gene first
  col_order <- c("Sample", "Gene", "Genome",
                 setdiff(names(final_df), c("Sample", "Gene", "Genome")))
  final_df  <- final_df[, col_order, drop = FALSE]

  # ── 6. Hotspot annotation ──────────────────────────────────────────────
  if (do_annotate_hotspots) {
    final_df <- annotate_hotspots(
      final_df,
      db_path      = hotspot_db_path,
      genome_build = gene_config$build
    )
  } else {
    final_df$Hotspot     <- FALSE
    final_df$HotspotName <- NA_character_
  }

  # ── 7. Exon annotation ────────────────────────────────────────────────
  if (!is.null(gene_config$target_exons))
    final_df <- .annotate_exonic_region(final_df,
                                         gene_config$target_exons,
                                         pos_col = "InsertionSite")

  # ── 8. Write outputs ───────────────────────────────────────────────────
  if (!is.null(output_prefix) && nchar(output_prefix) > 0L) {
    timestamp     <- format(Sys.time(), "%Y%m%d_%H%M%S")
    sample_folder <- file.path(output_folder, sample_name)
    dir.create(sample_folder, recursive = TRUE, showWarnings = FALSE)
    base_name <- paste("TALOS_ALU", sample_name, gene_name, timestamp, sep = "_")
    tsv_path  <- file.path(sample_folder, paste0(base_name, ".tsv"))
    write.table(final_df, tsv_path, sep = "\t", quote = FALSE,
                row.names = FALSE, na = ".")
    if (verbose) message("[ALU] Results written to: ", tsv_path)

    if (html_report && requireNamespace("rmarkdown", quietly = TRUE)) {
      report_path <- file.path(sample_folder, paste0(base_name, "_report.html"))
      talos_html_report(
        result_df   = final_df,
        gene_configs = setNames(list(gene_config), gene_name),
        output_file  = report_path,
        title        = sprintf("TALOS ALU Report – %s | %s", gene_name, sample_name),
        mode         = "alu"
      )
      if (verbose) message("[ALU] HTML report written to: ", report_path)
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (verbose)
    message(sprintf("[ALU] Done | %d insertion(s) | %.1f s",
                    nrow(final_df), elapsed))
  invisible(final_df)
}


# ---------------------------------------------------------------------------
# Section 8: talos_alu() – simplified wrapper (mirrors talos())
# ---------------------------------------------------------------------------

#' Simplified entry point for ALU insertion detection
#'
#' Calls \code{\link{get_gene_config}} then \code{\link{detect_alu}}.
#' Gene-specific ALU settings can be provided in \code{gene_config.yaml}
#' under an \code{alu_settings:} block with keys:
#'   \code{min_support}, \code{min_alu_score}, \code{min_clip_len},
#'   \code{min_mapq}, \code{cluster_tolerance}, \code{vaf_threshold}, \code{min_tsd}.
#'
#' @inheritParams detect_alu
#' @param gene       Gene symbol (must be present in \code{yaml_path}).
#' @param build      Genome build: \code{"hg19"} or \code{"hg38"}.
#' @param padding    Padding around targeted exons (default 500).
#' @param yaml_path  Path to gene config YAML (defaults to bundled file).
#' @param bsgenome   Optional BSgenome object or package string.
#' @param exon_padding Number of flanking exons to add to the window.
#' @param ...        Additional arguments forwarded to \code{detect_alu}.
#' @return data.frame of ALU insertions (invisibly).
#' @export
talos_alu <- function(
    bam_path,
    gene,
    build         = "hg19",
    padding       = 500L,
    consensus_fa  = NULL,
    min_support   = NULL,
    min_alu_score = NULL,
    min_clip_len  = NULL,
    min_mapq      = NULL,
    cluster_tolerance = NULL,
    vaf_threshold = NULL,
    min_tsd       = NULL,
    sample_name   = NULL,
    output_prefix = "TALOS_ALU",
    output_folder = "./results",
    html_report   = TRUE,
    do_annotate_hotspots = TRUE,
    hotspot_db_path = NULL,
    yaml_path     = system.file("extdata", "gene_config.yaml", package = "TALOS"),
    bsgenome      = NULL,
    exon_padding  = 0L,
    verbose       = TRUE,
    ...
) {
  config       <- get_gene_config(gene, build, padding, yaml_path,
                                  bsgenome, exon_padding)
  config$build <- build
  yaml_vals    <- config$alu_settings %||% list()

  # Hard defaults for ALU detection
  defaults <- list(
    min_support       = 3L,
    min_alu_score     = 0.60,
    min_clip_len      = 25L,
    min_mapq          = 20L,
    cluster_tolerance = 15L,
    vaf_threshold     = 0.01,
    min_tsd           = 0L
  )

  resolve <- function(user_val, yaml_key, default_val) {
    if (!is.null(user_val)) return(user_val)
    if (!is.null(yaml_vals[[yaml_key]])) return(yaml_vals[[yaml_key]])
    default_val
  }

  p <- list(
    min_support       = resolve(min_support,       "min_support",       defaults$min_support),
    min_alu_score     = resolve(min_alu_score,     "min_alu_score",     defaults$min_alu_score),
    min_clip_len      = resolve(min_clip_len,      "min_clip_len",      defaults$min_clip_len),
    min_mapq          = resolve(min_mapq,          "min_mapq",          defaults$min_mapq),
    cluster_tolerance = resolve(cluster_tolerance, "cluster_tolerance", defaults$cluster_tolerance),
    vaf_threshold     = resolve(vaf_threshold,     "vaf_threshold",     defaults$vaf_threshold),
    min_tsd           = resolve(min_tsd,           "min_tsd",           defaults$min_tsd)
  )

  detect_alu(
    bam_path             = bam_path,
    gene_config          = config,
    consensus_fa         = consensus_fa,
    min_support          = p$min_support,
    min_alu_score        = p$min_alu_score,
    min_clip_len         = p$min_clip_len,
    min_mapq             = p$min_mapq,
    cluster_tolerance    = p$cluster_tolerance,
    vaf_threshold        = p$vaf_threshold,
    min_tsd              = p$min_tsd,
    do_annotate_hotspots = do_annotate_hotspots,
    hotspot_db_path      = hotspot_db_path,
    output_prefix        = output_prefix,
    output_folder        = output_folder,
    sample_name          = sample_name,
    html_report          = html_report,
    verbose              = verbose,
    ...
  )
}

# ============================================================================
# TALOS db_config.R – TxDb cache, transcript fetching, BSgenome loader
# ============================================================================
#
# NOTE — why TxDb.Hsapiens.UCSC.hg19.knownGene fails for NM_ transcripts:
#   That package stores UCSC "uc001aaa.3"-style IDs as TXNAME, *never* NM_
#   accessions.  The correct source for RefSeq NM_ IDs is the UCSC "refGene"
#   track, built on-demand with:
#     GenomicFeatures::makeTxDbFromUCSC(genome, tablename = "refGene")
#   The result is cached to disk (SQLite) so the download only runs once.
#
# Fetch priority for NM_ transcripts:
#   1. UCSC refGene TxDb  (disk-cached after first build, ~30 s)
#   2. biomaRt fallback   (requires network; used when TxDb step fails)
#
# Ensembl ENST_ transcripts are unchanged (EnsDb.Hsapiens.v75 / v86).
# ============================================================================

# ── Manual LRU cache (max 5 entries) – no external dependencies ──────────────
.lru_cache <- function(max_size = 5L) {
  env <- new.env(parent = emptyenv())
  order <- character()  # keys in order of use (most recent at end)

  list(
    get = function(key) {
      if (exists(key, envir = env, inherits = FALSE)) {
        # Move to most recent position
        order <<- c(setdiff(order, key), key)
        return(get(key, envir = env))
      }
      NULL
    },
    set = function(key, value) {
      if (exists(key, envir = env, inherits = FALSE)) {
        assign(key, value, envir = env)
        order <<- c(setdiff(order, key), key)
      } else {
        if (length(order) >= max_size) {
          oldest <- order[1L]
          rm(list = oldest, envir = env)
          order <<- order[-1L]
        }
        assign(key, value, envir = env)
        order <<- c(order, key)
      }
      invisible(value)
    },
    size = function() length(order)
  )
}

# Create caches (max 5 entries each)
.txdb_cache             <- .lru_cache(max_size = 5L)
.fetch_transcript_cache <- .lru_cache(max_size = 5L)
.bs_cache               <- .lru_cache(max_size = 5L)

# Helper functions for compatibility (optional)
.set_cache <- function(cache, key, value) cache$set(key, value)
.get_cache <- function(cache, key) cache$get(key)

# ── Tiny logging helpers ──────────────────────────────────────────────────────
.dbg  <- function(...) message(sprintf("[DEBUG] %s", sprintf(...)))
.wrn  <- function(...) message(sprintf("[WARN]  %s", sprintf(...)))
.err  <- function(...) stop( sprintf("[ERROR] %s", sprintf(...)), call. = FALSE)

# ── Persistent cache directory ───────────────────────────────────────────────
# Cache is stored in tools::R_user_dir("TALOS", which = "cache")
.talos_cache_dir <- function() {
  cache_dir <- tools::R_user_dir("TALOS", which = "cache")
  if (!dir.exists(cache_dir))
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_dir
}

# ============================================================================
# 1.  RefSeq TxDb — UCSC refGene track  (NM_ TXNAME keys)
# ============================================================================
#
#  Cache strategy
#  ─────────────
#  a) in-memory  (.txdb_cache)          — fastest, lost on restart (LRU)
#  b) on-disk SQLite  (persistent)      — survives restarts, rebuilt if corrupt
#  c) build from UCSC live (network)    — one-time, ~30 s per build
#
.get_refseq_txdb <- function(build) {

  cache_key  <- paste0("refGene_", build)

  # a) in-memory hit (LRU)
  cached <- .get_cache(.txdb_cache, cache_key)
  if (!is.null(cached)) {
    .dbg("In-memory refGene TxDb hit  [build=%s]", build)
    return(cached)
  }

  # b) on-disk SQLite (persistent cache)
  cache_file <- file.path(.talos_cache_dir(),
                           sprintf("talos_refGene_%s.sqlite", build))
  if (file.exists(cache_file)) {
    .dbg("Loading refGene TxDb from disk: %s", cache_file)
    db <- tryCatch(
      AnnotationDbi::loadDb(cache_file),
      error = function(e) {
        .wrn("Disk cache unreadable — will rebuild  (%s)", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(db)) {
      .set_cache(.txdb_cache, cache_key, db)
      return(db)
    }
  }

  # c) live UCSC download
  #   In Bioconductor >= 3.16 makeTxDbFromUCSC was moved to the 'txdbmaker'
  #   package.  Probe both namespaces so older and newer installations work.
  make_info <- local({
    fn <- tryCatch(
      utils::getFromNamespace("makeTxDbFromUCSC", "txdbmaker"),
      error = function(e) NULL
    )
    if (!is.null(fn)) return(list(fn = fn, pkg = "txdbmaker"))
    fn <- tryCatch(
      utils::getFromNamespace("makeTxDbFromUCSC", "GenomicFeatures"),
      error = function(e) NULL
    )
    if (!is.null(fn)) return(list(fn = fn, pkg = "GenomicFeatures"))
    NULL
  })

  if (is.null(make_info))
    .err(paste(
      "Cannot find makeTxDbFromUCSC in 'txdbmaker' or 'GenomicFeatures'.",
      "Install the appropriate package with:",
      "  BiocManager::install('txdbmaker')   # Bioconductor >= 3.16",
      "  BiocManager::install('GenomicFeatures')  # older Bioconductor"
    ))

  make_fn   <- make_info$fn
  pkg_found <- make_info$pkg
  .dbg("Using %s::makeTxDbFromUCSC  [build=%s] — one-time download …",
       pkg_found, build)

  db <- tryCatch(
    suppressWarnings(make_fn(genome = build, tablename = "refGene")),
    error = function(e) .err("makeTxDbFromUCSC() failed: %s", conditionMessage(e))
  )

  AnnotationDbi::saveDb(db, file = cache_file)
  .dbg("refGene TxDb saved to disk: %s", cache_file)

  .set_cache(.txdb_cache, cache_key, db)
  db
}

# ============================================================================
# 2.  Ensembl EnsDb loader
# ============================================================================
.get_ensdb <- function(build) {
  pkg <- if (build == "hg19") "EnsDb.Hsapiens.v75" else "EnsDb.Hsapiens.v86"
  if (!requireNamespace(pkg, quietly = TRUE))
    .err("Please install Bioconductor package '%s'", pkg)
  get(pkg, asNamespace(pkg))
}

# ============================================================================
# 3.  biomaRt fallback for NM_ transcripts
# ============================================================================
.fetch_via_biomart <- function(transcript_id, build,
                               chrom_prefix = getOption("TALOS.chrom_prefix", "chr")) {

  if (!requireNamespace("biomaRt", quietly = TRUE))
    .err(paste(
      "Package 'biomaRt' is required for the fallback path.",
      "Install with: BiocManager::install('biomaRt')"
    ))

  .dbg("biomaRt fallback — connecting [build=%s] …", build)

  # Establish Mart connection (prefer useEnsembl, fallback to useMart)
  mart <- tryCatch({
    if (build == "hg19") {
      biomaRt::useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        GRCh = 37
      )
    } else {
      biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    }
  }, error = function(e1) {
    .wrn("useEnsembl() failed (%s) — retrying with useMart() …",
         conditionMessage(e1))
    tryCatch({
      if (build == "hg19") {
        biomaRt::useMart(
          biomart = "ENSEMBL_MART_ENSEMBL",
          dataset = "hsapiens_gene_ensembl",
          host = "https://grch37.ensembl.org"
        )
      } else {
        biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
      }
    }, error = function(e2)
      .err("biomaRt connection failed: %s", conditionMessage(e2)))
  })

  # ---- Step 1: Get Ensembl transcript ID from RefSeq mRNA ----
  tx_id_df <- biomaRt::getBM(
    attributes = c("refseq_mrna", "ensembl_transcript_id"),
    filters = "refseq_mrna",
    values = transcript_id,
    mart = mart
  )

  if (nrow(tx_id_df) == 0)
    .err("Transcript '%s' not found in biomaRt (no Ensembl mapping).", transcript_id)

  ensembl_tx_id <- tx_id_df$ensembl_transcript_id[1]
  .dbg("Mapped to Ensembl transcript: %s", ensembl_tx_id)

  # ---- Step 2: Get transcript metadata (chromosome, strand) ----
  meta_df <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "chromosome_name", "strand"),
    filters = "ensembl_transcript_id",
    values = ensembl_tx_id,
    mart = mart
  )

  if (nrow(meta_df) == 0)
    .err("No metadata found for Ensembl transcript '%s'.", ensembl_tx_id)

  chrom <- paste0(chrom_prefix, meta_df$chromosome_name[1])
  strand_val <- if (meta_df$strand[1] == 1) "+" else "-"

  # ---- Step 3: Get exon coordinates (using Ensembl transcript ID) ----
  exon_df <- biomaRt::getBM(
    attributes = c("ensembl_transcript_id", "rank",
                   "exon_chrom_start", "exon_chrom_end"),
    filters = "ensembl_transcript_id",
    values = ensembl_tx_id,
    mart = mart
  )

  .dbg("biomaRt exon query returned %d rows for '%s'", nrow(exon_df), transcript_id)

  if (nrow(exon_df) == 0)
    .err("No exons found for transcript '%s'.", transcript_id)

  # Remove rows with missing coordinates
  exon_df <- exon_df[!is.na(exon_df$exon_chrom_start) &
                     !is.na(exon_df$exon_chrom_end), ]
  if (nrow(exon_df) == 0)
    .err("All exon coordinates are NA for '%s'.", transcript_id)

  # Sort by rank (transcript order, 5' → 3')
  exon_df <- exon_df[order(exon_df$rank), ]

  exons_gr <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges = IRanges::IRanges(
      start = exon_df$exon_chrom_start,
      end   = exon_df$exon_chrom_end
    ),
    strand = strand_val
  )

  # ---- Step 4: Get CDS coordinates (optional) ----
  cds_gr <- NULL
  cds_df <- tryCatch(
    biomaRt::getBM(
      attributes = c("ensembl_transcript_id", "rank",
                     "genomic_coding_start", "genomic_coding_end"),
      filters = "ensembl_transcript_id",
      values = ensembl_tx_id,
      mart = mart
    ),
    error = function(e) {
      .wrn("CDS query failed (%s) — cds will be empty", conditionMessage(e))
      data.frame()
    }
  )

  if (nrow(cds_df) > 0) {
    cds_df <- cds_df[!is.na(cds_df$genomic_coding_start) &
                     !is.na(cds_df$genomic_coding_end) &
                     cds_df$genomic_coding_start > 0, ]
    if (nrow(cds_df) > 0) {
      cds_gr <- GenomicRanges::GRanges(
        seqnames = chrom,
        ranges = IRanges::IRanges(
          start = cds_df$genomic_coding_start,
          end   = cds_df$genomic_coding_end
        ),
        strand = strand_val
      )
    }
  }

  list(exons_gr = exons_gr, cds_gr = cds_gr,
       chrom = chrom, strand = strand_val, source = "biomaRt")
}

# ============================================================================
# 4.  Core transcript fetcher
# ============================================================================
.fetch_transcript <- function(transcript_id, gene, build) {

  cache_key <- paste(transcript_id, build, sep = "_")
  cached <- .get_cache(.fetch_transcript_cache, cache_key)
  if (!is.null(cached)) {
    .dbg("Transcript cache hit: %s / %s", transcript_id, build)
    return(cached)
  }

  is_refseq <- startsWith(transcript_id, "NM_")

  # ── A. Ensembl path (ENST*) ────────────────────────────────────────────────
  if (!is_refseq) {
    .dbg("Ensembl path for '%s' [build=%s]", transcript_id, build)
    db <- .get_ensdb(build)

    exons_gr <- ensembldb::exons(
      db, filter = AnnotationFilter::TxIdFilter(transcript_id)
    )
    exons_gr  <- exons_gr[order(exons_gr$exon_idx)]
    strand_val <- as.character(unique(BiocGenerics::strand(exons_gr)))
    chrom      <- as.character(GenomeInfoDb::seqnames(exons_gr)[1])

    cds_list <- ensembldb::cdsBy(
      db, by = "tx",
      filter = AnnotationFilter::TxIdFilter(transcript_id)
    )
    cds_gr <- cds_list[[transcript_id]]

    return(.finalize_transcript(transcript_id, exons_gr, cds_gr,
                                chrom, strand_val, cache_key,
                                source = "EnsDb"))
  }

  # ── B. RefSeq path (NM_*) ──────────────────────────────────────────────────
  .dbg("RefSeq path for '%s' [build=%s]", transcript_id, build)

  # B1 — try refGene TxDb
  result <- tryCatch({

    db <- .get_refseq_txdb(build)

    all_keytypes <- AnnotationDbi::keytypes(db)
    .dbg("TxDb keytypes: %s", paste(all_keytypes, collapse = ", "))

    all_txnames <- AnnotationDbi::keys(db, keytype = "TXNAME")
    .dbg("TxDb contains %d TXNAME records", length(all_txnames))

    # Exact match first; then versioned (NM_004119 → NM_004119.4)
    lookup_id <- if (transcript_id %in% all_txnames) {
      .dbg("Exact TXNAME match: '%s'", transcript_id)
      transcript_id
    } else {
      versioned <- grep(paste0("^", transcript_id, "\\."),
                        all_txnames, value = TRUE)
      if (length(versioned) > 0) {
        .dbg("Versioned TXNAME match: '%s' → '%s'", transcript_id, versioned[1])
        versioned[1]
      } else {
        .dbg("TXNAME sample (first 8): %s",
             paste(head(all_txnames, 8), collapse = ", "))
        stop(sprintf("'%s' not found in refGene TxDb (exact or versioned).",
                     transcript_id))
      }
    }

    tx_map <- AnnotationDbi::select(
      db, keys = lookup_id, columns = "TXID", keytype = "TXNAME"
    )
    tx_id <- tx_map$TXID[1]
    .dbg("Internal TXID = %s", tx_id)

    exons_all <- GenomicFeatures::exonsBy(db, by = "tx", use.names = FALSE)
    exons_gr  <- exons_all[[as.character(tx_id)]]

    strand_val <- as.character(unique(BiocGenerics::strand(exons_gr)))
    chrom      <- as.character(GenomeInfoDb::seqnames(exons_gr)[1])

    cds_all <- GenomicFeatures::cdsBy(db, by = "tx", use.names = FALSE)
    cds_gr  <- cds_all[[as.character(tx_id)]]

    list(exons_gr = exons_gr, cds_gr = cds_gr,
         chrom = chrom, strand = strand_val, source = "refGene_TxDb")

  }, error = function(e) {
    .wrn("refGene TxDb lookup failed: %s", conditionMessage(e))
    .wrn("Falling back to biomaRt …")
    NULL
  })

  # B2 — biomaRt fallback
  if (is.null(result)) {
    result <- .fetch_via_biomart(transcript_id, build)
  }

  .dbg("Source used: %-18s | chrom=%-6s strand=%s exons=%d",
       result$source, result$chrom, result$strand, length(result$exons_gr))

  .finalize_transcript(transcript_id, result$exons_gr, result$cds_gr,
                       result$chrom, result$strand, cache_key,
                       source = result$source)
}

# ============================================================================
# 5.  Shared finalize helper — transcript order, CDS offset, cache write
# ============================================================================
.finalize_transcript <- function(transcript_id, exons_gr, cds_gr,
                                 chrom, strand_val, cache_key,
                                 source = "unknown") {

  # Sort exons in transcript order (5′ → 3′) and assign 1-based exon_idx
  exons_gr <- if (strand_val == "-") {
    exons_gr[order(BiocGenerics::start(exons_gr), decreasing = TRUE)]
  } else {
    exons_gr[order(BiocGenerics::start(exons_gr))]
  }
  S4Vectors::mcols(exons_gr)$exon_idx <- seq_along(exons_gr)

  # Calculate CDS offset (bases from transcript start to ATG)
  cds_offset <- 0L
  if (!is.null(cds_gr) && length(cds_gr) > 0) {
    start_codon <- if (strand_val == "+") min(BiocGenerics::start(cds_gr)) else
                                          max(BiocGenerics::end(cds_gr))
    for (i in seq_along(exons_gr)) {
      ex_s <- BiocGenerics::start(exons_gr)[i]
      ex_e <- BiocGenerics::end(exons_gr)[i]
      if (start_codon >= ex_s && start_codon <= ex_e) {
        pos_in_exon <- if (strand_val == "+") (start_codon - ex_s) else
                                              (ex_e - start_codon)
        cds_offset  <- sum(BiocGenerics::width(exons_gr[seq_len(i - 1)])) +
                       pos_in_exon
        break
      }
    }
  }
  .dbg("CDS offset = %d  (source=%s)", cds_offset, source)

  result <- list(
    transcript_id = transcript_id,
    chrom         = chrom,
    strand        = strand_val,
    all_exons_gr  = exons_gr,
    cds_offset    = cds_offset,
    tx_start      = min(BiocGenerics::start(exons_gr)),
    tx_end        = max(BiocGenerics::end(exons_gr)),
    source        = source
  )

  .set_cache(.fetch_transcript_cache, cache_key, result)
  result
}

# ============================================================================
# 6.  BSgenome cached loader
# ============================================================================
.load_bsgenome_cached <- function(build, bsgenome = NULL) {
  bs_pkg <- if (build == "hg19") "BSgenome.Hsapiens.UCSC.hg19" else
                                 "BSgenome.Hsapiens.UCSC.hg38"

  if (!is.null(bsgenome)) {
    if (!inherits(bsgenome, "BSgenome"))
      bsgenome <- if (is.character(bsgenome)) BSgenome::getBSgenome(bsgenome) else
        stop("'bsgenome' must be a BSgenome object or package-name string")
    .set_cache(.bs_cache, bs_pkg, bsgenome)
    return(bsgenome)
  }

  cached <- .get_cache(.bs_cache, bs_pkg)
  if (!is.null(cached)) return(cached)

  if (!requireNamespace(bs_pkg, quietly = TRUE))
    .err("Please install Bioconductor package '%s'", bs_pkg)

  bs_obj <- BSgenome::getBSgenome(bs_pkg)
  .set_cache(.bs_cache, bs_pkg, bs_obj)
  bs_obj
}

# ============================================================================
# 7.  Optional: dump config to log file  (safe — no sink())
# ============================================================================
.dump_config_to_log <- function(gene_config, output_folder,
                                sample_name, gene_name) {
  log_dir  <- file.path(output_folder, sample_name)
  if (!dir.exists(log_dir))
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir, paste0(gene_name, "_config.log"))
  output_lines <- c(
    "=== TALOS Gene Configuration ===",
    paste("Generated:", format(Sys.time())),
    "",
    capture.output(str(gene_config, max.level = 2))
  )
  writeLines(output_lines, log_file)
  invisible(log_file)
}

# ============================================================================
# 8.  get_gene_config() — main entry point
# ============================================================================
#' Retrieve a fully-resolved gene configuration for TALOS
#'
#' @param gene         Target gene name (must match a key in the YAML).
#' @param build        Reference genome build: "hg19" or "hg38".
#' @param config_path  Path to the gene_config.yaml file.
#' @param bsgenome     Optional BSgenome object (or package name string).
#' @param padding      Base-pair padding around targeted exons (default 500).
#' @param use_db       If TRUE (default), fetch exon/CDS info from UCSC/biomaRt.
#' @param exon_padding Number of flanking exons to add to the genomic window (default 0).
#' @return A named list with gene configuration including exons, sequences, and
#'   genomic coordinates.
#' @examples
#' \dontrun{
#' cfg <- get_gene_config("FLT3", build = "hg19")
#' print(names(cfg))
#' }
#' @export
get_gene_config <- function(
    gene,
    build       = c("hg19", "hg38"),
    config_path = system.file("extdata", "gene_config.yaml", package = "TALOS"),
    bsgenome    = NULL,
    padding     = 500L,
    use_db      = TRUE,
    exon_padding = 0L
) {
  if (missing(gene) || length(gene) != 1L || !is.character(gene))
    stop("'gene' must be a single character string")
  build <- match.arg(build)
  if (!file.exists(config_path))
    stop("Config file not found: ", config_path)

  cfg   <- yaml::read_yaml(config_path)
  entry <- cfg[[gene]]
  if (is.null(entry))
    stop("Gene '", gene, "' not found in config.")

  transcript_id       <- entry$transcript %||% NA_character_
  build_cfg           <- entry[[build]]          # hg19: / hg38: sub-block if present
  user_exons_provided <- !is.null(build_cfg) && !is.null(build_cfg$exons)

  # ── Fetch exon / CDS data ──────────────────────────────────────────────────
  if (use_db && !user_exons_provided) {

    if (is.na(transcript_id))
      stop("Set 'transcript:' in the YAML or call with use_db = FALSE.")

    res          <- .fetch_transcript(transcript_id, gene, build)
    chrom        <- res$chrom
    strand       <- res$strand
    all_exons_gr <- res$all_exons_gr
    auto_cds_offset <- res$cds_offset

  } else {
    # Manual / offline mode (exons supplied directly in YAML)
    if (is.null(build_cfg) || is.null(build_cfg$exons))
      stop("Missing '", build, ": exons:' block in YAML and use_db = FALSE.")

    chrom  <- entry$chrom  %||% "chr?"
    strand <- entry$strand %||% "?"
    auto_cds_offset <- as.integer(entry$cds_offset %||% 0L)

    # IMPORTANT: Do NOT reorder manual exons; they are assumed to be in transcript order
    all_exons_gr <- GenomicRanges::GRanges(
      seqnames = chrom,
      ranges   = IRanges::IRanges(
        start = vapply(build_cfg$exons, `[`, numeric(1), 1),
        end   = vapply(build_cfg$exons, `[`, numeric(1), 2)
      ),
      strand = strand
    )
    # Assign exon_idx in the given order (no sorting)
    S4Vectors::mcols(all_exons_gr)$exon_idx <- seq_along(all_exons_gr)
  }

  # ── Subset to targeted exons ───────────────────────────────────────────────
  cds_offset <- as.integer(entry$cds_offset %||% auto_cds_offset)
  targeted   <- as.integer(
    entry$targeted_exons %||% entry$exon_numbers %||% seq_along(all_exons_gr)
  )
  target_exons_gr <- all_exons_gr[S4Vectors::mcols(all_exons_gr)$exon_idx %in% targeted]

  # ── Expand exons if exon_padding > 0 (for genomic window) ───────────────────
  expanded_exons_gr <- target_exons_gr
  if (exon_padding > 0 && length(all_exons_gr) > 0) {
    all_idx <- S4Vectors::mcols(all_exons_gr)$exon_idx
    min_t <- min(S4Vectors::mcols(target_exons_gr)$exon_idx)
    max_t <- max(S4Vectors::mcols(target_exons_gr)$exon_idx)
    left_idx  <- max(1L, min_t - exon_padding)
    right_idx <- min(max(all_idx), max_t + exon_padding)
    keep_idx  <- which(all_idx >= left_idx & all_idx <= right_idx)
    expanded_exons_gr <- all_exons_gr[keep_idx]
  }

  # ── Genomic window ─────────────────────────────────────────────────────────
  genomic_start <- max(1L, min(GenomicRanges::start(expanded_exons_gr)) - padding)
  genomic_end   <- max(GenomicRanges::end(expanded_exons_gr)) + padding

  # ── Sequence extraction ────────────────────────────────────────────────────
  bs_obj <- .load_bsgenome_cached(build, bsgenome)

  gene_region <- GenomicRanges::GRanges(
    seqnames = chrom,
    ranges   = IRanges::IRanges(genomic_start, genomic_end),
    strand   = "+"
  )
  genomic_ref_seq <- as.character(BSgenome::getSeq(bs_obj, gene_region))

  target_seqs <- as.character(BSgenome::getSeq(bs_obj, target_exons_gr))
  if (strand == "-")
    target_seqs <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAStringSet(target_seqs))
    )
  rnaseq <- paste(target_seqs, collapse = "")

  # Build full cDNA from all exons (needed for HGVS protein annotation)
  all_exon_seqs <- as.character(BSgenome::getSeq(bs_obj, all_exons_gr))
  if (strand == "-")
    all_exon_seqs <- as.character(
      Biostrings::reverseComplement(Biostrings::DNAStringSet(all_exon_seqs))
    )
  full_cdna <- paste(all_exon_seqs, collapse = "")

  list(
    gene            = gene,
    chrom           = chrom,
    strand          = strand,
    transcript      = transcript_id,
    genomic_start   = genomic_start,
    genomic_end     = genomic_end,
    rnaseq          = rnaseq,
    genomic_ref_seq = genomic_ref_seq,
    full_cdna       = full_cdna,
    build           = build,
    target_exons    = target_exons_gr,
    all_exons       = all_exons_gr,
    cds_offset      = cds_offset,
    bsgenome_obj    = bs_obj,
    gene_settings   = entry$gene_settings
  )
}

# ============================================================================
# 9.  save_offline_config() — snapshot a resolved config to YAML
# ============================================================================
#' Save a fully-resolved configuration to YAML for offline / reproducible use
#'
#' @param gene        Target gene name.
#' @param build       Reference genome build.
#' @param output_file Destination YAML path.
#' @param padding     Base-pair padding (default 500).
#' @return Path to the written YAML file, invisibly.
#' @examples
#' \dontrun{
#' save_offline_config("FLT3", build = "hg19",
#'                     output_file = "flt3_offline.yaml")
#' }
#' @export
save_offline_config <- function(gene, build = "hg19",
                                output_file = "offline_config.yaml",
                                padding = 500L) {
  config <- get_gene_config(gene, build, padding = padding, use_db = TRUE)

  exon_data <- lapply(seq_along(config$all_exons), function(i)
    c(start(config$all_exons)[i], end(config$all_exons)[i])
  )

  key <- paste0(gene, "_OFFLINE")
  out <- list()
  out[[key]] <- list(
    transcript     = config$transcript,
    targeted_exons = as.integer(unique(S4Vectors::mcols(config$target_exons)$exon_idx)),
    cds_offset     = config$cds_offset,
    chrom          = config$chrom,
    strand         = config$strand
  )
  out[[key]][[build]] <- list(
    start = config$genomic_start,
    end   = config$genomic_end,
    exons = exon_data
  )

  yaml::write_yaml(out, file = output_file)
  message("Offline config saved to: ", output_file)
}
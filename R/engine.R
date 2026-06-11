# ============================================================================
# TALOS – K-mer helpers and optional-use capability flags
#
# Sources in order:
#   engine_bam.R        – BAM loading & CIGAR pre-filter
#   engine_candidates.R – candidate extraction & clustering
#   engine_metrics.R    – variant metrics & filter application
# ============================================================================

.has_fastmatch  <- requireNamespace("fastmatch",  quietly = TRUE)
.has_cigarillo  <- requireNamespace("cigarillo",  quietly = TRUE)
.has_pwalign    <- requireNamespace("pwalign",    quietly = TRUE)
.has_biostrings <- requireNamespace("Biostrings", quietly = TRUE)

# ---------------------------------------------------------------------------
# Prepare reference k-mer table from a DNA string
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Hash-accelerated k-mer lookup (fastmatch when available)
# ---------------------------------------------------------------------------
.kmer_match <- function(query, table) {
  if (.has_fastmatch) fastmatch::fmatch(query, table) else match(query, table)
}
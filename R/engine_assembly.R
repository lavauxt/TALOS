#' Kernel Density Estimation for Breakpoint Refinement
#' @param positions Numeric vector of breakpoint positions
#' @return Refined breakpoint (mode of density)
.kde_breakpoint <- function(positions) {
  positions <- positions[!is.na(positions)]
  if (length(positions) == 0) return(NA_integer_)
  # density() needs a non-zero bandwidth: identical positions are the answer.
  if (length(unique(positions)) == 1L) return(as.integer(positions[1]))
  d <- stats::density(positions, bw = "nrd")
  as.integer(d$x[which.max(d$y)])
}

#' Build Local de Bruijn Graph from Soft-Clip Sequences
#'
#' @param seqs Character vector of soft-clip sequences
#' @param k K-mer size (default 31)
#' @param min_coverage Minimum k-mer coverage to keep
#' @param max_steps Hard cap on walk length (safety bound against cycles)
#' @return Assembled sequence, or NA if assembly isn't possible
.build_local_debruijn <- function(seqs, k = 31, min_coverage = 2, max_steps = 2000L) {
  if (length(seqs) == 0) return(NA_character_)
  all_seqs <- seqs[!is.na(seqs)]
  # Fallbacks get every input sequence: the k-filtered set was passed before,
  # which is empty on this branch and always returned NA.
  fallback <- function() .simple_overlap_assembly(all_seqs, min_overlap = min(k, 10L))
  seqs <- all_seqs[nchar(all_seqs) >= k]
  if (length(seqs) == 0) return(fallback())
  kmer_env <- new.env(hash = TRUE, parent = emptyenv())
  for (s in seqs) {
    n_k <- nchar(s) - k + 1L
    if (n_k < 1L) next
    for (i in seq_len(n_k)) {
      kmer <- substr(s, i, i + k - 1L)
      cur  <- kmer_env[[kmer]]
      kmer_env[[kmer]] <- if (is.null(cur)) 1L else cur + 1L
    }
  }
  kmer_names <- ls(kmer_env)
  if (length(kmer_names) == 0L) return(fallback())
  coverage <- vapply(kmer_names, function(kn) kmer_env[[kn]], integer(1L))
  names(coverage) <- kmer_names
  kmer_names <- kmer_names[coverage >= min_coverage]
  if (length(kmer_names) == 0L) return(fallback())
  coverage <- coverage[kmer_names]
  prefix_idx <- new.env(hash = TRUE, parent = emptyenv())
  for (kn in kmer_names) {
    px <- substr(kn, 1L, k - 1L)
    bucket <- prefix_idx[[px]]
    prefix_idx[[px]] <- if (is.null(bucket)) kn else c(bucket, kn)
  }
  successors_of <- function(kmer) prefix_idx[[substr(kmer, 2L, k)]]
  indeg_env <- new.env(hash = TRUE, parent = emptyenv())
  for (kn in kmer_names) {
    for (nx in successors_of(kn)) {
      if (identical(nx, kn)) next
      cur <- indeg_env[[nx]]
      indeg_env[[nx]] <- if (is.null(cur)) 1L else cur + 1L
    }
  }
  indeg <- vapply(kmer_names, function(kn) {
    v <- indeg_env[[kn]]; if (is.null(v)) 0L else v
  }, integer(1L))
  names(indeg) <- kmer_names
  zero_indeg <- kmer_names[indeg == 0L]
  seed <- if (length(zero_indeg) > 0L) zero_indeg[which.max(coverage[zero_indeg])]
  else kmer_names[which.max(coverage)]
  path <- seed; visited <- seed; current <- seed
  for (step in seq_len(max_steps)) {
    cand <- setdiff(successors_of(current), visited)
    if (length(cand) == 0L) break
    nxt <- cand[which.max(coverage[cand])]
    path <- c(path, nxt); visited <- c(visited, nxt); current <- nxt
  }
  if (length(path) < 2L) return(path[1])
  paste0(path[1], paste(substr(path[-1], k, k), collapse = ""))
}

#' Simple Overlap Assembly (fallback for short/sparse soft-clip sets)
#' @param seqs Character vector of sequences
#' @param min_overlap Minimum overlap length to merge two sequences
#' @return Assembled contig, or NA
.simple_overlap_assembly <- function(seqs, min_overlap = 10) {
  if (length(seqs) == 0) return(NA_character_)
  seqs <- seqs[!is.na(seqs) & nchar(seqs) >= min_overlap]
  if (length(seqs) == 0) return(NA_character_)
  contig <- seqs[which.max(nchar(seqs))]
  seqs <- seqs[-which.max(nchar(seqs))]
  changed <- TRUE
  while (changed && length(seqs) > 0) {
    changed <- FALSE
    for (i in seq_along(seqs)) {
      s <- seqs[i]
      if (grepl(s, contig, fixed = TRUE)) {
        seqs <- seqs[-i]; changed <- TRUE; break
      }
      for (ov in seq(min_overlap, nchar(s), by = 1)) {
        if (substr(s, 1, ov) == substr(contig, nchar(contig) - ov + 1, nchar(contig))) {
          contig <- paste0(contig, substr(s, ov + 1, nchar(s)))
          seqs <- seqs[-i]; changed <- TRUE; break
        }
        if (substr(contig, 1, ov) == substr(s, nchar(s) - ov + 1, nchar(s))) {
          contig <- paste0(substr(s, 1, nchar(s) - ov), contig)
          seqs <- seqs[-i]; changed <- TRUE; break
        }
      }
      if (changed) break
    }
  }
  contig
}

#' [BUG 6 FIX] Weighted Consensus from Soft-Clip Sequences
#'
#' Per-position weighted majority vote. Weights (e.g. MAPQ) are now actually
#' used: each base vote is weighted by its read's weight. Falls back to
#' unweighted voting when weights are NULL, malformed, or all zero.
#'
#' @param seqs Character vector of sequences
#' @param weights Numeric vector of weights, one per sequence (optional)
#' @return Consensus sequence
.assemble_weighted_consensus <- function(seqs, weights = NULL) {
  if (length(seqs) == 0L) return(NA_character_)
  if (is.null(weights) || length(weights) != length(seqs)) {
    weights <- rep(1, length(seqs))
  }
  weights <- suppressWarnings(as.numeric(weights))
  weights[is.na(weights) | !is.finite(weights) | weights < 0] <- 0
  if (sum(weights) <= 0) weights <- rep(1, length(seqs))
  max_len <- max(nchar(seqs), na.rm = TRUE)
  if (max_len == 0L) return(NA_character_)
  mat <- matrix(NA_character_, nrow = length(seqs), ncol = max_len)
  for (i in seq_along(seqs)) {
    chars <- strsplit(seqs[i], "")[[1]]
    mat[i, seq_len(length(chars))] <- chars
  }
  consensus <- vapply(seq_len(max_len), function(pos) {
    ok  <- !is.na(mat[, pos])
    col <- mat[ok, pos]
    if (length(col) == 0L) return("N")
    w <- weights[ok]
    score <- tapply(w, col, sum)
    names(score)[which.max(score)]
  }, character(1L))
  paste(consensus, collapse = "")
}

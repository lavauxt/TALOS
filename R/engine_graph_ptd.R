# Consensus coordinate: mode of precise (split-read) positions, otherwise the
# median of imprecise (discordant-pair) positions.
.consensus_coord <- function(x, precise) {
  xp <- x[precise & !is.na(x)]
  if (length(xp)) return(as.integer(.mode(xp)))
  as.integer(round(stats::median(x, na.rm = TRUE)))
}

#' Map tandem-duplication evidence onto annotated exon-graph edges
#'
#' Each event is a duplicated interval [S, E]. The duplicated exons are those
#' lying inside [S, E] (within \code{exon_tol} bp); the back-edge runs from
#' the last duplicated exon to the first one, in transcript order, so it is
#' strand independent. The previous "nearest exon to each coordinate" rule
#' mis-assigned intronic breakpoints (nearest exon of a breakpoint in intron
#' k/k+1 can be either exon) and dropped single-exon duplications as
#' self-loops.
#'
#' @param evidence_dt data.frame(qname, start, end, prob, evidence_type, precise)
#' @param exon_granges Exons (GRanges with exon_idx, or data.frame start/end/exon_idx)
#' @param max_distance Maximum distance (bp) between each breakpoint and the
#'   nearest exon; events further away are ignored. KMT2A introns run into
#'   several kb, hence the default.
#' @param exon_tol Slack (bp) when testing that an exon lies inside [S, E].
#' @return Aggregated edge data.frame (one row per (from, to) exon pair) with
#'   joint probability, unique-read support and consensus S / E.
.map_evidence_to_exon_edges <- function(evidence_dt, exon_granges, max_distance = 5000L,
                                        exon_tol = 10L) {
  ex <- .exons_df(exon_granges)
  if (is.null(evidence_dt) || nrow(evidence_dt) == 0L || is.null(ex)) return(data.frame())

  near <- function(p) vapply(p, function(v) {
    min(pmax(0L, ex$start - v, v - ex$end))
  }, numeric(1L))
  ev <- evidence_dt[near(evidence_dt$start) <= max_distance &
                      near(evidence_dt$end) <= max_distance, , drop = FALSE]
  if (nrow(ev) == 0L) return(data.frame())

  inside <- lapply(seq_len(nrow(ev)), function(i) {
    ex$exon_idx[ex$start >= ev$start[i] - exon_tol & ex$end <= ev$end[i] + exon_tol]
  })
  has_exon <- lengths(inside) > 0L
  if (!any(has_exon)) return(data.frame())
  ev <- ev[has_exon, , drop = FALSE]
  inside <- inside[has_exon]
  ev$from <- paste0("Exon_", vapply(inside, max, numeric(1L)))
  ev$to <- paste0("Exon_", vapply(inside, min, numeric(1L)))

  key <- paste(ev$from, ev$to, sep = "->")
  agg <- lapply(split(seq_len(nrow(ev)), key), function(rows) {
    sub <- ev[rows, , drop = FALSE]
    best <- tapply(sub$prob, sub$qname, max)
    data.frame(
      from = sub$from[1L],
      to = sub$to[1L],
      prob = 1 - prod(1 - best),
      support_count = length(best),
      n_split = length(unique(sub$qname[sub$evidence_type == "split_read"])),
      n_pairs = length(unique(sub$qname[sub$evidence_type == "discordant_pair"])),
      qnames = paste(names(best), collapse = ";"),
      evidence_type = paste(unique(sub$evidence_type), collapse = "+"),
      dup_start = .consensus_coord(sub$start, sub$precise),
      dup_end = .consensus_coord(sub$end, sub$precise),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, agg)
  rownames(out) <- NULL
  out
}

#' Minimal Dijkstra shortest path over a small named-node edge list
#'
#' Exon graphs here are tiny (dozens of nodes at most), so a plain O(V^2)
#' Dijkstra in base R is fast enough and avoids an igraph dependency.
#'
#' @param edges data.frame(from, to, weight)
#' @param nodes character vector of all node names
#' @param from,to start/end node names
#' @return character vector of node names on the shortest path (empty if unreachable)
.shortest_path_dijkstra <- function(edges, nodes, from, to) {
  n <- length(nodes)
  idx <- stats::setNames(seq_len(n), nodes)
  if (!from %in% nodes || !to %in% nodes) return(character(0))
  adj <- vector("list", n)
  for (i in seq_len(nrow(edges))) {
    if (!edges$from[i] %in% nodes || !edges$to[i] %in% nodes) next
    a <- idx[[edges$from[i]]]
    adj[[a]] <- rbind(adj[[a]], c(idx[[edges$to[i]]], edges$weight[i]))
  }
  dist <- rep(Inf, n); prev <- rep(NA_integer_, n); visited <- rep(FALSE, n)
  dist[idx[[from]]] <- 0
  for (iter in seq_len(n)) {
    remaining <- which(!visited)
    if (length(remaining) == 0L) break
    u <- remaining[which.min(dist[remaining])]
    if (!is.finite(dist[u])) break
    visited[u] <- TRUE
    if (nodes[u] == to) break
    if (!is.null(adj[[u]])) {
      mat <- matrix(adj[[u]], ncol = 2)
      for (r in seq_len(nrow(mat))) {
        v <- mat[r, 1]; w <- mat[r, 2]
        if (dist[u] + w < dist[v]) { dist[v] <- dist[u] + w; prev[v] <- u }
      }
    }
  }
  ti <- idx[[to]]
  if (!is.finite(dist[ti])) return(character(0))
  path <- ti; cur <- ti
  while (!is.na(prev[cur])) { cur <- prev[cur]; path <- c(cur, path) }
  nodes[path]
}

#' Identify PTD back-edges and score their duplication cycles
#'
#' [BUG 7] The score now uses the back-edge itself: \code{nll_score} is the
#' cycle cost (forward path + back-edge, each \code{-log(p)}) and
#' \code{evidence_score = log1p(support) - nll_score} is the sortable
#' quantity. Single-exon duplications (back-edge onto the same exon) are
#' valid cycles. Results are sorted by evidence_score, best first.
#'
#' @param edges data.frame(from, to, prob, weight, support_count, qnames,
#'   evidence_type, ...)
#' @param nodes character vector of exon node names, in transcript order
#' @return list of PTD candidates, or NULL
.detect_ptd_cycles <- function(edges, nodes) {
  if (is.null(edges) || nrow(edges) == 0L) return(NULL)
  from_idx <- as.numeric(sub("Exon_", "", edges$from))
  to_idx <- as.numeric(sub("Exon_", "", edges$to))
  back_mask <- (from_idx >= to_idx) & (edges$evidence_type != "canonical")
  back_edges <- edges[back_mask, , drop = FALSE]
  if (nrow(back_edges) == 0L) return(NULL)
  forward <- edges[!back_mask, , drop = FALSE]
  edge_pair_key <- paste(forward$from, forward$to, sep = ",")
  col_or_na <- function(df, col, i) if (col %in% names(df)) df[[col]][i] else NA

  results <- list()
  for (i in seq_len(nrow(back_edges))) {
    start_node <- back_edges$to[i]
    end_node <- back_edges$from[i]
    node_sequence <- .shortest_path_dijkstra(forward, nodes, start_node, end_node)
    if (length(node_sequence) == 0L) next
    path_pairs <- if (length(node_sequence) > 1L) {
      paste(node_sequence[-length(node_sequence)], node_sequence[-1L], sep = ",")
    } else {
      character(0)
    }
    path_w <- sum(forward$weight[edge_pair_key %in% path_pairs])
    back_w <- -log(back_edges$prob[i] + 1e-9)
    support <- back_edges$support_count[i]
    dup_start <- col_or_na(back_edges, "dup_start", i)
    dup_end <- col_or_na(back_edges, "dup_end", i)
    results[[length(results) + 1L]] <- list(
      ptd_signature = paste(end_node, "->", start_node),
      dup_start_exon = start_node,
      dup_end_exon = end_node,
      internal_path = node_sequence,
      nll_score = path_w + back_w,
      evidence_score = log1p(support) - (path_w + back_w),
      support_reads = support,
      support_qnames = back_edges$qnames[i],
      evidence_type = back_edges$evidence_type[i],
      n_split = col_or_na(back_edges, "n_split", i),
      n_pairs = col_or_na(back_edges, "n_pairs", i),
      dup_start = dup_start,
      dup_end = dup_end,
      dup_length = if (is.na(dup_start) || is.na(dup_end)) NA_integer_ else as.integer(dup_end - dup_start + 1L)
    )
  }
  if (length(results) == 0L) return(NULL)
  scores <- vapply(results, function(r) r$evidence_score, numeric(1L))
  results[order(scores, decreasing = TRUE)]
}

#' Build the exon graph from DUP evidence and return scored PTD cycles
#'
#' Pure (no BAM access) so it can be unit-tested; see extract_graph_candidates().
#' @param evidence data.frame from .parse_sa_tags() / .extract_discordant_and_shadow()
#' @param exons Exons (GRanges or data.frame start/end/exon_idx)
#' @param max_exon_distance,min_support See extract_graph_candidates().
#' @return list of PTD cycles (possibly empty)
.graph_cycles_from_evidence <- function(evidence, exons, max_exon_distance = 5000L,
                                        min_support = 2L) {
  ex <- .exons_df(exons)
  if (is.null(ex) || nrow(ex) < 1L || is.null(evidence) || nrow(evidence) == 0L) return(list())
  ptd_edges <- .map_evidence_to_exon_edges(evidence, ex, max_distance = max_exon_distance)
  if (nrow(ptd_edges) == 0L) return(list())
  exon_names <- paste0("Exon_", ex$exon_idx)
  canonical <- if (nrow(ex) >= 2L) {
    data.frame(from = exon_names[-length(exon_names)], to = exon_names[-1L],
               prob = 0.99, support_count = 100L, qnames = "", evidence_type = "canonical",
               stringsAsFactors = FALSE)
  } else NULL
  common <- c("from", "to", "prob", "support_count", "qnames", "evidence_type")
  extra <- setdiff(names(ptd_edges), common)
  if (!is.null(canonical)) for (cl in extra) canonical[[cl]] <- NA
  combined <- rbind(canonical, ptd_edges[, c(common, extra), drop = FALSE])
  combined$weight <- -log(combined$prob + 1e-9)
  cycles <- .detect_ptd_cycles(combined, exon_names)
  if (is.null(cycles)) return(list())
  # [BUG 7 gate] keep cycles whose evidence outweighs their cost and that
  # have enough distinct reads.
  Filter(function(r) r$evidence_score > 0 && r$support_reads >= min_support, cycles)
}

#' Main Graph Candidate Extraction Workflow
#'
#' @param bam_file Path to BAM file
#' @param gene_config Resolved gene config list (must include $all_exons)
#' @param max_exon_distance Max distance (bp) between a breakpoint and the
#'   nearest exon. Default gene_config$graph_settings$max_exon_distance, else 5000.
#' @param buffer_bp Flanking BAM-scan window in bp. Default
#'   gene_config$graph_settings$buffer_bp, else 15000.
#' @param min_support Minimum distinct reads per cycle. Default
#'   gene_config$graph_settings$min_support, else 2.
#' @param expected_isize_max Minimum mate separation for a discordant pair.
#'   Default gene_config$graph_settings$expected_isize_max, else 1000.
#' @return list of PTD cycles (possibly empty). Each carries the duplicated
#'   exons, the consensus duplicated interval (dup_start, dup_end,
#'   dup_length) and the supporting read names.
extract_graph_candidates <- function(bam_file, gene_config, max_exon_distance = NULL,
                                     buffer_bp = NULL, min_support = NULL,
                                     expected_isize_max = NULL) {
  gs <- gene_config$graph_settings
  if (is.null(gene_config$all_exons) || length(gene_config$all_exons) < 1L) return(list())
  max_exon_distance <- max_exon_distance %||% gs$max_exon_distance %||% 5000L
  buffer_bp <- buffer_bp %||% gs$buffer_bp %||% 15000L
  min_support <- min_support %||% gs$min_support %||% 2L
  expected_isize_max <- expected_isize_max %||% gs$expected_isize_max %||% 1000L

  bam_dt <- .get_evidence_dt(bam_file, gene_config, buffer_bp = buffer_bp)
  if (nrow(bam_dt) == 0L) return(list())
  evidence <- rbind(
    .parse_sa_tags(bam_dt),
    .extract_discordant_and_shadow(bam_dt, expected_isize_max = expected_isize_max)
  )
  .graph_cycles_from_evidence(evidence, gene_config$all_exons,
                              max_exon_distance = max_exon_distance,
                              min_support = min_support)
}

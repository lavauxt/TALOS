=========================================================================

#' Map Raw Genomic Evidence Coordinates to Annotated Exon Graph Edges
#'
#' @param evidence_dt data.frame of split reads, discordant mates, shadows
#'   (from .parse_sa_tags() / .extract_discordant_and_shadow())
#' @param exon_granges Annotated GRanges for exons (must carry an exon_idx column)
#' @param max_distance Maximum genomic distance (bp) to assign a boundary
#'   mapping. NOTE: this was hardcoded to 250bp in the original version, which
#'   silently drops most evidence for deep-intronic PTD breakpoints -- KMT2A
#'   introns run into the multiple-kb range, and two of the three fragments in
#'   the sample this module was built for were intronic. Callers should pass a
#'   per-gene value (see extract_graph_candidates()) rather than relying on a
#'   short-intron-tuned default.
#' @return Aggregated edge data.frame with joint probabilities (one row per
#'   distinct (from,to) exon pair)
.map_evidence_to_exon_edges <- function(evidence_dt, exon_granges, max_distance = 5000L) {
  if (is.null(evidence_dt) || nrow(evidence_dt) == 0L) return(data.frame())

  chr_name <- as.character(GenomicRanges::seqnames(exon_granges)[1])
  from_gr  <- GenomicRanges::GRanges(chr_name, IRanges::IRanges(evidence_dt$from_coord, width = 1L))
  to_gr    <- GenomicRanges::GRanges(chr_name, IRanges::IRanges(evidence_dt$to_coord,   width = 1L))

  from_hits <- GenomicRanges::distanceToNearest(from_gr, exon_granges)
  to_hits   <- GenomicRanges::distanceToNearest(to_gr, exon_granges)

  from_dist <- rep(Inf, length(from_gr));           from_dist[S4Vectors::queryHits(from_hits)] <- S4Vectors::mcols(from_hits)$distance
  to_dist   <- rep(Inf, length(to_gr));             to_dist[S4Vectors::queryHits(to_hits)]     <- S4Vectors::mcols(to_hits)$distance
  from_exon <- rep(NA_integer_, length(from_gr));   from_exon[S4Vectors::queryHits(from_hits)]  <- exon_granges$exon_idx[S4Vectors::subjectHits(from_hits)]
  to_exon   <- rep(NA_integer_, length(to_gr));     to_exon[S4Vectors::queryHits(to_hits)]      <- exon_granges$exon_idx[S4Vectors::subjectHits(to_hits)]

  valid_idx <- which(from_dist <= max_distance & to_dist <= max_distance &
                      !is.na(from_exon) & !is.na(to_exon))
  if (length(valid_idx) == 0L) return(data.frame())

  mapped <- data.frame(
    from_exon = paste0("Exon_", from_exon[valid_idx]),
    to_exon   = paste0("Exon_", to_exon[valid_idx]),
    prob      = evidence_dt$prob[valid_idx],
    qname     = evidence_dt$qname[valid_idx],
    type      = evidence_dt$evidence_type[valid_idx],
    stringsAsFactors = FALSE
  )
  mapped <- mapped[mapped$from_exon != mapped$to_exon, , drop = FALSE]
  if (nrow(mapped) == 0L) return(data.frame())

  edge_key <- paste(mapped$from_exon, mapped$to_exon, sep = "->")
  agg <- lapply(split(seq_len(nrow(mapped)), edge_key), function(rows) {
    sub <- mapped[rows, , drop = FALSE]
    data.frame(
      from = sub$from_exon[1], to = sub$to_exon[1],
      prob = 1 - prod(1 - sub$prob),
      support_count = nrow(sub),
      qnames = paste(unique(sub$qname), collapse = ";"),
      evidence_type = paste(unique(sub$type), collapse = "+"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, agg)
}

#' Minimal Dijkstra shortest path over a small named-node edge list
#'
#' Exon graphs here are tiny (dozens of nodes at most), so a plain O(V^2)
#' Dijkstra in base R is more than fast enough and avoids an igraph
#' dependency for a single-purpose lookup.
#'
#' @param edges data.frame(from, to, weight)
#' @param nodes character vector of all node names
#' @param from,to start/end node names
#' @return character vector of node names on the shortest path (empty if unreachable)
.shortest_path_dijkstra <- function(edges, nodes, from, to) {
  n <- length(nodes)
  idx <- stats::setNames(seq_len(n), nodes)
  adj <- vector("list", n)
  for (i in seq_len(nrow(edges))) {
    a <- idx[[edges$from[i]]]; b <- idx[[edges$to[i]]]
    if (is.null(a) || is.null(b)) next
    adj[[a]] <- rbind(adj[[a]], c(b, edges$weight[i]))
  }
  dist <- rep(Inf, n); prev <- rep(NA_integer_, n); visited <- rep(FALSE, n)
  from_i <- idx[[from]]
  if (is.null(from_i) || is.null(idx[[to]])) return(character(0))
  dist[from_i] <- 0

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

#' Identify PTD Back-Edges and Solve Maximum-Likelihood Duplication Cycles
#'
#' @param edges data.frame(from, to, weight, support_count, qnames, evidence_type)
#' @param nodes character vector of exon node names, in transcript order
#' @return list of detected PTD candidates (one entry per back-edge), or NULL
.detect_ptd_cycles <- function(edges, nodes) {
  if (is.null(edges) || nrow(edges) == 0L) return(NULL)

  from_idx <- as.numeric(gsub("Exon_", "", edges$from))
  to_idx   <- as.numeric(gsub("Exon_", "", edges$to))
  back_mask  <- (from_idx >= to_idx) & (edges$evidence_type != "canonical")
  back_edges <- edges[back_mask, , drop = FALSE]
  if (nrow(back_edges) == 0L) return(NULL)

  edge_pair_key <- paste(edges$from, edges$to, sep = ",")

  results <- list()
  for (i in seq_len(nrow(back_edges))) {
    start_node <- back_edges$to[i]   
    end_node   <- back_edges$from[i]  

    node_sequence <- .shortest_path_dijkstra(edges, nodes, start_node, end_node)
    if (length(node_sequence) == 0L) next

    path_pairs <- if (length(node_sequence) > 1L)
      paste(node_sequence[-length(node_sequence)], node_sequence[-1], sep = ",") else character(0)
    confidence_score <- sum(edges$weight[edge_pair_key %in% path_pairs])

    results[[length(results) + 1L]] <- list(
      ptd_signature  = paste(end_node, "->", start_node),
      dup_start_exon = start_node,
      dup_end_exon   = end_node,
      internal_path  = node_sequence,
      nll_score      = confidence_score,
      support_reads  = back_edges$support_count[i],
      support_qnames = back_edges$qnames[i],
      evidence_type  = back_edges$evidence_type[i]
    )
  }
  results
}

#' Main Graph Candidate Extraction Workflow
#'
#' @param bam_file Path to BAM file
#' @param gene_config Resolved gene config list (must include $all_exons with
#'   an exon_idx column for this to run)
#' @param max_exon_distance Max genomic distance (bp) for mapping evidence to
#'   the nearest exon. Defaults to gene_config$graph_settings$max_exon_distance
#'   if set in the YAML, otherwise 5000L.
#' @param buffer_bp Flanking BAM-scan window in bp (see .get_evidence_dt).
#'   Defaults to gene_config$graph_settings$buffer_bp if set, otherwise 15000L.
#' @return list of candidate PTD cycle structures (possibly empty list())
extract_graph_candidates <- function(bam_file, gene_config,
                                      max_exon_distance = NULL,
                                      buffer_bp = NULL) {
  if (is.null(gene_config$all_exons) || length(gene_config$all_exons) < 2L) return(list())

  if (is.null(max_exon_distance))
    max_exon_distance <- gene_config$graph_settings$max_exon_distance %||% 5000L
  if (is.null(buffer_bp))
    buffer_bp <- gene_config$graph_settings$buffer_bp %||% 15000L

  bam_dt <- .get_evidence_dt(bam_file, gene_config, buffer_bp = buffer_bp)
  if (nrow(bam_dt) == 0L) return(list())

  split_dt <- .parse_sa_tags(bam_dt)
  pair_dt  <- .extract_discordant_and_shadow(bam_dt)
  all_ev   <- rbind(
    if (nrow(split_dt) > 0L) split_dt else NULL,
    if (nrow(pair_dt)  > 0L) pair_dt  else NULL
  )
  if (is.null(all_ev) || nrow(all_ev) == 0L) return(list())

  exons_gr   <- gene_config$all_exons
  ptd_edges  <- .map_evidence_to_exon_edges(all_ev, exons_gr, max_distance = max_exon_distance)

  num_exons  <- length(exons_gr)
  exon_names <- paste0("Exon_", exons_gr$exon_idx)

  canonical_edges <- data.frame(
    from = exon_names[1:(num_exons - 1L)], to = exon_names[2:num_exons],
    prob = 0.99, support_count = 100L, qnames = "", evidence_type = "canonical",
    stringsAsFactors = FALSE
  )

  combined <- rbind(canonical_edges, if (!is.null(ptd_edges) && nrow(ptd_edges) > 0L) ptd_edges else NULL)
  combined$weight <- -log(combined$prob + 1e-9)

  .detect_ptd_cycles(combined, exon_names)
}

# Regression tests for engine_graph_ptd.R helper functions
#
# Covers:
#   - .consensus_coord() precise vs imprecise, all-NA guard
#   - .shortest_path_dijkstra() basic pathfinding
#   - .detect_ptd_cycles() back-edge detection, single-exon duplication
#   - .graph_cycles_from_evidence() NA parameter guard
#   - .map_evidence_to_exon_edges() exon-in-interval mapping

library(testthat)

get_fn <- function(name) {
  if (exists(name, envir = globalenv(), inherits = FALSE)) {
    get(name, envir = globalenv())
  } else if (requireNamespace("TALOS", quietly = TRUE)) {
    getFromNamespace(name, "TALOS")
  } else {
    stop("Function ", name, " not found. Load the package with devtools::load_all().")
  }
}

# ---------------------------------------------------------------------------
# .consensus_coord
# ---------------------------------------------------------------------------
describe(".consensus_coord", {
  cc <- get_fn(".consensus_coord")

  it("returns the mode of precise positions when available", {
    x <- c(100, 101, 100, 103, 100)
    precise <- c(TRUE, TRUE, TRUE, TRUE, TRUE)
    expect_equal(cc(x, precise), 100L)
  })

  it("falls back to rounded median of all positions when no precise", {
    x <- c(100, 101, 102, 103)
    precise <- c(FALSE, FALSE, FALSE, FALSE)
    expect_equal(cc(x, precise), 102L)  # round(median(100:103)) = round(101.5) = 102
  })

  it("returns NA when all values are NA", {
    x <- c(NA_integer_, NA_integer_, NA_integer_)
    precise <- c(TRUE, TRUE, TRUE)
    expect_true(is.na(cc(x, precise)))
  })

  it("returns NA for empty input", {
    expect_true(is.na(cc(integer(0), logical(0))))
  })

  it("handles a single value", {
    expect_equal(cc(42L, TRUE), 42L)
  })

  it("uses only non-NA values for the imprecise fallback", {
    x <- c(NA, 100, 102, NA)
    precise <- c(FALSE, FALSE, FALSE, FALSE)
    expect_equal(cc(x, precise), 101L)  # round(median(100, 102)) = 101
  })
})

# ---------------------------------------------------------------------------
# .shortest_path_dijkstra
# ---------------------------------------------------------------------------
describe(".shortest_path_dijkstra", {
  sp <- get_fn(".shortest_path_dijkstra")

  it("finds a direct path", {
    edges <- data.frame(
      from = c("A", "B"),
      to = c("B", "C"),
      weight = c(1, 1),
      stringsAsFactors = FALSE
    )
    expect_equal(sp(edges, c("A", "B", "C"), "A", "C"), c("A", "B", "C"))
  })

  it("returns empty when no path exists", {
    edges <- data.frame(
      from = "A", to = "B", weight = 1,
      stringsAsFactors = FALSE
    )
    expect_equal(sp(edges, c("A", "B", "C"), "A", "C"), character(0))
  })

  it("returns single node when from == to", {
    edges <- data.frame(
      from = "A", to = "B", weight = 1,
      stringsAsFactors = FALSE
    )
    expect_equal(sp(edges, c("A", "B"), "A", "A"), "A")
  })

  it("finds the shortest of multiple paths", {
    edges <- data.frame(
      from = c("A", "A", "B"),
      to = c("B", "C", "C"),
      weight = c(1, 5, 1),
      stringsAsFactors = FALSE
    )
    # A->B->C has cost 2, A->C has cost 5
    expect_equal(sp(edges, c("A", "B", "C"), "A", "C"), c("A", "B", "C"))
  })

  it("returns empty when from or to is not in nodes", {
    edges <- data.frame(from = "A", to = "B", weight = 1, stringsAsFactors = FALSE)
    expect_equal(sp(edges, c("A", "B"), "X", "B"), character(0))
    expect_equal(sp(edges, c("A", "B"), "A", "X"), character(0))
  })

  it("handles exon-named nodes", {
    edges <- data.frame(
      from = c("Exon_1", "Exon_2", "Exon_3"),
      to = c("Exon_2", "Exon_3", "Exon_4"),
      weight = c(0.01, 0.01, 0.01),
      stringsAsFactors = FALSE
    )
    nodes <- c("Exon_1", "Exon_2", "Exon_3", "Exon_4")
    expect_equal(sp(edges, nodes, "Exon_1", "Exon_4"),
                 c("Exon_1", "Exon_2", "Exon_3", "Exon_4"))
  })
})

# ---------------------------------------------------------------------------
# .map_evidence_to_exon_edges
# ---------------------------------------------------------------------------
describe(".map_evidence_to_exon_edges", {
  map_fn <- get_fn(".map_evidence_to_exon_edges")

  exons <- data.frame(
    start = c(100, 300, 500, 700),
    end = c(199, 399, 599, 799),
    exon_idx = 1:4,
    stringsAsFactors = FALSE
  )

  it("maps a duplication spanning exons 2-3 to a back-edge", {
    ev <- data.frame(
      qname = "r1",
      start = 250L,
      end = 650L,
      prob = 0.9,
      evidence_type = "split_read",
      precise = TRUE,
      stringsAsFactors = FALSE
    )
    result <- map_fn(ev, exons)
    expect_equal(nrow(result), 1L)
    expect_equal(result$from, "Exon_3")  # last duplicated exon
    expect_equal(result$to, "Exon_2")    # first duplicated exon
  })

  it("maps a single-exon duplication as a self-loop", {
    ev <- data.frame(
      qname = "r1",
      start = 310L,
      end = 390L,
      prob = 0.9,
      evidence_type = "split_read",
      precise = TRUE,
      stringsAsFactors = FALSE
    )
    result <- map_fn(ev, exons)
    expect_equal(nrow(result), 1L)
    expect_equal(result$from, "Exon_2")
    expect_equal(result$to, "Exon_2")
  })

  it("returns empty when no exon is inside the interval", {
    ev <- data.frame(
      qname = "r1",
      start = 210L,
      end = 290L,  # between exon 1 and 2, no exon fully inside
      prob = 0.9,
      evidence_type = "split_read",
      precise = TRUE,
      stringsAsFactors = FALSE
    )
    result <- map_fn(ev, exons)
    expect_equal(nrow(result), 0L)
  })

  it("returns empty when evidence is too far from any exon", {
    ev <- data.frame(
      qname = "r1",
      start = 10000L,
      end = 10500L,
      prob = 0.9,
      evidence_type = "split_read",
      precise = TRUE,
      stringsAsFactors = FALSE
    )
    result <- map_fn(ev, exons, max_distance = 1000L)
    expect_equal(nrow(result), 0L)
  })

  it("aggregates multiple reads for the same exon pair", {
    ev <- data.frame(
      qname = c("r1", "r2", "r3"),
      start = c(250L, 255L, 260L),
      end = c(650L, 655L, 660L),
      prob = c(0.8, 0.9, 0.7),
      evidence_type = c("split_read", "split_read", "discordant_pair"),
      precise = c(TRUE, TRUE, FALSE),
      stringsAsFactors = FALSE
    )
    result <- map_fn(ev, exons)
    expect_equal(nrow(result), 1L)
    expect_equal(result$support_count, 3L)
    expect_equal(result$n_split, 2L)
    expect_equal(result$n_pairs, 1L)
  })

  it("returns empty for empty evidence", {
    result <- map_fn(data.frame(), exons)
    expect_equal(nrow(result), 0L)
  })

  it("returns empty for NULL evidence", {
    result <- map_fn(NULL, exons)
    expect_equal(nrow(result), 0L)
  })
})

# ---------------------------------------------------------------------------
# .detect_ptd_cycles
# ---------------------------------------------------------------------------
describe(".detect_ptd_cycles", {
  detect_fn <- get_fn(".detect_ptd_cycles")
  nodes <- c("Exon_1", "Exon_2", "Exon_3", "Exon_4")

  it("detects a back-edge from Exon_3 to Exon_2", {
    edges <- data.frame(
      from = c("Exon_1", "Exon_2", "Exon_3", "Exon_3"),
      to = c("Exon_2", "Exon_3", "Exon_4", "Exon_2"),
      prob = c(0.99, 0.99, 0.99, 0.8),
      weight = c(-log(0.99 + 1e-9), -log(0.99 + 1e-9), -log(0.99 + 1e-9), -log(0.8 + 1e-9)),
      support_count = c(100L, 100L, 100L, 3L),
      qnames = c("", "", "", "r1;r2;r3"),
      evidence_type = c("canonical", "canonical", "canonical", "split_read+discordant_pair"),
      stringsAsFactors = FALSE
    )
    result <- detect_fn(edges, nodes)
    expect_length(result, 1L)
    expect_equal(result[[1]]$dup_start_exon, "Exon_2")
    expect_equal(result[[1]]$dup_end_exon, "Exon_3")
    expect_equal(result[[1]]$internal_path, c("Exon_2", "Exon_3"))
  })

  it("handles single-exon duplication (self-loop back-edge)", {
    edges <- data.frame(
      from = c("Exon_1", "Exon_2", "Exon_2"),
      to = c("Exon_2", "Exon_3", "Exon_2"),
      prob = c(0.99, 0.99, 0.8),
      weight = c(-log(0.99 + 1e-9), -log(0.99 + 1e-9), -log(0.8 + 1e-9)),
      support_count = c(100L, 100L, 5L),
      qnames = c("", "", "r1;r2;r3;r4;r5"),
      evidence_type = c("canonical", "canonical", "split_read"),
      stringsAsFactors = FALSE
    )
    result <- detect_fn(edges, nodes)
    expect_length(result, 1L)
    expect_equal(result[[1]]$dup_start_exon, "Exon_2")
    expect_equal(result[[1]]$dup_end_exon, "Exon_2")
    expect_equal(result[[1]]$internal_path, "Exon_2")
  })

  it("returns NULL when no back-edges exist", {
    edges <- data.frame(
      from = c("Exon_1", "Exon_2"),
      to = c("Exon_2", "Exon_3"),
      prob = c(0.99, 0.99),
      weight = c(0.01, 0.01),
      support_count = c(100L, 100L),
      qnames = c("", ""),
      evidence_type = c("canonical", "canonical"),
      stringsAsFactors = FALSE
    )
    result <- detect_fn(edges, nodes)
    expect_null(result)
  })

  it("sorts results by evidence_score (best first)", {
    edges <- data.frame(
      from = c("Exon_1", "Exon_2", "Exon_3", "Exon_3", "Exon_4", "Exon_4"),
      to = c("Exon_2", "Exon_3", "Exon_4", "Exon_2", "Exon_4", "Exon_1"),
      prob = c(0.99, 0.99, 0.99, 0.9, 0.99, 0.5),
      weight = c(-log(0.99 + 1e-9), -log(0.99 + 1e-9), -log(0.99 + 1e-9),
                 -log(0.9 + 1e-9), -log(0.99 + 1e-9), -log(0.5 + 1e-9)),
      support_count = c(100L, 100L, 100L, 10L, 100L, 2L),
      qnames = c("", "", "", "r1;r2;r3;r4;r5;r6;r7;r8;r9;r10", "", "r1;r2"),
      evidence_type = c("canonical", "canonical", "canonical", "split_read",
                        "canonical", "discordant_pair"),
      stringsAsFactors = FALSE
    )
    result <- detect_fn(edges, nodes)
    expect_gte(length(result), 1L)
    # Verify sorted by evidence_score descending
    scores <- vapply(result, function(r) r$evidence_score, numeric(1L))
    expect_equal(scores, sort(scores, decreasing = TRUE))
  })
})

# ---------------------------------------------------------------------------
# .graph_cycles_from_evidence (NA parameter guard)
# ---------------------------------------------------------------------------
describe(".graph_cycles_from_evidence (NA parameter guard)", {
  gcf <- get_fn(".graph_cycles_from_evidence")

  exons <- data.frame(
    start = c(100, 300),
    end = c(199, 399),
    exon_idx = 1:2,
    stringsAsFactors = FALSE
  )

  ev <- data.frame(
    qname = c("r1", "r2", "r3"),
    start = c(150L, 155L, 160L),
    end = c(350L, 355L, 360L),
    prob = c(0.9, 0.8, 0.85),
    evidence_type = c("split_read", "split_read", "split_read"),
    precise = c(TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )

  it("handles NA max_exon_distance by falling back to default", {
    result <- gcf(ev, exons, max_exon_distance = NA_integer_, min_support = 2L)
    # Should not error; may return cycles or empty list
    expect_true(is.list(result))
  })

  it("returns empty list for NULL evidence", {
    result <- gcf(NULL, exons)
    expect_length(result, 0L)
  })

  it("returns empty list for empty exons", {
    result <- gcf(ev, data.frame(start = integer(0), end = integer(0),
                                  exon_idx = integer(0)))
    expect_length(result, 0L)
  })
})

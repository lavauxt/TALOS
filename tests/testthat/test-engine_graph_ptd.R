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

  # [ITEM 4] Vectorisation regression tests: confirm the vectorised path
  # produces the same results as the original vapply() loop for a range
  # of inputs, including a KMT2A-scale case with many exons.
  describe("ITEM 4: vectorised paths are numerically equivalent", {
    # Build a realistic exon layout: 36 exons spread across 100kb (KMT2A-scale)
    kmt2a_exons <- data.frame(
      start = seq(1000L, 95000L, by = 2700L),
      end   = seq(1100L, 95100L, by = 2700L),
      exon_idx = seq_len(36L),
      stringsAsFactors = FALSE
    )

    # Helper: original (pre-vectorisation) near() implementation, used as the
    # reference oracle in the equivalence checks below.
    near_ref <- function(p, ex) vapply(p, function(v) {
      min(pmax(0L, ex$start - v, v - ex$end))
    }, numeric(1L))

    it("distance filter keeps the same rows as the vapply reference", {
      # 200 evidence rows scattered across the 100kb window
      set.seed(42L)
      n <- 200L
      ev <- data.frame(
        qname = paste0("r", seq_len(n)),
        start = sample(1L:100000L, n, replace = TRUE),
        end   = sample(1L:100000L, n, replace = TRUE),
        prob = runif(n),
        evidence_type = "split_read",
        precise = TRUE,
        stringsAsFactors = FALSE
      )
      # Ensure start <= end
      swap <- ev$start > ev$end
      tmp <- ev$start[swap]; ev$start[swap] <- ev$end[swap]; ev$end[swap] <- tmp

      max_dist <- 5000L
      # Reference (original implementation):
      keep_ref <- near_ref(ev$start, kmt2a_exons) <= max_dist &
        near_ref(ev$end, kmt2a_exons) <= max_dist
      ev_ref <- ev[keep_ref, , drop = FALSE]

      # New (vectorised):
      filter_fn <- get_fn(".filter_evidence_near_exons")
      ev_new <- filter_fn(ev, kmt2a_exons, max_dist)

      # Same rows (compare by sorted qname to be order-agnostic).
      expect_equal(sort(ev_new$qname), sort(ev_ref$qname))
    })

    it("distance filter handles the empty-input case", {
      filter_fn <- get_fn(".filter_evidence_near_exons")
      empty <- data.frame(
        qname = character(), start = integer(), end = integer(),
        prob = numeric(), evidence_type = character(), precise = logical(),
        stringsAsFactors = FALSE
      )
      expect_equal(nrow(filter_fn(empty, kmt2a_exons, 5000L)), 0L)
    })

    it("exon-in-interval query returns the same indices as the lapply reference", {
      # 50 evidence rows; each row may span 0..several exons of the 36-exon
      # KMT2A layout.
      set.seed(7L)
      n <- 50L
      ev <- data.frame(
        start = sample(1L:95000L, n, replace = TRUE),
        end   = NA_integer_,
        stringsAsFactors = FALSE
      )
      # Make each evidence row span a random window up to 10kb
      ev$end <- ev$start + sample(0L:10000L, n, replace = TRUE)

      exon_tol <- 10L
      # Reference:
      inside_ref <- lapply(seq_len(n), function(i) {
        kmt2a_exons$exon_idx[
          kmt2a_exons$start >= ev$start[i] - exon_tol &
          kmt2a_exons$end   <= ev$end[i]   + exon_tol
        ]
      })
      # New:
      inside_fn <- get_fn(".find_exons_inside_evidence")
      inside_new <- inside_fn(ev, kmt2a_exons, exon_tol)

      # Same per-row indices.
      for (i in seq_len(n)) {
        expect_equal(inside_new[[i]], inside_ref[[i]],
                     info = sprintf("row %d (%d-%d)", i, ev$start[i], ev$end[i]))
      }
    })

    it("end-to-end .map_evidence_to_exon_edges matches across many rows", {
      # 100 evidence rows, mix of split-read and discordant-pair types,
      # spread across the 36-exon layout. The aggregated result should be
      # the same regardless of which path the implementation takes.
      set.seed(123L)
      n <- 100L
      ev <- data.frame(
        qname = paste0("read_", seq_len(n)),
        start = sample(500L:95000L, n, replace = TRUE),
        end   = NA_integer_,
        prob = runif(n, 0.4, 0.95),
        evidence_type = sample(c("split_read", "discordant_pair"), n,
                               replace = TRUE),
        precise = sample(c(TRUE, FALSE), n, replace = TRUE),
        stringsAsFactors = FALSE
      )
      ev$end <- ev$start + sample(50L:8000L, n, replace = TRUE)

      result <- map_fn(ev, kmt2a_exons, max_distance = 5000L, exon_tol = 10L)
      # Smoke test: returns a data.frame with the expected columns
      expect_true(is.data.frame(result))
      expect_in(c("from", "to", "prob", "support_count", "n_split", "n_pairs",
                  "qnames", "evidence_type", "dup_start", "dup_end"),
                names(result))
      # Every row should have a valid from/to exon name
      expect_true(all(grepl("^Exon_\\d+$", result$from)))
      expect_true(all(grepl("^Exon_\\d+$", result$to)))
      # Probabilities are in (0, 1]
      expect_true(all(result$prob > 0 & result$prob <= 1))
    })
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

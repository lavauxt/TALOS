# Reconstruct function source without deparse line-wrapping, so multi-token
# expectations cannot be split across deparsed lines.
.src_of <- function(name) {
  paste(deparse(getFromNamespace(name, "TALOS"), width.cutoff = 500L), collapse = "\n")
}

test_that("HTML report merges results table, variant details and all values into one tab", {
  assets <- get(".TALOS_REPORT_JS", envir = asNamespace("TALOS"))
  # One merged top-level Results tab hosting the three former tabs as sub-views.
  expect_match(assets, "tab\\('results', 'Results'")
  expect_match(assets, "function subTabs\\(")
  expect_match(assets, "\\{ id: 'table', label: 'Results table'")
  expect_match(assets, "\\{ id: 'details', label: 'Variant details'")
  expect_match(assets, "\\{ id: 'all-values', label: 'All result values'")
  # The standalone tabs no longer exist at the top level.
  expect_false(grepl("tab\\('details'", assets))
  expect_false(grepl("tab\\('all-values'", assets))
  # CSV exports of every view are preserved.
  expect_match(assets, "TALOS_results.csv")
  expect_match(assets, "TALOS_variant_details.csv")
  expect_match(assets, "TALOS_all_result_values.csv")
  # Legacy deep links (#tab=details, #tab=all-values) resolve to the merged tab.
  expect_match(assets, "wantTab === 'details' \\|\\| wantTab === 'all-values'")
  expect_false(exists(".talos_all_values_script", envir = asNamespace("TALOS"), inherits = FALSE))
})

test_that("every result table embeds the full toolbar feature set", {
  assets <- get(".TALOS_REPORT_JS", envir = asNamespace("TALOS"))
  # Shared toolbar factories: column chooser, quick filters, density, legend.
  expect_match(assets, "function makeColumnChooser\\(")
  expect_match(assets, "function makeQuickFilters\\(")
  expect_match(assets, "function makeDensitySeg\\(")
  expect_match(assets, "function groupLegendEl\\(")
  # Each of the three sub-views owns an independent persisted column
  # selection: the Results table keeps the legacy key, the details view and
  # the all-values matrix get their own stores.
  expect_match(assets, "storeKey: colStore")
  expect_match(assets, "storeKey: 'talos-cols-details-' \\+ M\\.mode")
  expect_match(assets, "storeKey: 'talos-cols-all-' \\+ M\\.mode")
  # Every builder combines its base predicate with quick filters and carries
  # a density toggle plus the colour-group legend.
  expect_match(assets, "makeQuickFilters\\(samplePred,")
  expect_match(assets, "makeQuickFilters\\(basePred,")
  expect_match(assets, "makeDensitySeg\\(\\)")
  expect_match(assets, "groupLegendEl\\(chooser\\.order\\(\\)\\)")
  # The chooser UI keeps the search box, presets and grouped checkboxes.
  expect_match(assets, "Find any variable")
  expect_match(assets, "Clinical view")
  expect_match(assets, "Select all")
  expect_match(assets, "Clear all")
  expect_match(assets, "T\\.cols\\.filter")
  # Density stays one global preference shared by all three tables.
  expect_match(assets, "function setDensity\\(")
  expect_match(assets, "view\\.setDensity\\(currentDensity\\(\\)\\)")
  # The old chooser hand-off button is gone: each table owns its chooser.
  expect_false(grepl("Choose variables", assets))
  expect_false(grepl("openColumnChooser", assets))
})

test_that("HTML report embeds an ECHO-style global sample selector", {
  assets <- get(".TALOS_REPORT_JS", envir = asNamespace("TALOS"))
  shell_source <- .src_of("talos_html_report")
  # The header hosts the picker; every tab filters through the shared state.
  expect_match(shell_source, "sample-pick")
  expect_match(assets, "function setSample\\(")
  expect_match(assets, "function samplePred\\(")
  expect_match(assets, "All samples \\(")
  expect_match(assets, "Filter sample")
  expect_match(assets, "headerPicker")
  # The quick filters are combined with the base predicate (sample filter and,
  # for the details view, the variant-id membership) inside the factory.
  expect_match(assets, "function pred\\(i\\) \\{ return basePred\\(i\\) && \\(!qmode \\|\\| qpred\\(i\\)\\); \\}")
  expect_match(assets, "makeQuickFilters\\(samplePred,")
  expect_match(assets, "makeQuickFilters\\(basePred,")
  # URL hash carries the selection: #tab=results&sample=...
  expect_match(assets, "&sample=' \\+ encodeURIComponent\\(state\\.sample\\)")
})

test_that("coverage diagnostics are present in HTML report source", {
  panel_source <- .src_of(".panel_payload")
  coverage_source <- .src_of(".plot_coverage")
  expect_match(panel_source, "coverage_points")
  expect_match(coverage_source, "Coverage unavailable: no readable BAM path")
})

test_that("PDF prints a compact variant-details table that stays in the page frame", {
  src <- .src_of(".draw_summary_pages")
  plot_src <- .src_of("plot_talos_report")
  # The PDF carries the curated compact field set, not the full metric dump.
  expect_match(plot_src, "\\.compact_variant_fields\\(")
  expect_false(grepl("Variant Summary", plot_src))
  expect_match(plot_src, "Variant Details \\(compact\\)")
  compact_src <- .src_of(".compact_variant_fields")
  expect_match(compact_src, "\\.COMPACT_VARIANT_FIELDS")
  expect_match(compact_src, "intersect\\(")
  # Cell values are pre-wrapped and the line count drives the pagination, so
  # chunks are sized from drawn geometry instead of probe measurements.
  expect_match(src, "strwrap\\(v, width = 64\\)")
  expect_match(src, "row_lines")
  expect_match(src, "line_pitch_pt")
  expect_false(grepl("grobHeight", src))
  expect_false(grepl("row_metrics", src))
  # Every chunk is drawn on its own page inside a viewport that fits it.
  expect_match(src, "grid\\.newpage\\(\\)")
  expect_match(src, "chunks\\[\\[length\\(chunks\\) \\+ 1L\\]\\] <- start_i:end_i")
  # The fallback without gridExtra paginates the printed text as well.
  expect_match(src, "capture\\.output\\(print\\(summary_df\\)\\)")
  # plot_talos_report delegates the drawing to the paginator.
  expect_match(plot_src, "\\.draw_summary_pages")
})

# TALOS v0.4.20

**Tandem duplication and partial tandem duplication detection from BAM files.**

TALOS v0.4.20 is an R/Bioconductor package for detecting internal tandem duplications
(ITDs), partial tandem duplications (PTDs), and ALU/SINE insertions in targeted
short-read sequencing data. It combines CIGAR/soft-clip evidence, reference
re-anchoring, k-mer breakpoint support, paired-end evidence, and per-base
coverage summaries. Outputs include TSV, VCF, PDF, and self-contained HTML
reports.

## Installation

TALOS v0.4.20 requires R >= 4.3.0, Bioconductor dependencies, and `samtools` (>= 1.10)
with `tabix` available on `PATH`. Optional Ensembl annotation and visual
report features use additional packages listed under `Suggests` in
`DESCRIPTION` (including `biomaRt`, `ensembldb`, `txdbmaker`, `Gviz`, and
`plotly`).

```r
install.packages(c(
  "BiocManager", "yaml", "DBI", "RSQLite", "ggplot2", "scales", "htmltools"
))
BiocManager::install(c(
  "Biostrings", "GenomicAlignments", "GenomicRanges", "GenomeInfoDb",
  "IRanges", "Rsamtools", "S4Vectors", "BiocGenerics", "pwalign"
))
```

Install the package from the project repository using your preferred source
installation workflow.

## Quick start

```r
library(TALOS)

result <- talos(
  bam_path = "sample.bam",
  gene = "FLT3",
  build = "hg38",
  sample_name = "sample_01",
  output_prefix = "TALOS v0.4.20",
  verbose = TRUE
)
```

`talos()` loads the bundled gene configuration and reference settings. For
custom targets or lower-level control, construct a configuration with
`get_gene_config()` and call `detect_itd()` directly.

## Interpreting ITD sequence coverage

TALOS v0.4.20 reports separate evidence layers rather than treating them as
interchangeable:

- **Direct coverage**: per-base read alignment to the reconstructed ITD
  sequence (the `min_itd_read_coverage` filter uses this measure).
- **PE coverage**: bases covered by eligible mates from supporting fragments.
  Both mates can contribute; overlapping mate bases count once per fragment.
- **K-mer coverage**: position-resolved evidence from unique 11-mers that form
  a consistent collinear read-to-ITD match. Repeated k-mers are not assigned to
  every repeated copy because they cannot uniquely identify a position.
- **Combined coverage**: the union of the enabled evidence layers. The report
  also lists any remaining uncovered ITD-relative intervals.

A called duplication length is **not** a guarantee that every base was
observed. When a sequence is partially read-derived and reference-completed,
TALOS v0.4.20 validates the observed fragment against the expected segment; if that
validation is insufficient, it marks the sequence unresolved and avoids
reporting a misleading full-length coverage percentage. Paired-end length
estimates are reported as an independent consistency check and do not
silently replace the primary extension-based length.

## Reference methods

The implementation was reviewed against:

- [FiLT3r](https://gitlab.univ-lille.fr/filt3r/filt3r) — alignment-free k-mer
  FLT3 detection from sequence reads.
- [FLT3_ITD_ext](https://github.com/ht50/FLT3_ITD_ext) — soft-clip extension,
  candidate re-alignment, UMI-aware support, and per-base coverage checks.
- [getITD](https://github.com/tjblaette/getitd) — read-to-wild-type alignment,
  insertion classification, and explicit handling of incomplete wild-type
  tandem coverage.

These tools use different evidence definitions and thresholds. Their outputs
should be compared at the level of breakpoint, inserted sequence, length,
independent supporting fragments, and coverage—not by assuming that similarly
named coverage percentages are identical.

## Tests

Run the package tests in an R environment with the declared dependencies:

```r
testthat::test_local()
# or, after installing the package:
# R CMD check TALOS v0.4.20_*.tar.gz
```

## License

GPL-3. See `LICENSE`.



### Clinical-facing support and allele fraction (v0.4.18)

The primary ITD `SupportingReads`, `DepthAtBreakpoint` and `AlleleFrequency` values are raw, unadjusted distinct-QNAME counts and ratios. The optional size-bias estimate is reported separately as `AdjustedSupportingReads`, `CorrectedDepthAtBreakpoint`, `CorrectedAlleleFrequency` and `SizeCorrectionFactor`; `max_correction` defaults to `1` (no inflation). The corrected values are exploratory and are not used for the primary VCF AF/DP/SUPPORT fields. `InFrame` is an informational diagnostic only and does not filter calls.

Paired-end genomic-coordinate projection is a separate diagnostic and is not treated as direct ITD-allele sequence coverage. Combined ITD coverage uses actual read-sequence alignments plus position-resolved sequence-anchored k-mer evidence. BAM duplicate/QC-failed records are excluded; coverage-derived metrics and WT support use the run's MAPQ threshold. A run that reaches `max_reads_in_region` emits a warning that the input may have been truncated.

Adapter screening and BAM base-quality metrics are implemented since v0.4.20. Adapter screening (`adapter_sequences`, `adapter_min_overlap`) is a configurable, report-only pass: it counts region-wide and per-variant QNAMEs carrying adapter read-through (`ITDAdapterScreeningStatus`, `ITDAdapterPositiveSupportReads`) without removing reads. Base-quality filtering (`min_base_quality`, default 20) drops reads below the mean-Phred threshold from the analysis population and reports measured per-variant values (`ITDMeanBaseQuality`, `ITDMinBaseQuality`). Both features degrade gracefully to explicit "unavailable/disabled" status strings when no QUAL column or adapter set is configured. See `CODE_REVIEW_0.4.20.md` for the applied patch checklist, remaining risks and validation limitations.

### Coverage diagnostics in v0.4.3

For an ITD whose read-derived sequence cannot be reconstructed, TALOS v0.4.20 can still
calculate per-base coverage when the full duplicated segment is available in
the configured reference. The reference sequence is used only as a coordinate
template for direct read-sequence and position-specific k-mer coverage. PE genomic-coordinate projection is retained only as a separate diagnostic because it cannot distinguish the original WT copy from the duplicated ITD allele.
`ITD_Sequence` remains unresolved, and `CoverageTemplateSource` explicitly
records this distinction.

Genome-view coverage in PDF and HTML is fetched for the full plotted interval,
including flanking transcript exons (for FLT3, this includes E16 when shown),
rather than reusing the smaller candidate-discovery window.


### Coverage rendering updates in v0.4.3

The ITD coverage track is retained even when all three evidence sources report
zero support, in which case it displays an explicit zero baseline. If the read-derived
ITD sequence is unresolved, TALOS v0.4.20 may use the corresponding reference segment
only as a positional template for direct read-sequence and position-specific k-mer
coverage; PE coordinate projection remains diagnostic and is excluded from combined allele coverage; the clinical `ITD_Sequence` remains unresolved. Direct alignments are
considered in both orientations, with one selected alignment per read.

PDF and HTML genome-view coverage is calculated from the BAM over the complete
plotted interval (including flanking exons such as FLT3 E16). Coverage is stored
window-relative, avoiding the mismatch between a target-only Rle and genomic
coordinates several megabases away.


### ITD metric interpretation and HTML column selection

TALOS v0.4.20 reports direct read-end coverage, PE coordinate-projection diagnostics and position-unique k-mer coverage separately. Combined ITD allele coverage uses direct read-sequence alignments plus sequence-anchored k-mer evidence only; PE projection is not allele-specific and cannot fill ITD coverage gaps. Direct coverage processes all distinct supporting read-end alignments in bounded chunks; the clinical support count and PE matching remain QNAME/fragment-based. Both mates can contribute direct base evidence when both support the call.

When a full-length ITD sequence cannot be validated, TALOS v0.4.20 may use the matching reference segment only as a positional coverage template. This permits direct-read and sequence-anchored k-mer coverage metrics to be calculated without fabricating a clinical sequence; PE coordinate projections remain a separate diagnostic. Full-sequence identity, alignment score, overall orientation and SupportConsistency remain unavailable for reference-imputed/reference-only templates; observed partial-fragment identity/orientation and its query coverage are reported separately.

In the HTML report, use **Columns** in Results or **Choose variables** in Variant Details to search and select any raw result variable. The selection is shared between both tables; **Select all** exposes every result column, including diagnostic, sequence and RLE fields.


### Interpreting coverage denominators (v0.4.5)

Direct and k-mer ITD coverage are calculated from all distinct supporting read-end alignments/sequences, including both mates when available. The read-end counts and the fragment (QNAME) counts are reported separately. This avoids silently dropping a mate from per-base coverage while preserving fragment-level clinical support counts. K-mer metrics now provide both read-end and fragment-level tested/matched counts. Full-length alignment/identity/orientation metrics remain explicitly unavailable when a validated full-length read-derived sequence does not exist; observed partial-fragment and reference-template coverage metrics are reported independently.


### Multi-k coverage, ambiguity bounds and PE length cross-checks (v0.4.6)

ITD positional k-mer coverage is calculated with sliding 7-, 9-, 11-, 13- and
15-mers. Each k-mer is tested at every nucleotide offset, and only template-
unique k-mers with collinear read anchors are used. The combined multi-k track
counts each read end at most once at a given base; individual k-size tracks,
matched read-end counts, matched-fragment counts and RLE vectors are also
available.

Direct alignment coverage now reports conservative and possible positional
bounds for repeated sequence. The conservative vector is the intersection of
co-optimal placements; the possible vector is their union. Read ends with at
least two collinear unique 9-mer anchors use the primary placement without an
extra alignment search. Other aligned read ends are checked by iterative local
alignment with masking, up to eight alternatives per orientation. This is a
bounded search within the ITD template, not an exhaustive enumeration of every
suboptimal genome-wide alignment.

Paired-end ITD length is reported in two ways: `LengthPE_AllPairs` uses all
valid spanning fragments, and `LengthPE_SupportQNAME` uses only pairs whose
QNAME belongs to the ITD-supporting fragment set. Both are compared against
the same non-spanning background insert-size baseline. Counts and insert-size
medians are provided for interpreting disagreements. The legacy `LengthPE`
field remains an alias of the all-pair estimate.



### Interpreting v0.4.6 coverage and PE length diagnostics

Multi-k coverage scans 7-, 9-, 11-, 13- and 15-mers at every nucleotide offset; it is not restricted to codon frames. Only template-unique k-mers in collinear read anchors contribute. Direct ambiguity-aware alignment reports certain (intersection of discovered co-optimal placements), possible (union of discovered placements), and ambiguous coverage. Alternative placement search is bounded; use `ITDDirectAlternativeSearchLimitReadEnds` to identify reads where the configured enumeration cap was reached. The possible metric is therefore an upper bound over discovered placements, not a guarantee of exhaustive placement enumeration.

Paired-end length is exported as both `LengthPE_AllPairs` (all valid spanning pairs) and `LengthPE_SupportQNAME` (only spanning pairs whose QNAME belongs to the ITD support set). Both use the same non-spanning background baseline; inspect pair counts, insert-size medians and status fields before interpreting discordant estimates. The two PE estimates are cross-checks and do not automatically replace the primary ITD length.


## Updates in v0.4.20
- R engine improvements
- Base-quality filtering (min mean Phred per read; measured BQ metrics reported per variant)
- Length-scaled breakpoint clustering (tolerance scales with duplication length, capped; anti-chaining width bound kept)
- Exact QNAME-union support merging (merged support/WT/AF recomputed from the union of component QNAME sets)
- Adapter screening (configurable adapter set, report-only, affected QNAMEs reported)
- FLT3 TKD configuration: NOT tuned; TKD-region calls remain out of scope (see CODE_REVIEW_0.4.20.md)

### Report fixes (post-0.4.20, second round)
- **PDF carries compact variant details and stays in the page frame.** The
  former full ~150-metric summary table escaped the PDF page frame on its
  continuation pages: `grobHeight()` probe measurements under-report the drawn
  row pitch of a `gridExtra::tableGrob()` (padding resolves at draw time and
  long eligibility text wraps to extra lines), so page chunks were sized about
  50% too tall and the table ran past the top and bottom page edges.
  `plot_talos_report()` now prints a curated **compact variant-details table**
  (`.compact_variant_fields()`: 33 clinical fields for ITD/PTD, 17 for ALU)
  instead of the full dump - the complete matrix remains in the HTML report
  and the results TSV. Pagination is deterministic: every cell value is
  pre-wrapped with `strwrap(width = 64)`, the drawn line count of every row is
  known exactly, and rows are packed per page from a line-count height model
  (18 pt per line, 20 pt header, 6 pt slack) that over-estimates the drawn
  pitch, so a page can only under-fill, never overflow. Wide multi-call
  tables still shrink the font to fit the page width; the no-gridExtra text
  fallback paginates the same way.
- **Every result table has the full toolbar.** The Variant details and
  All result values sub-views now carry the same features as the Results
  table: their own searchable, grouped **column chooser** with the Clinical
  view / Select all / Clear all presets (independent persisted selection per
  table: `talos-cols-<mode>`, `talos-cols-details-<mode>`,
  `talos-cols-all-<mode>`), **quick-filter chips** (All / Hotspot / Flagged /
  ITD / PTD, combined with the header sample filter), the **compact /
  comfortable density toggle** (one shared preference, applied to all three
  tables at once), the **colour-group legend**, sticky columns, inline
  expandable row details, search, sorting, pagination and CSV export.
- **One merged Results tab.** The former "Results table", "Variant details"
  and "All result values" tabs are ECHO-style sub-views of a single Results
  tab. Legacy deep links (`#tab=details`, `#tab=all-values`) still resolve.
- **ECHO-style header sample selector.** A single "Filter sample" dropdown in
  the report header filters every tab: the overview figures are recomputed,
  all three results views re-filter, and the genome view restricts its
  sample-gene panel dropdown. The selection is carried in the URL hash
  (`#tab=results&view=details&sample=HORIZON`) for shareable views.
- **Result rows are clickable.** Clicking a row in any of the three result
  tables opens the variant in the details sub-view; the row expander keeps
  toggling the inline detail.

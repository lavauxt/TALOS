# TALOS

<img src="assets/logo.png" align="right" height="139" />

Tandem duplications identified via Alignment‑free Lookup Of Skipped k‑mers

[License: GPL v3](https://www.gnu.org/licenses/gpl-3.0)

TALOS is an R package for detecting internal tandem duplications (ITDs) and
partial tandem duplications (PTDs) from BAM files. It uses a k‑mer breakpoint
detection algorithm and provides HGVS‑compliant cDNA and protein annotations.
Hotspot annotation is performed using a bundled CSV file.

## Features

- **k‑mer based detection** – identifies duplication boundaries without relying
  on split‑read alignment, using both direct adjacency and missing k‑mer block
  methods
- **HGVS annotations** – automatic `c.` and `p.` notation for exonic
  duplications, including frameshift notation for non‑multiples of 3
- **Bundled hotspot CSV** – `hotspots.csv` ships with the package and is used
  by default; no internet connection or database build is required. Pass a
  custom CSV via `hotspot_db_path`.
- **Flexible gene configuration** – transcript-based exon lookup via UCSC
  refGene / EnsDb with disk‑cached TxDb; static offline blocks supported in
  `gene_config.yaml`
- **CIGAR pre‑filtering** – optional fast pre‑filter retaining only
  soft‑clipped and insertion‑bearing reads before k‑mer analysis
- **Size‑bias correction** – corrects supporting read counts for
  length‑dependent detection bias
- **Output formats** – TSV and VCF (v4.2) export for downstream analysis
- **Paired-end support metrics** – optional read-pair diagnostics including
  softclip-restricted PE support, event-size/long-span PE counts, and FR/RF/FF/RR
  orientation summaries
- **Visualisation** – publication‑ready PDF plot and interactive plotly widget
  per detected event, plus a self‑contained HTML report; reports use padded exon
  context, compressed intronic display, genomic coverage, ITD coverage, and ITD
  interval tracks

## Installation

Install directly from GitHub using remotes:

    remotes::install_github("lavauxt/TALOS")

Alternatively, build from source:

    devtools::build()
    devtools::install()

## Requirements

- R ≥ 4.1.0 (the native pipe `|>` is used)
- Bioconductor packages (required):
  `Rsamtools`, `Biostrings`, `GenomicRanges`, `GenomicAlignments`, `IRanges`,
  `S4Vectors`, `BiocGenerics`, `GenomeInfoDb`, `GenomicFeatures`,
  `AnnotationDbi`, `BSgenome`
- CRAN packages (required):
  `yaml`
- Bioconductor packages (required for TxDb building):
  `txdbmaker` (Bioconductor ≥ 3.16) or `GenomicFeatures` (older Bioconductor)
- Optional – RefSeq NM_ transcript fallback: `biomaRt`
- Optional – Ensembl ENST transcript lookup: `ensembldb`, `AnnotationFilter`
- Optional – BSgenome sequence data:
  `BSgenome.Hsapiens.UCSC.hg19` or `BSgenome.Hsapiens.UCSC.hg38`
- Optional – offline EnsDb exon lookup:
  `EnsDb.Hsapiens.v75` (hg19) or `EnsDb.Hsapiens.v86` (hg38)
- Optional – accelerated k‑mer matching: `fastmatch`
- Optional – accelerated CIGAR parsing: `cigarillo`
- Optional – alignment-based orientation / consistency scoring: `pwalign`
- Optional – publication PDF plots: `Gviz`
- Optional – interactive plots: `plotly`, `htmlwidgets`
- Optional – PDF report tables: `gridExtra`
- Optional – HTML report: `rmarkdown`, `DT`
- Optional – hotspot SQLite database: `DBI`, `RSQLite`

Install Bioconductor dependencies:

    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c(
      "Rsamtools", "Biostrings", "GenomicRanges", "GenomicAlignments",
      "IRanges", "S4Vectors", "BiocGenerics", "GenomeInfoDb",
      "GenomicFeatures", "ensembldb", "AnnotationFilter",
      "AnnotationDbi", "BSgenome", "txdbmaker"
    ))

## Basic usage

    library(TALOS)

    results <- talos(
        bam_path      = "path/to/sample.bam",
        gene          = "FLT3",
        build         = "hg19",
        min_support   = 10,
        min_size      = 30,
        output_prefix = "FLT3_ITD",
        output_folder = "./results",
        plot          = TRUE
    )

    print(results)

For KMT2A partial tandem duplication (PTD) detection:

    results_kmt2a <- talos(
        bam_path      = "path/to/sample.bam",
        gene          = "KMT2A",
        build         = "hg19",
        output_prefix = "KMT2A_PTD",
        output_folder = "./results",
        plot          = TRUE
    )
    # Gene-specific settings from gene_config.yaml are applied automatically:
    # min_support=3, max_missing_kmers=0.8, max_itd_length=30000,
    # convert_long_to_ptd=TRUE, merge_ptd_intervals=TRUE.
    # Threshold relaxations (min_support, softclip counts) are applied by talos()
    # for KMT2A regardless of ptd_mode.

On the first call for a given genome build, TALOS downloads the UCSC refGene
annotation (~30 s) and caches it to disk so subsequent calls are instant.
A BSgenome data package must be installed for sequence extraction.

## Hotspot annotation

TALOS annotates detected ITDs against a set of known hotspot regions. By
default the package uses the bundled static CSV file (`hotspots.csv`) that is
installed alongside TALOS. This file contains pre‑computed hotspot intervals
for the supported genes and builds (hg19/hg38). **No internet connection is
required.**

If you wish to use a custom hotspot set, provide your own CSV file via the
`hotspot_db_path` parameter.

### CSV format for custom hotspots

Your custom CSV must contain the following exact column names
(case‑sensitive):

| Column | Type      | Description |
|--------|-----------|-------------|
| Gene   | character | Gene symbol (e.g., `"FLT3"`, `"KMT2A"`). Must match the gene names used in TALOS. |
| Build  | character | Genome build: `"hg19"` or `"hg38"`. |
| Start  | integer   | Genomic start coordinate of the hotspot region (1‑based, inclusive). |
| End    | integer   | Genomic end coordinate of the hotspot region (inclusive). |
| Name   | character | Identifier for the hotspot (e.g., `"FLT3_TKD_hotspot_1"`). |

All columns are required. Extra columns are ignored.

Example snippet (`my_hotspots.csv`):

    Gene,Build,Start,End,Name
    FLT3,hg19,123456,123470,FLT3_hotspot_A
    FLT3,hg38,124567,124581,FLT3_hotspot_A
    KMT2A,hg19,987654,987670,KMT2A_hotspot_1

The hotspot matching logic: an ITD is considered to fall into a hotspot if any
part of the duplicated segment overlaps the interval `[Start, End]`:

    ITD_start <= End  AND  ITD_end >= Start

where `ITD_end = GenomicPosition + Length - 1`.

If multiple hotspots match, all names are concatenated with semicolons (`;`)
in the `HotspotName` column.

### Using a custom CSV

    results <- talos(..., hotspot_db_path = "my_hotspots.csv")

## Configuration modes

TALOS uses a YAML configuration file (`inst/extdata/gene_config.yaml`) that
can contain both static coordinates (under `hg19`/`hg38` blocks) and
transcript identifiers for each gene.

| Mode      | Behavior |
|-----------|----------|
| `full`    | Uses static genomic coordinates, exon boundaries, and cDNA sequence from the config file. Requires `hg19:` / `hg38:` block with `exons:` list in YAML. |
| `minimal` | Ignores static data; fetches exon coordinates from EnsDb database using `transcript` and `targeted_exons`. |
| `auto` (default) | Attempts static coordinates first; automatically falls back to transcript‑based lookup. |

The current bundled `gene_config.yaml` does **not** include static coordinate
blocks; all genes therefore use transcript‑based lookup (`use_db = TRUE`)
and require a BSgenome package to be installed. To enable fully offline
operation, use `save_offline_config()` to snapshot a resolved configuration:

    save_offline_config("FLT3", build = "hg19", output_file = "flt3_offline.yaml")

## Supported genes

The bundled `gene_config.yaml` includes transcript identifiers for:

| Gene   | Transcript     | Targeted exons | Notes |
|--------|----------------|----------------|-------|
| FLT3   | NM_004119      | 14–15          | JM domain ITD |
| KMT2A  | NM_005933      | 2–12           | Partial tandem duplication |
| BCOR   | NM_001123385   | 15             | Internal tandem duplication |
| UBTF   | NM_014233      | 13             | Tandem duplication |

Additional genes can be supported by extending the YAML file – provide static
coordinate blocks, a transcript ID, or both.

## Algorithm overview

1. **BAM reading** – reads are fetched from the genomic window defined in the
   gene config.
2. **CIGAR pre‑filter** (optional) – retains only unmapped, soft‑clipped, or
   net‑insertion reads to reduce the search space.
3. **k‑mer indexing** – a hash‑based index is built from the reference
   sequence.
4. **Breakpoint detection per read** – two complementary strategies:
   - *Direct adjacency*: consecutive k‑mers whose positions overlap in the
     reference indicate a duplicated region.
   - *Missing k‑mer block*: a contiguous run of unmatched k‑mers flanked by
     backward‑jumping positions signals a duplication.
5. **Breakpoint clustering** – nearby breakpoints (within
   `cluster_tolerance` bp) are merged; the consensus is taken as the median.
6. **Size‑bias correction** – raw supporting read count is scaled to
   compensate for the lower detection probability of short duplications.
7. **HGVS annotation** – the genomic breakpoint is mapped to cDNA coordinates
   via cumulative exon lengths.
8. **Hotspot annotation** – duplicated segments are intersected with the
   bundled (or custom) hotspot CSV.
9. **Output** – TSV and optional VCF are written; PDF and interactive HTML
   plots are generated per detected ITD.

## Functions

| Function                     | Description |
|------------------------------|-------------|
| `talos()`                    | Simplified entry point – calls `get_gene_config()` then `detect_itd()` |
| `detect_itd()`               | Core detection, clustering, correction, and output |
| `get_gene_config()`          | Loads gene reference from YAML; implements auto/full/minimal logic |
| `annotate_hotspots()`        | Annotates ITDs using a DBI connection or CSV file |
| `compute_hgvs_annotations()` | Maps genomic breakpoint and duplicated sequence to HGVS `c.` and `p.` notation |
| `plot_talos_report()`        | Generates a publication‑ready PDF with coverage track, gene model, and ITD arcs |
| `plot_talos_interactive()`   | Generates a plotly widget for interactive exploration |
| `talos_html_report()`        | Renders a self‑contained HTML report for one or more samples |
| `talos_batch()`              | Runs `talos()` over multiple BAM files and/or genes |
| `talos_batch_dir()`          | Discovers BAM files in a directory and runs `talos_batch()` |
| `save_offline_config()`      | Snapshots a resolved gene config to YAML for reproducible offline use |
| `detect_orientation()`       | Detects whether a duplication is tandem or inverted |
| `build_vcf_header()`         | Constructs VCF v4.2 header lines |
| `build_vcf_records()`        | Formats ITD results as VCF data lines |
| `write_itd_vcf()`            | Writes a complete VCF file from TALOS results |

## Parameters

### `talos()` – simplified interface

| Parameter                  | Default      | Description |
|----------------------------|--------------|-------------|
| `bam_path`                 | –            | Path to an indexed BAM file |
| `gene`                     | –            | Gene name: `"FLT3"`, `"BCOR"`, `"KMT2A"`, or `"UBTF"` |
| `build`                    | `"hg19"`     | Genome build: `"hg19"` or `"hg38"` |
| `min_support`              | `10`         | Minimum bias‑corrected supporting reads required to report a call |
| `min_size`                 | `15`         | Minimum duplication length (bp) |
| `max_correction`           | `2.0`        | Cap on the size‑bias correction factor |
| `plot`                     | `TRUE`       | Generate a PDF plot for each detected ITD (requires `Gviz`) |
| `output_prefix`            | `"TALOS"`    | Base name for output files; set to `NULL` to suppress file writing |
| `output_folder`            | `"./results"` | Directory for TSV, VCF, and PDF output |
| `sample_name`              | `NULL`       | Sample identifier; derived from the BAM filename if `NULL` |
| `hotspot_db_path`          | `NULL`       | Optional path to a custom CSV file (bypasses the bundled CSV) |
| `yaml_path`                | (bundled)    | Path to the gene config YAML; defaults to the installed `gene_config.yaml` |
| `padding`                  | `500`        | Base‑pair padding around targeted exons |
| `bsgenome`                 | `NULL`       | Optional BSgenome object or package‑name string to override the default |
| `html_report`              | `TRUE`       | Write a self‑contained HTML report |
| `vaf_threshold`            | `0.01`       | Minimum allele frequency to report (default 1%) |
| `filter_intronic`          | `FALSE`      | Drop breakpoints not overlapping a targeted exon |
| `annotate_hotspots`        | `TRUE`       | Run hotspot annotation on results |
| `detect_orientation`       | `TRUE`       | Detect tandem vs inverted orientation |
| `compute_alignment_score`  | `TRUE`       | Compute alignment score between ITD and reference |
| `compute_support_bases`    | `TRUE`       | Count total soft-clip bases from supporting reads |
| `compute_consistency`      | `TRUE`       | Estimate fraction of supporting reads containing the ITD |
| `compute_itd_coverage`     | `TRUE`       | Estimate per-position coverage across the ITD sequence |
| `compute_coverage_drop`    | `TRUE`       | Compute fold-change in coverage at the breakpoint |
| `compute_microhomology`    | `TRUE`       | Compute median microhomology length |
| `compute_repeat_entropy`   | `TRUE`       | Compute dinucleotide Shannon entropy around the breakpoint |
| `compute_discordant_ratio` | `TRUE`       | Compute discordant read-pair ratio (requires a second BAM pass) |
| `compute_hgvs`             | `TRUE`       | Generate HGVS `c.` and `p.` notation |
| `add_config_to_report`     | `FALSE`      | Include gene configuration details in PDF and HTML reports |
| `max_pairwise_alignments`  | `30`         | Maximum reads subsampled for pairwise consistency alignment |

### Additional `detect_itd()` parameters

| Parameter                  | Default   | Description |
|----------------------------|-----------|-------------|
| `k`                        | `11`      | k‑mer length |
| `max_missing_kmers`        | `0.5`     | Maximum proportion of unmatched k‑mers before a read is skipped |
| `cluster_tolerance`        | `10`      | Maximum distance (bp) between breakpoints to be merged into one cluster |
| `prefilter`                | `TRUE`    | Enable CIGAR‑based pre‑filtering |
| `min_ins_filter`           | `3`       | Minimum net insertion length for a non‑clipped read to pass pre‑filtering |
| `min_mapq`                 | `20`      | Minimum MAPQ for reads used in wildtype depth calculation |
| `write_vcf`                | `TRUE`    | Write a VCF file in addition to TSV |
| `verbose`                  | `TRUE`    | Print progress messages |
| `debug`                    | `FALSE`   | Print detailed exon‑mapping diagnostics |
| `nominal_read_len`         | `150`     | Expected read length, used in size‑bias correction |
| `annotate_hotspots`        | `TRUE`    | Run hotspot annotation on results |
| `detect_orientation`       | `TRUE`    | Detect whether the duplication is tandem or inverted |
| `filter_intronic`          | `FALSE`   | Silently drop breakpoints not overlapping a targeted exon |
| `ptd_mode`                 | `FALSE`   | Enable soft-clip–only parsing for partial tandem duplications |
| `use_cigar_bp`             | `TRUE`    | Refine PTD breakpoints from CIGAR rather than k‑mer trace |
| `refine_bp`                | `FALSE`   | Enable secondary k‑mer refinement loop for PTD breakpoints |
| `compute_alignment_score`  | `TRUE`    | Compute per-base alignment score between ITD and reference |
| `compute_support_bases`    | `TRUE`    | Count total soft-clip bases from supporting reads |
| `compute_consistency`      | `TRUE`    | Estimate fraction of supporting reads that contain the ITD sequence |
| `compute_itd_coverage`     | `TRUE`    | Estimate per-position coverage across the ITD sequence |
| `compute_coverage_drop`    | `TRUE`    | Compute fold-change in coverage across the breakpoint |
| `compute_microhomology`    | `TRUE`    | Compute median microhomology length at the breakpoint |
| `compute_repeat_entropy`   | `TRUE`    | Compute dinucleotide Shannon entropy around the breakpoint |
| `compute_discordant_ratio` | `TRUE`    | Compute discordant read-pair ratio spanning the breakpoint |
| `compute_hgvs`             | `TRUE`    | Generate HGVS `c.` and `p.` notation (requires transcript) |
| `max_pairwise_alignments`  | `30`      | Maximum reads subsampled for pairwise consistency alignment |
| `max_reads_in_region`      | `200000`  | Hard cap on reads loaded from the BAM; safety guard for large regions |

## Output

`talos()` and `detect_itd()` return a data frame with one row per detected
ITD:

| Column               | Description |
|----------------------|-------------|
| `Sample`             | Sample identifier |
| `Gene`               | Gene symbol |
| `Genome`             | Genome build used |
| `GenomicPosition`    | Breakpoint coordinate (1‑based, start of duplication) |
| `ITD_Sequence`       | Duplicated nucleotide sequence |
| `Length`             | Duplication length (bp) |
| `SupportingReads`    | Bias‑corrected mutant fragment count |
| `WildtypeReads`      | Wild‑type fragment count spanning the breakpoint |
| `DepthAtBreakpoint`  | Total depth (`SupportingReads + WildtypeReads`) |
| `AlleleFrequency`    | Mutant allele frequency (rounded to 4 decimal places) |
| `HGVS_cDNA`          | HGVS cDNA notation, e.g. `c.1705_1837dup` |
| `HGVS_Protein`       | HGVS protein notation, e.g. `p.Tyr569_Gly613dup` |
| `StrandBias`         | Fraction of supporting reads on the reverse strand |
| `MeanSupportMAPQ`    | Mean MAPQ of supporting reads |
| `BreakpointSpread`   | Range of breakpoint positions within cluster (bp) |
| `SoftclipFraction`   | Fraction of supporting reads detected via soft‑clip |
| `UniqueBreakpoints`  | Number of unique breakpoint positions in cluster |
| `CoverageDrop`       | Coverage fold‑change ratio across the breakpoint |
| `MedianMicrohomology`| Median microhomology length (bp) |
| `DiscordantRatio`    | Discordant pair ratio spanning the breakpoint |
| `RepeatEntropy`      | Dinucleotide Shannon entropy around the breakpoint |
| `SequenceImputed`    | `TRUE` if the ITD sequence was imputed from reference |
| `SequencePartial`    | `TRUE` if the ITD sequence is partially read + partially imputed |
| `SequenceSource`     | One of `"observed"`, `"partial_read+ref"`, or `"ref_imputed"` |
| `RefMatch_Observed`  | % identity between the observed (read-derived) portion of the ITD and the reference |
| `RefMatch_Total`     | % identity between the full ITD sequence (including imputed bases) and the reference |
| `ITDReadCoverage`    | % of ITD sequence positions covered by at least one supporting read |
| `ITDCoverageRLE`     | Run-length encoding of per-position read depth across the ITD sequence |
| `SupportConsistency` | % of supporting reads containing the ITD sequence (≥ 90 % local alignment identity) |
| `AlignmentScore`     | Per-base alignment score between the full ITD sequence and reference (0–1) |
| `TotalSupportBases`  | Total soft‑clip bases across all supporting reads |
| `Orientation`        | `"tandem"` or `"inverted"` |
| `Hotspot`            | `TRUE` if the duplication overlaps a known hotspot region |
| `HotspotName`        | Name(s) of the matched hotspot(s) |
| `Region`             | `"exonic"` or `"intronic"` |
| `ExonNumber`         | Exon number of the breakpoint (if exonic) |
| `TranscriptRef`      | Transcript accession used for HGVS annotation |

When `output_prefix` is provided, the following files are written to
`output_folder/<sample_name>/`:

- `<prefix>_<gene>_<timestamp>.tsv` – tab‑separated results table
- `<prefix>_<gene>_<timestamp>.vcf` – VCF v4.2 with SVLEN, GENE, CDS, AA,
  AF, DP, and SUPPORT INFO fields
- `<sample>_<gene>.pdf` – one PDF per sample/gene (when `plot = TRUE`)
- `<sample>_<gene>_interactive.html` – plotly widget (if `plotly` installed)
- `TALOS_report_<gene>_<timestamp>.html` – summary HTML report (if
  `html_report = TRUE`)

## Notes

- An allele frequency below 1% (`AF < 0.01`, the default `vaf_threshold`)
  is silently filtered before output. Lower this threshold to recover weak
  signals.
- Intronic breakpoints produce `NA` in the HGVS columns. Filter with
  `!is.na(HGVS_cDNA)` to retain only exonic calls.
- For duplications whose length is not a multiple of 3, the protein notation
  uses frameshift format: `p.AaaXXXnnnfsTer?`.
- The bundled `hotspots.csv` is used by default. To use a custom CSV, pass its path to `hotspot_db_path`.
- The UCSC refGene TxDb is built on the first `use_db = TRUE` call and
  cached to disk via `tools::R_user_dir("TALOS", "cache")`. Subsequent calls
  load the cached SQLite file instantly. To force a rebuild, delete the cache
  file or call `get_gene_config()` after removing it.
- If internet is unavailable, provide static coordinate blocks in
  `gene_config.yaml` and run with `use_db = FALSE`, or pre-build the cache
  with `save_offline_config()`.
- **`detect_itd()` vs `talos()` defaults** – `min_size` defaults to `15` in
  both `detect_itd()` and `talos()`. Other parameters may differ; consult
  the parameter tables above when calling `detect_itd()` directly.
- **Per-gene YAML overrides** – the gene config YAML supports an optional
  `gene_settings:` block per gene entry. Any key listed there is used as a
  per-gene default that sits between the package default and an explicit
  argument passed by the user. This is useful for applying gene-specific
  sensitivity thresholds without modifying call sites.
- **`save_offline_config` key naming** – `save_offline_config("FLT3", ...)`
  writes the resolved configuration under the key `FLT3_OFFLINE` in the output
  YAML, not `FLT3`. To use the file with `talos()`, rename the top-level key
  back to the plain gene symbol (`FLT3`) before passing the path via
  `yaml_path`.
- **`compute_discordant_ratio = TRUE`** triggers a second BAM read pass to
  load read pairs. Set it to `FALSE` when processing many samples in parallel
  to halve I/O and memory usage if discordant-pair evidence is not required.
- **KMT2A PTD detection strategy** – KMT2A uses `ptd_mode = FALSE` in the
  bundled YAML (k-mer backward-jump path) combined with
  `convert_long_to_ptd = TRUE` and `max_itd_length = 30000`. This means
  duplications up to 30 kb receive a size estimate from k-mer geometry, while
  larger ones are promoted to PTD records (Length = 0). The `min_length = 1000`
  filter applies only to sized calls (Length > 0) and does not discard PTD-
  converted candidates. `talos()` automatically relaxes support and softclip
  thresholds for KMT2A regardless of `ptd_mode`.

## License

GPL (≥ 3)
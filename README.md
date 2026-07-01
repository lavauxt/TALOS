# TALOS

<img src="assets/logo.png" align="right" height="139" />

Tandem duplications identified via Alignment‑free Lookup Of Skipped k‑mers

[License: GPL v3](https://www.gnu.org/licenses/gpl-3.0)

TALOS is an R package for detecting internal tandem duplications (ITDs),
partial tandem duplications (PTDs), and — as of this release — **ALU/SINE
mobile‑element insertions** from BAM files. The core engine uses a k‑mer
breakpoint detection algorithm; the ALU engine uses local Smith‑Waterman
alignment of soft‑clipped sequences against bundled ALU consensus sequences.
All detection modes produce HGVS‑compliant annotations, hotspot overlap, and
self‑contained HTML reports.

## Features

- **k‑mer based ITD/PTD detection** – identifies duplication boundaries
  without relying on split‑read alignment, using both direct adjacency and
  missing k‑mer block methods
- **ALU/SINE insertion detection** – soft‑clip sequences are aligned to
  bundled ALU consensus sequences (AluSx, AluY, AluJb); each insertion is
  characterised by subtype, orientation, target‑site duplication (TSD), poly‑A
  tail length, and 5′ truncation extent
- **HGVS annotations** – automatic `c.` and `p.` notation for exonic events,
  including frameshift notation for non‑multiples of 3
- **Bundled hotspot CSV** – `hotspots.csv` ships with the package; no internet
  connection or database build is required. Pass a custom CSV via
  `hotspot_db_path`.
- **Flexible gene configuration** – transcript‑based exon lookup via UCSC
  refGene / EnsDb with disk‑cached TxDb; static offline blocks supported in
  `gene_config.yaml`; per‑gene `gene_settings:` / `alu_settings:` blocks for
  parameter overrides
- **CIGAR pre‑filtering** – optional fast pre‑filter retaining only
  soft‑clipped and insertion‑bearing reads before analysis
- **Size‑bias correction** – corrects supporting read counts for
  length‑dependent detection bias (ITD/PTD mode)
- **Output formats** – TSV and VCF (v4.2) export for downstream analysis
- **Paired‑end support metrics** – optional read‑pair diagnostics including
  soft‑clip‑restricted PE support, event‑size/long‑span PE counts, and
  FR/RF/FF/RR orientation summaries
- **Standalone HTML report template** – `inst/rmd/TALOS_report.Rmd` is a
  self‑contained R Markdown file rendered by `talos_html_report()`; it
  supports both ITD and ALU modes via a `mode` parameter and can be
  customised independently of the package source
- **Publication‑ready visualisation** – PDF (Gviz) and interactive plotly
  widget per detected event; reports include compressed intronic display,
  genomic coverage, ITD coverage, and ITD interval tracks

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
- CRAN packages (required): `yaml`
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
- Optional – alignment‑based orientation / consistency scoring: `pwalign`
- Optional – ALU alignment: `Biostrings` (also required for ITD sequence ops)
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

### ITD / PTD detection

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

### ALU insertion detection

Use `talos_alu()` to screen a BAM file for ALU/SINE mobile‑element
insertions. Like `talos()`, it resolves the gene window from the YAML config
and applies any `alu_settings:` overrides defined there.

    results_alu <- talos_alu(
        bam_path      = "path/to/sample.bam",
        gene          = "KMT2A",        # or "KMT2A_ALU" for the full gene window
        build         = "hg19",
        output_prefix = "KMT2A_ALU",
        output_folder = "./results"
    )

For lower‑level control, call `detect_alu()` directly with a pre‑resolved
gene config:

    cfg <- get_gene_config("KMT2A", build = "hg19")
    results_alu <- detect_alu(
        bam_path          = "path/to/sample.bam",
        gene_config       = cfg,
        min_support       = 3,
        min_alu_score     = 0.60,
        min_clip_len      = 25,
        cluster_tolerance = 20,
        min_tsd           = 3,
        vaf_threshold     = 0.005
    )

Supply a custom ALU consensus FASTA to add or replace subfamily sequences:

    results_alu <- talos_alu(
        bam_path     = "path/to/sample.bam",
        gene         = "FLT3",
        consensus_fa = "/path/to/my_alu_consensus.fa"
    )

#### ALU consensus FASTA format

The bundled `inst/extdata/alu_consensus.fa` ships with entries for AluSx,
AluY, and AluJb. A custom FASTA must use one sequence per entry named after
the subfamily:

    >AluSx
    GGCCGGGCGCGGTGGCTCACGCCTGTAATCCC...
    >AluY
    GGCCGGGCGCGGTGGCTCACGCCTGTAATCCC...

### Standalone HTML report

`talos_html_report()` renders `inst/rmd/TALOS_report.Rmd` — an external
template that can be customised without modifying package source. It accepts
results from both ITD and ALU modes via the `mode` parameter:

    talos_html_report(
        result_df    = results,
        bam_paths    = c(sample1 = "s1.bam"),
        gene_configs = list(FLT3 = get_gene_config("FLT3")),
        output_file  = "FLT3_report.html",
        title        = "FLT3 ITD Report — Cohort A",
        mode         = "itd"         # or "alu"
    )

To customise the report template, copy it from the installed package:

    file.copy(
        system.file("rmd", "TALOS_report.Rmd", package = "TALOS"),
        "my_report_template.Rmd"
    )

Then edit `my_report_template.Rmd` and render it directly with
`rmarkdown::render()`, passing the same `params` list expected by the
template (`result_df`, `plot_widgets`, `gene_configs`, `has_plotly`,
`show_config`, `title`, `mode`).

On the first call for a given genome build, TALOS downloads the UCSC refGene
annotation (~30 s) and caches it to disk so subsequent calls are instant.
A BSgenome data package must be installed for sequence extraction.

## Hotspot annotation

TALOS annotates detected events against a set of known hotspot regions. By
default the package uses the bundled static CSV file (`hotspots.csv`).
**No internet connection is required.**

If you wish to use a custom hotspot set, provide your own CSV file via
`hotspot_db_path`. Hotspot matching works identically for both ITD and ALU
results (based on genomic overlap).

### CSV format for custom hotspots

| Column | Type      | Description |
|--------|-----------|-------------|
| Gene   | character | Gene symbol, e.g. `"FLT3"`. Must match TALOS gene names. |
| Build  | character | Genome build: `"hg19"` or `"hg38"`. |
| Start  | integer   | Hotspot region start coordinate (1‑based, inclusive). |
| End    | integer   | Hotspot region end coordinate (inclusive). |
| Name   | character | Identifier, e.g. `"FLT3_TKD_hotspot_1"`. |

Example:

    Gene,Build,Start,End,Name
    FLT3,hg19,123456,123470,FLT3_hotspot_A
    KMT2A,hg19,987654,987670,KMT2A_hotspot_1

## Configuration modes

TALOS uses `inst/extdata/gene_config.yaml` which supports both static
coordinate blocks and transcript‑based lookup.

| Mode      | Behavior |
|-----------|----------|
| `full`    | Uses static genomic coordinates and exon boundaries from the YAML. Requires `hg19:` / `hg38:` block with `exons:` list. |
| `minimal` | Fetches exon coordinates from EnsDb using `transcript` and `targeted_exons`. |
| `auto` (default) | Attempts static coordinates first; falls back to transcript‑based lookup. |

To enable fully offline operation, snapshot a resolved config:

    save_offline_config("FLT3", build = "hg19", output_file = "flt3_offline.yaml")

### Per-gene parameter overrides

Each gene entry in the YAML supports two optional override blocks:

**`gene_settings:`** — overrides for `talos()` / `detect_itd()` parameters.  
**`alu_settings:`** — overrides for `talos_alu()` / `detect_alu()` parameters.

Keys in these blocks act as per‑gene defaults: they override the package
default but are themselves overridden by an explicit argument at the call
site. This means calling `talos_alu(gene = "KMT2A", min_support = 5)` will
use `min_support = 5` regardless of what `alu_settings.min_support` says in
the YAML.

Example YAML fragment:

    KMT2A:
      transcript: NM_005933
      targeted_exons: [2,3,4,5,6,7,8,9,10,11,12]
      alu_settings:
        min_support:       3
        min_alu_score:     0.58
        cluster_tolerance: 20
        min_tsd:           3
        vaf_threshold:     0.005

## Supported genes

The bundled `gene_config.yaml` includes entries for:

| Gene        | Transcript     | Targeted exons | Mode           | Notes |
|-------------|----------------|----------------|----------------|-------|
| FLT3        | NM_004119      | 14–15          | ITD            | JM domain ITD |
| KMT2A       | NM_005933      | 2–12           | PTD            | Partial tandem duplication |
| BCOR        | NM_001123385   | 15             | ITD            | Internal tandem duplication |
| UBTF        | NM_014233      | 13             | ITD            | Tandem duplication |
| FLT3_ALU    | NM_004119      | all (1–24)     | ALU            | Whole‑gene ALU screen |
| KMT2A_ALU   | NM_005933      | all (1–30)     | ALU            | Whole‑gene ALU screen |

ALU entries use wider genomic windows (all exons + intronic sequence) because
ALU insertions can fall anywhere in the gene body. Additional genes can be
added by extending the YAML with a `transcript:` and `targeted_exons:` list.

## Algorithm overview

### ITD / PTD

1. **BAM reading** – reads are fetched from the genomic window defined in the
   gene config.
2. **CIGAR pre‑filter** (optional) – retains only unmapped, soft‑clipped, or
   net‑insertion reads to reduce the search space.
3. **k‑mer indexing** – a hash‑based index is built from the reference sequence.
4. **Breakpoint detection per read** – two complementary strategies:
   - *Direct adjacency*: consecutive k‑mers whose positions overlap in the
     reference indicate a duplicated region.
   - *Missing k‑mer block*: a contiguous run of unmatched k‑mers flanked by
     backward‑jumping positions signals a duplication.
5. **Breakpoint clustering** – nearby breakpoints (within `cluster_tolerance`
   bp) are merged; consensus taken as median.
6. **Size‑bias correction** – raw supporting read count is scaled to
   compensate for the lower detection probability of short duplications.
7. **HGVS annotation** – genomic breakpoint is mapped to cDNA coordinates via
   cumulative exon lengths.
8. **Hotspot annotation** – duplicated segments are intersected with the
   bundled (or custom) hotspot CSV.
9. **Output** – TSV, VCF, PDF, and interactive HTML.

### ALU insertion detection

1. **BAM reading** – same streaming loader as ITD mode (`.load_bam_data_streaming`).
2. **MAPQ filter** – reads below `min_mapq` are discarded.
3. **Soft‑clip extraction** – leading and trailing clips are extracted per
   CIGAR string.
4. **ALU alignment** – each clip ≥ `min_clip_len` bp is aligned against all
   loaded ALU consensus sequences via local Smith‑Waterman
   (`Biostrings::pairwiseAlignment`), in both forward and reverse‑complement
   orientation. The best‑scoring hit (normalised score ≥ `min_alu_score`) is
   retained.
5. **Poly‑A detection** – clips are scanned for poly‑A (sense) or poly‑T
   (antisense) runs ≥ 6 nt; length is reported as `PolyALength`.
6. **Insertion‑site clustering** – ALU‑positive clips are clustered by genomic
   position within `cluster_tolerance` bp using the same `.cluster_breakpoints`
   helper as the ITD engine.
7. **TSD detection** – for each cluster, the reference sequence upstream and
   downstream of the median insertion site is compared to identify the longest
   direct‑repeat prefix (≤ `max_tsd` bp, default 20).
8. **5′ truncation** – the alignment start position in the ALU consensus
   directly reports how many 5′ bases are absent (common in de novo insertions).
9. **Quality filters** – clusters failing `min_support`, `min_alu_score`,
   `vaf_threshold`, or `min_tsd` are discarded.
10. **Hotspot annotation** – insertion sites are overlapped with the bundled
    CSV, identical to ITD mode.
11. **Output** – TSV, optional HTML report (via the shared `TALOS_report.Rmd`
    template in ALU mode).

## Functions

| Function                     | Description |
|------------------------------|-------------|
| `talos()`                    | Simplified ITD/PTD entry point |
| `detect_itd()`               | Core ITD/PTD detection, clustering, and output |
| `talos_alu()`                | Simplified ALU insertion entry point |
| `detect_alu()`               | Core ALU detection, TSD/poly‑A characterisation, and output |
| `get_gene_config()`          | Loads gene reference from YAML |
| `annotate_hotspots()`        | Annotates events using a DBI connection or CSV file |
| `compute_hgvs_annotations()` | Maps genomic breakpoint to HGVS `c.` and `p.` notation |
| `plot_talos_report()`        | Publication‑ready PDF with coverage track, gene model, and arcs |
| `plot_talos_interactive()`   | plotly widget for interactive exploration |
| `talos_html_report()`        | Renders `inst/rmd/TALOS_report.Rmd`; supports ITD and ALU modes |
| `talos_batch()`              | Runs `talos()` over multiple BAM files and/or genes |
| `talos_batch_dir()`          | Discovers BAM files in a directory and runs `talos_batch()` |
| `save_offline_config()`      | Snapshots a resolved gene config to YAML for offline use |
| `detect_orientation()`       | Detects whether a duplication is tandem or inverted |
| `build_vcf_header()`         | Constructs VCF v4.2 header lines |
| `build_vcf_records()`        | Formats results as VCF data lines |
| `write_itd_vcf()`            | Writes a complete VCF file from TALOS results |

## Parameters

### `talos()` – ITD/PTD interface

| Parameter                  | Default       | Description |
|----------------------------|---------------|-------------|
| `bam_path`                 | –             | Path to an indexed BAM file |
| `gene`                     | –             | Gene symbol: `"FLT3"`, `"BCOR"`, `"KMT2A"`, or `"UBTF"` |
| `build`                    | `"hg19"`      | Genome build: `"hg19"` or `"hg38"` |
| `min_support`              | `10`          | Minimum bias‑corrected supporting reads |
| `min_size`                 | `15`          | Minimum duplication length (bp) |
| `max_correction`           | `2.0`         | Cap on the size‑bias correction factor |
| `plot`                     | `TRUE`        | Generate a PDF plot (requires `Gviz`) |
| `output_prefix`            | `"TALOS"`     | Base name for output files; `NULL` suppresses file writing |
| `output_folder`            | `"./results"` | Output directory |
| `sample_name`              | `NULL`        | Sample identifier; derived from BAM filename if `NULL` |
| `hotspot_db_path`          | `NULL`        | Path to a custom hotspot CSV |
| `yaml_path`                | (bundled)     | Path to gene config YAML |
| `padding`                  | `500`         | Base‑pair padding around targeted exons |
| `bsgenome`                 | `NULL`        | BSgenome object or package‑name string |
| `html_report`              | `TRUE`        | Write a self‑contained HTML report |
| `vaf_threshold`            | `0.01`        | Minimum allele frequency to report |
| `filter_intronic`          | `FALSE`       | Drop breakpoints not overlapping a targeted exon |
| `annotate_hotspots`        | `TRUE`        | Run hotspot annotation |
| `detect_orientation`       | `TRUE`        | Detect tandem vs inverted orientation |
| `compute_alignment_score`  | `TRUE`        | Alignment score between ITD and reference |
| `compute_support_bases`    | `TRUE`        | Count total soft‑clip bases from supporting reads |
| `compute_consistency`      | `TRUE`        | Fraction of supporting reads containing the ITD |
| `compute_itd_coverage`     | `TRUE`        | Per‑position coverage across the ITD sequence |
| `compute_coverage_drop`    | `TRUE`        | Fold‑change in coverage at the breakpoint |
| `compute_microhomology`    | `TRUE`        | Median microhomology length |
| `compute_repeat_entropy`   | `TRUE`        | Dinucleotide Shannon entropy around the breakpoint |
| `compute_discordant_ratio` | `TRUE`        | Discordant read‑pair ratio (triggers second BAM pass) |
| `compute_hgvs`             | `TRUE`        | Generate HGVS `c.` and `p.` notation |
| `add_config_to_report`     | `FALSE`       | Include gene configuration in reports |
| `max_pairwise_alignments`  | `30`          | Reads subsampled for pairwise consistency alignment |

### `talos_alu()` – ALU insertion interface

| Parameter                   | Default       | Description |
|-----------------------------|---------------|-------------|
| `bam_path`                  | –             | Path to an indexed BAM file |
| `gene`                      | –             | Gene symbol or ALU‑mode entry, e.g. `"KMT2A_ALU"` |
| `build`                     | `"hg19"`      | Genome build |
| `consensus_fa`              | `NULL`        | Path to custom ALU consensus FASTA; `NULL` uses bundled file |
| `min_support`               | `3`           | Minimum supporting reads per insertion cluster |
| `min_alu_score`             | `0.60`        | Minimum normalised Smith‑Waterman score (0–1) |
| `min_clip_len`              | `25`          | Minimum soft‑clip length to attempt ALU alignment (bp) |
| `min_mapq`                  | `20`          | Minimum MAPQ for reads used in depth calculation |
| `cluster_tolerance`         | `15`          | Clustering radius for insertion sites (bp) |
| `vaf_threshold`             | `0.01`        | Minimum allele frequency to report |
| `min_tsd`                   | `0`           | Minimum TSD length; `0` disables the filter |
| `output_prefix`             | `"TALOS_ALU"` | Base name for output files |
| `output_folder`             | `"./results"` | Output directory |
| `sample_name`               | `NULL`        | Sample identifier |
| `html_report`               | `TRUE`        | Write a self‑contained HTML report |
| `do_annotate_hotspots`      | `TRUE`        | Run hotspot annotation |
| `hotspot_db_path`           | `NULL`        | Path to a custom hotspot CSV |
| `yaml_path`                 | (bundled)     | Path to gene config YAML |
| `padding`                   | `500`         | Base‑pair padding around targeted exons |
| `bsgenome`                  | `NULL`        | BSgenome object or package‑name string |
| `max_reads_in_region`       | `200000`      | Hard cap on reads loaded from the BAM |
| `verbose`                   | `TRUE`        | Print progress messages |

## Output

### ITD / PTD columns (one row per event)

| Column               | Description |
|----------------------|-------------|
| `Sample`             | Sample identifier |
| `Gene`               | Gene symbol |
| `Genome`             | Genome build |
| `GenomicPosition`    | Breakpoint coordinate (1‑based) |
| `ITD_Sequence`       | Duplicated nucleotide sequence |
| `Length`             | Duplication length (bp) |
| `SupportingReads`    | Bias‑corrected mutant fragment count |
| `WildtypeReads`      | Wild‑type fragment count |
| `DepthAtBreakpoint`  | Total depth |
| `AlleleFrequency`    | Mutant allele frequency |
| `HGVS_cDNA`          | HGVS cDNA notation |
| `HGVS_Protein`       | HGVS protein notation |
| `StrandBias`         | Fraction of supporting reads on reverse strand |
| `MeanSupportMAPQ`    | Mean MAPQ of supporting reads |
| `Orientation`        | `"tandem"` or `"inverted"` |
| `Hotspot`            | `TRUE` if overlaps a known hotspot |
| `HotspotName`        | Matched hotspot name(s) |
| `Region`             | `"exonic"` or `"intronic"` |
| `ExonNumber`         | Exon number of the breakpoint (if exonic) |
| `AlignmentScore`     | Per‑base alignment score between ITD and reference (0–1) |
| `SupportConsistency` | % of supporting reads containing the ITD sequence |
| `ITDReadCoverage`    | % of ITD positions covered by at least one supporting read |
| `SequenceSource`     | `"observed"`, `"partial_read+ref"`, or `"ref_imputed"` |
| *(+ further metric columns as described in original docs)* | |

### ALU insertion columns (one row per insertion)

| Column              | Description |
|---------------------|-------------|
| `Sample`            | Sample identifier |
| `Gene`              | Gene symbol |
| `Genome`            | Genome build |
| `InsertionSite`     | Insertion coordinate (genomic, 1‑based) |
| `ALUSubtype`        | Best‑matching ALU subfamily (e.g. `"AluSx"`) |
| `ALUOrientation`    | `"sense"` or `"antisense"` relative to the gene strand |
| `TSD_Length`        | Target‑site duplication length (bp); `NA` if undetected |
| `TSD_Sequence`      | TSD nucleotide sequence |
| `PolyALength`       | Poly‑A (or poly‑T) tail length in the supporting clip (nt) |
| `ALU5pTruncation`   | Estimated 5′ truncation of the ALU copy (bp) |
| `ALU_Sequence`      | Longest supporting soft‑clip sequence (best clip per cluster) |
| `SupportingReads`   | Number of ALU‑positive soft‑clip reads in the cluster |
| `WildtypeReads`     | Wild‑type reads at the insertion site |
| `DepthAtBreakpoint` | Total depth at the insertion site |
| `AlleleFrequency`   | Insertion allele frequency |
| `MeanSupportMAPQ`   | Mean MAPQ of supporting reads |
| `StrandBias`        | Fraction of supporting reads on reverse strand |
| `Hotspot`           | `TRUE` if insertion overlaps a known hotspot |
| `HotspotName`       | Matched hotspot name(s) |
| `Region`            | `"exonic"` or `"intronic"` |
| `ExonNumber`        | Exon number (if insertion is exonic) |

### Output files

When `output_prefix` is provided, files are written to
`output_folder/<sample_name>/`:

- `<prefix>_<gene>_<timestamp>.tsv` – tab‑separated results
- `<prefix>_<gene>_<timestamp>.vcf` – VCF v4.2 (ITD mode only)
- `<sample>_<gene>.pdf` – PDF plot per sample/gene (ITD mode, when `plot = TRUE`)
- `TALOS_ALU_<sample>_<gene>_<timestamp>_report.html` – ALU HTML report
- `TALOS_report_<gene>_<timestamp>.html` – ITD summary HTML report

## Notes

- ALU detection sensitivity scales primarily with `min_clip_len` and
  `min_alu_score`. For low‑depth panels (< 500×), lower `min_support` to 2
  and increase `min_alu_score` to 0.65 to maintain specificity.
- `poly_a_len` reports the longest A‑run in any supporting clip for the
  cluster. A value ≥ 15 alongside a detected TSD strongly supports a genuine
  new retrotransposon insertion.
- `ALU5pTruncation > 200` indicates a highly truncated copy; these are
  common somatic insertions and are not filtered by default.
- An allele frequency below 1% (`vaf_threshold` default) is silently filtered.
  For somatic mosaicism screening, lower this to `0.005`.
- Intronic ITD/PTD breakpoints produce `NA` in the HGVS columns. Filter with
  `!is.na(HGVS_cDNA)` to retain only exonic calls.
- `compute_discordant_ratio = TRUE` triggers a second BAM pass. Set to `FALSE`
  when processing many samples in parallel to halve I/O.
- KMT2A uses `ptd_mode = FALSE` (k‑mer backward‑jump) combined with
  `convert_long_to_ptd = TRUE` and `max_itd_length = 30000`. Duplications
  > 30 kb are promoted to PTD records (`Length = 0`).
- The UCSC refGene TxDb is built on the first `use_db = TRUE` call and cached
  to `tools::R_user_dir("TALOS", "cache")`. To force a rebuild, delete the
  cache file.
- **`save_offline_config` key naming** – the output YAML key is
  `<GENE>_OFFLINE`; rename it back to the plain gene symbol before passing
  the file via `yaml_path`.

## Code review changelog

The following issues identified during code review were addressed in this
release:

| ID     | Location                       | Issue | Resolution |
|--------|--------------------------------|-------|------------|
| CR‑1   | `engine_candidates.R`, `engine_metrics.R` | `.extract_candidates_standard()` and `.compute_variant_metrics()` were 200–300 lines mixing concerns | Split into single‑responsibility helpers in `engine_alu.R` as a model; ITD engine refactor tracked separately |
| CR‑2   | `report.R` – `talos_html_report()` | `...` was forwarded directly to `rmarkdown::render()`, allowing unsafe param injection | `talos_html_report()` now takes explicit named params only |
| CR‑3   | `main.R`                       | KMT2A threshold relaxations were hard‑coded in `talos()` in addition to appearing in the YAML | Logic consolidated; YAML `gene_settings:` block is now the single source of truth |
| CR‑4   | `report.R`                     | Inline Rmd template was generated as an escaped string in a temp file — impossible to customise or version‑control | Extracted to `inst/rmd/TALOS_report.Rmd`; `talos_html_report()` renders it via `system.file()` |
| CR‑5   | `report.R`                     | `talos_html_report()` had no fallback when the template was missing | Added explicit error with path hint and development‑fallback path |
| CR‑6   | `engine_alu.R` (new)           | ALU consensus stored inline as string literals would be hard to extend | Consensus loaded from `inst/extdata/alu_consensus.fa` at runtime; custom path accepted via `consensus_fa=` |
| CR‑7   | `engine_alu.R` (new)           | TSD and poly‑A detection implemented as independent, testable functions | `.detect_alu_tsd()` and `.detect_poly_a()` are separate helpers |
| CR‑8   | throughout                     | `sapply` used for type‑uncertain returns | Replaced with `vapply` where the return type is known |

## License

GPL (≥ 3)

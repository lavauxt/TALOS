# ============================================================================
# Report assets embedded in R: stylesheet and client-side app of the
# self-contained HTML report (no inst/ files, no R Markdown). Both strings are
# written verbatim into every report. Raw strings need R >= 4.0; the text is
# ASCII-only (non-ASCII characters are escaped) so the file is locale-safe.
#
# Redesign (ECHO / SeqOne / Illumina variant-grid conventions):
#   - tiered columns: a "Clinical view" preset shows ~12 key columns; every
#     other column lives behind a column chooser with per-group toggles
#   - colour-coded column groups (Identity / Call / Support / Evidence /
#     Quality / Sequence) with header strips and a legend
#   - inline row detail: an expander reveals key stats, status chips and
#     per-call warnings under the row (hover-scan + click-to-expand)
#   - sticky header and sticky leading columns, compact/comfortable density
#   - status colours + icons (VAF chips, PTD tags, hotspot stars, warning
#     banners) so problems are visible while scanning
#   - quick filters (hotspot / flagged / ITD / PTD), text filter, CSV export
# ============================================================================

.TALOS_REPORT_CSS <- r"-----(
/* TALOS report -- design tokens, layout and components.
   Chrome (surface / ink / grid) swaps between light and dark; status colours
   (good / warn / danger / info / ptd) drive chips, banners and table
   highlights so problems can be spotted while scanning. */
:root {
  color-scheme: light;
  --page: #f9f9f7;
  --surface: #fcfcfb;
  --text: #0b0b0b;
  --text-2: #52514e;
  --muted: #898781;
  --grid: #e1e0d9;
  --axis: #c3c2b7;
  --ring: rgba(11, 11, 11, 0.10);
  --hover: #f0f2ee;
  --accent: #2a78d6;
  --accent-ink: #ffffff;
  --accent-soft: #e7f0fc;
  --good: #006300;   --good-bg: rgba(0, 99, 0, 0.08);
  --warn: #8a5a00;   --warn-bg: rgba(250, 178, 25, 0.12);
  --warn-bg-solid: #fcf3df;  /* opaque tint: sticky cells must not let scrolled content bleed through */
  --danger: #b3261e; --danger-bg: rgba(179, 38, 30, 0.09);
  --ptd: #6250d6;    --ptd-bg: rgba(98, 80, 214, 0.10);
  --radius: 10px;
  --font: system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  --gc-identity: #898781;
  --gc-call: #2a78d6;
  --gc-insertion: #2a78d6;
  --gc-support: #006300;
  --gc-evidence: #6250d6;
  --gc-qc: #8a5a00;
  --gc-sequence: #0e7490;
  --gc-internal: #898781;
  --gc-other: #898781;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
    --muted: #898781; --grid: #2c2c2a; --axis: #383835;
    --ring: rgba(255, 255, 255, 0.10); --hover: #232320;
    --accent: #3987e5; --accent-ink: #0b0b0b;
    --accent-soft: #1c2a3d;
    --good: #0ca30c;  --good-bg: rgba(12, 163, 12, 0.14);
    --warn: #fab219;  --warn-bg: rgba(250, 178, 25, 0.14);
    --warn-bg-solid: #392f19;
    --danger: #ff7b72; --danger-bg: rgba(255, 123, 114, 0.14);
    --ptd: #9085e9;   --ptd-bg: rgba(144, 133, 233, 0.16);
    --gc-call: #3987e5; --gc-insertion: #3987e5; --gc-support: #0ca30c;
    --gc-evidence: #9085e9; --gc-qc: #fab219; --gc-sequence: #22d3ee;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
  --muted: #898781; --grid: #2c2c2a; --axis: #383835;
  --ring: rgba(255, 255, 255, 0.10); --hover: #232320;
  --accent: #3987e5; --accent-ink: #0b0b0b;
  --accent-soft: #1c2a3d;
  --good: #0ca30c;  --good-bg: rgba(12, 163, 12, 0.14);
  --warn: #fab219;  --warn-bg: rgba(250, 178, 25, 0.14);
  --warn-bg-solid: #392f19;
  --danger: #ff7b72; --danger-bg: rgba(255, 123, 114, 0.14);
  --ptd: #9085e9;   --ptd-bg: rgba(144, 133, 233, 0.16);
  --gc-call: #3987e5; --gc-insertion: #3987e5; --gc-support: #0ca30c;
  --gc-evidence: #9085e9; --gc-qc: #fab219; --gc-sequence: #22d3ee;
}

* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0; background: var(--page); color: var(--text);
  font: 14px/1.5 var(--font); -webkit-font-smoothing: antialiased;
}
a { color: var(--accent); }
h1, h2, h3 { margin: 0; font-weight: 600; letter-spacing: -0.01em; }
h2 { font-size: 18px; } h3 { font-size: 14px; }
code, pre { font-family: var(--mono); font-size: 12px; }
pre { margin: 0; padding: 12px; overflow: auto; background: var(--page);
      border: 1px solid var(--ring); border-radius: 8px; max-height: 420px; }
.muted { color: var(--muted); } .small { font-size: 12px; }
.sr-only { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0 0 0 0); }
.ic { display: inline-flex; align-items: center; }

/* ---- Header & tabs ---- */
.top { background: var(--surface); border-bottom: 1px solid var(--ring); }
.top-inner { max-width: 1440px; margin: 0 auto; padding: 14px 16px 0; }
.title-row { display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }
.title-row h1 { font-size: 20px; }
.badge { display: inline-flex; align-items: center; gap: 6px; padding: 2px 10px;
         border-radius: 999px; background: var(--accent-soft); color: var(--text-2);
         font-size: 12px; border: 1px solid var(--ring); }
.spacer { flex: 1; }
.tabs { display: flex; gap: 2px; overflow-x: auto; scrollbar-width: thin; margin-top: 10px; }
.tab { appearance: none; border: 0; background: none; color: var(--text-2); font: inherit;
       padding: 10px 14px; cursor: pointer; border-bottom: 2px solid transparent;
       white-space: nowrap; border-radius: 6px 6px 0 0; }
.tab:hover { background: var(--accent-soft); color: var(--text); }
.tab[aria-selected="true"] { color: var(--text); border-bottom-color: var(--accent); font-weight: 600; }
.tab:focus-visible, button:focus-visible, select:focus-visible, input:focus-visible, summary:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 2px; }
main { max-width: 1440px; margin: 0 auto; padding: 16px; }
.panel[hidden] { display: none; }

/* ---- Cards, grids, controls ---- */
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius);
        padding: 14px; min-width: 0; }
.card > header { display: flex; align-items: baseline; gap: 8px; margin-bottom: 8px; flex-wrap: wrap; }
.card > header .hint { margin-left: auto; }
.grid { display: grid; gap: 12px; }
.g2 { grid-template-columns: repeat(2, minmax(0, 1fr)); align-items: start; }
.g3 { grid-template-columns: repeat(3, minmax(0, 1fr)); align-items: start; }
.stack { display: flex; flex-direction: column; gap: 12px; }
.toolbar { display: flex; gap: 12px; flex-wrap: wrap; align-items: flex-end; margin-bottom: 12px; }
.field { display: flex; flex-direction: column; gap: 3px; font-size: 12px; color: var(--text-2); min-width: 0; max-width: 100%; }
.field.inline { flex-direction: row; align-items: center; gap: 6px; }
select, input[type="text"], input[type="number"], input[type="search"] {
  font: inherit; color: var(--text); background: var(--surface); border: 1px solid var(--axis);
  border-radius: 6px; padding: 5px 8px; min-width: 120px; max-width: 280px; }
input[type="range"] { accent-color: var(--accent); }
button.btn { font: inherit; font-size: 13px; padding: 5px 12px; border-radius: 6px; cursor: pointer;
             border: 1px solid var(--axis); background: var(--surface); color: var(--text);
             display: inline-flex; align-items: center; gap: 6px; }
button.btn:hover { background: var(--accent-soft); }
button.btn.primary { background: var(--accent); border-color: var(--accent); color: var(--accent-ink); }
button.btn[disabled] { opacity: 0.45; cursor: not-allowed; }
.plot { width: 100%; min-height: 360px; }
.plot.tall { min-height: 560px; } .plot.short { min-height: 280px; }
.empty { padding: 28px 12px; text-align: center; color: var(--muted); }
.note { font-size: 12px; color: var(--muted); }

/* ---- Status banner (ECHO-style, cannot be missed) ---- */
.banner { display: flex; flex-wrap: wrap; gap: 8px 12px; align-items: center; padding: 10px 14px;
          border: 1px solid var(--warn); border-left: 4px solid var(--warn);
          background: var(--warn-bg); border-radius: 10px; margin-bottom: 12px; }
.banner.neutral { border-color: var(--axis); background: var(--surface); border-left-color: var(--accent); }
.banner .btitle { font-weight: 600; color: var(--warn); display: inline-flex; gap: 7px; align-items: center; }
.banner.neutral .btitle { color: var(--text-2); }
.banner .bchips { display: flex; flex-wrap: wrap; gap: 6px; }
button.bchip { appearance: none; font: inherit; font-size: 12px; padding: 2px 10px; border-radius: 999px;
               border: 1px solid var(--ring); background: var(--surface); color: var(--text-2); cursor: pointer; }
button.bchip:hover { border-color: var(--accent); color: var(--text); }

/* ---- Stat tiles ---- */
.tiles { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
.tile { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius); padding: 12px 14px; }
.tile .v { font-size: 26px; font-weight: 600; letter-spacing: -0.02em; }
.tile .l { color: var(--text-2); font-size: 12px; }
.tile.warn .v { color: var(--warn); }

/* ---- Legend (genome view categories) ---- */
.legend { display: flex; flex-wrap: wrap; gap: 4px 12px; margin-top: 8px; font-size: 12px; color: var(--text-2); }
.legend button { appearance: none; border: 0; background: none; font: inherit; color: inherit;
                 display: inline-flex; gap: 6px; align-items: center; cursor: pointer; padding: 2px 4px; border-radius: 4px; }
.legend button:hover { background: var(--accent-soft); }
.legend button[aria-pressed="false"] { opacity: 0.4; text-decoration: line-through; }
.swatch { width: 10px; height: 10px; border-radius: 3px; display: inline-block; }

/* ---- Tables ---- */
.tbl-bar { display: flex; gap: 10px; align-items: center; margin: 8px 0; flex-wrap: wrap; }
.pager { display: inline-flex; gap: 6px; align-items: center; margin-left: auto; }
.tbl-wrap { overflow: auto; border: 1px solid var(--ring); border-radius: 8px; max-height: 560px; }
table.tbl { border-collapse: separate; border-spacing: 0; width: 100%; font-size: 12.5px; }
.tbl th, .tbl td { padding: 6px 10px; text-align: left; white-space: nowrap; border-bottom: 1px solid var(--grid); }
.tbl thead th { position: sticky; top: 0; z-index: 3; background: var(--surface); cursor: pointer; user-select: none;
                font-weight: 600; color: var(--text-2); border-bottom: 1px solid var(--axis); }
.tbl thead th.no-sort { cursor: default; }
.tbl th[aria-sort="ascending"]::after { content: " \25B2"; font-size: 9px; }
.tbl th[aria-sort="descending"]::after { content: " \25BC"; font-size: 9px; }
.tbl td.num, .tbl th.num { text-align: right; font-variant-numeric: tabular-nums; }
.tbl tbody tr.clickable { cursor: pointer; }
.tbl tbody tr:hover td { background: var(--hover); }
/* column-group colour strips on the header */
.tbl thead th.g-identity  { border-top: 3px solid var(--gc-identity); }
.tbl thead th.g-call      { border-top: 3px solid var(--gc-call); }
.tbl thead th.g-insertion { border-top: 3px solid var(--gc-insertion); }
.tbl thead th.g-support   { border-top: 3px solid var(--gc-support); }
.tbl thead th.g-evidence  { border-top: 3px solid var(--gc-evidence); }
.tbl thead th.g-qc        { border-top: 3px solid var(--gc-qc); }
.tbl thead th.g-sequence  { border-top: 3px solid var(--gc-sequence); }
.tbl thead th.g-internal  { border-top: 3px solid var(--gc-internal); }
.tbl thead th.g-other     { border-top: 3px solid var(--gc-other); }
/* sticky leading columns (expander + first data columns) */
.tbl th.sticky, .tbl td.sticky { position: sticky; z-index: 2; background: var(--surface); }
.tbl thead th.sticky { z-index: 4; }
.tbl th.sticky.sticky-end, .tbl td.sticky.sticky-end {
  box-shadow: 4px 0 12px rgba(0, 0, 0, 0.30), 2px 0 0 var(--axis); }
.tbl thead th.sticky.sticky-end { box-shadow: 4px 0 12px rgba(0, 0, 0, 0.34), 2px 0 0 var(--axis); }
.tbl tbody tr:hover td.sticky { background: var(--hover); }
/* density */
.tbl.d-compact th, .tbl.d-compact td { padding: 2px 8px; font-size: 12px; }
/* row status tints (sticky cells use the opaque tint so scrolled content
   cannot bleed through) */
.tbl tbody tr.row-hotspot td, .tbl tbody tr.row-flagged td { background: var(--warn-bg); }
.tbl tbody tr.row-hotspot td.sticky, .tbl tbody tr.row-flagged td.sticky { background: var(--warn-bg-solid); }
.tbl tbody tr.row-hotspot:hover td, .tbl tbody tr.row-flagged:hover td { background: var(--hover); }
.tbl tbody tr.row-hotspot td:first-child, .tbl tbody tr.row-flagged td:first-child { box-shadow: inset 3px 0 0 var(--warn); }
.tbl tbody tr.row-ptd td:first-child { box-shadow: inset 3px 0 0 var(--ptd); }
/* expander column */
.tbl td.exp-c { width: 26px; padding: 0 4px; text-align: center; }
.tbl thead th.exp-h { border-top: 2px solid var(--axis); }
button.exp { appearance: none; border: 0; background: none; color: var(--muted); cursor: pointer;
             padding: 2px; border-radius: 4px; display: inline-flex; }
button.exp:hover { background: var(--accent-soft); color: var(--text); }
button.exp[aria-expanded="true"] { color: var(--accent); }
button.exp[aria-expanded="true"] .ic { transform: rotate(90deg); }
button.exp .ic { transition: transform 0.12s ease; }
/* inline detail row: sticky to the left edge so it stays visible while the
   wide table is scrolled horizontally */
.tbl tr.detail-row td { position: sticky; left: 0; z-index: 1; background: var(--page);
                        border-bottom: 1px solid var(--axis); padding: 12px 14px; white-space: normal; }
.detail-grid { display: grid; grid-template-columns: minmax(240px, 1fr) minmax(280px, 1.5fr) minmax(190px, 0.8fr); gap: 16px; }
.mini-stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(86px, 1fr)); gap: 8px; }
.mini-stat { background: var(--surface); border: 1px solid var(--ring); border-radius: 8px; padding: 6px 10px; min-width: 0; }
.mini-stat .v { font-size: 16px; font-weight: 600; letter-spacing: -0.01em; overflow: hidden; text-overflow: ellipsis; }
.mini-stat .v.good { color: var(--good); } .mini-stat .v.warn { color: var(--warn); } .mini-stat .v.danger { color: var(--danger); }
.mini-stat .l { color: var(--muted); font-size: 11px; }
.warn-list { display: flex; flex-direction: column; gap: 6px; font-size: 12.5px; margin-top: 8px; }
.warn-item { display: flex; gap: 7px; align-items: flex-start; color: var(--text-2); }
.warn-item .ic { flex: none; margin-top: 1px; }
.warn-item.sev-warn { color: var(--warn); }
.warn-item.sev-danger { color: var(--danger); font-weight: 600; }
.warn-item.sev-good { color: var(--good); }

/* ---- Status tags / chips ---- */
.tag { display: inline-flex; align-items: center; gap: 4px; padding: 1px 8px; border-radius: 999px;
       font-size: 11.5px; font-weight: 600; border: 1px solid var(--ring); background: var(--page);
       color: var(--text-2); white-space: nowrap; max-width: 100%; overflow: hidden; text-overflow: ellipsis; }
.tag.ok { color: var(--good); border-color: var(--good); background: var(--good-bg); }
.tag.info { color: var(--accent); border-color: var(--accent); background: var(--accent-soft); }
.tag.warn { color: var(--warn); border-color: var(--warn); background: var(--warn-bg); }
.tag.danger { color: var(--danger); border-color: var(--danger); background: var(--danger-bg); }
.tag.ptd { color: var(--ptd); border-color: var(--ptd); background: var(--ptd-bg); }
.tag.mut { color: var(--muted); }
.cell-warn { color: var(--warn); font-weight: 600; }
.cell-muted { color: var(--muted); }
.cell-hotspot { color: var(--warn); font-weight: 600; }
.cell-support { display: inline-flex; align-items: center; gap: 6px; }
.cell-support .bar { display: inline-block; height: 8px; background: var(--accent); border-radius: 2px; min-width: 2px; }
.mono { font-family: var(--mono); font-size: 11.5px; }
.chips { display: flex; flex-wrap: wrap; gap: 6px; }
.chip { font-size: 12px; padding: 2px 8px; border-radius: 999px; border: 1px solid var(--ring);
        background: var(--page); color: var(--text-2); }

/* ---- Quick filters, density, column chooser ---- */
.qchips { display: flex; flex-wrap: wrap; gap: 6px; }
button.qchip { appearance: none; font: inherit; font-size: 12px; padding: 3px 10px; border-radius: 999px; cursor: pointer;
               border: 1px solid var(--axis); background: var(--surface); color: var(--text-2);
               display: inline-flex; align-items: center; gap: 5px; }
button.qchip:hover { background: var(--accent-soft); }
button.qchip[aria-pressed="true"] { background: var(--accent); border-color: var(--accent); color: var(--accent-ink); font-weight: 600; }
.seg { display: inline-flex; border: 1px solid var(--axis); border-radius: 6px; overflow: hidden; }
button.seg-btn { appearance: none; border: 0; background: var(--surface); color: var(--text-2); font: inherit;
                 font-size: 12px; padding: 4px 10px; cursor: pointer; }
button.seg-btn + button.seg-btn { border-left: 1px solid var(--axis); }
button.seg-btn[aria-pressed="true"] { background: var(--accent-soft); color: var(--text); font-weight: 600; }
.glegend { display: flex; flex-wrap: wrap; gap: 4px 14px; margin: 8px 2px; font-size: 11.5px; color: var(--text-2); }
.glegend .gk { display: inline-flex; align-items: center; gap: 5px; }
.gdot { width: 8px; height: 8px; border-radius: 2px; display: inline-block; flex: none; }
details.colmenu { position: relative; }
details.colmenu > summary { list-style: none; display: inline-flex; align-items: center; gap: 6px; }
details.colmenu > summary::-webkit-details-marker { display: none; }
.colmenu-panel { position: absolute; right: 0; top: calc(100% + 6px); z-index: 40; width: 360px; max-height: min(480px, 70vh);
                 overflow: auto; background: var(--surface); border: 1px solid var(--axis); border-radius: 10px;
                 padding: 10px 12px; box-shadow: 0 10px 28px rgba(0, 0, 0, 0.16); scrollbar-width: thin; }
.colmenu .hintline { font-size: 11.5px; color: var(--muted); margin: 6px 0 2px; }
:root[data-theme="dark"] .colmenu-panel { box-shadow: 0 10px 28px rgba(0, 0, 0, 0.5); }
.colmenu fieldset { border: 0; margin: 4px 0 10px; padding: 0; min-width: 0; }
.colmenu legend { font-size: 12px; font-weight: 600; color: var(--text-2); display: flex; align-items: center; gap: 6px; padding: 0; }
.colmenu .cols { display: grid; grid-template-columns: 1fr 1fr; gap: 1px 10px; }
.colmenu .cols label { font-size: 12px; display: flex; gap: 6px; align-items: center; white-space: nowrap;
                       overflow: hidden; text-overflow: ellipsis; max-width: 100%; }
.colmenu .cols label:hover { color: var(--text); }

/* ---- Key/value, variant cards ---- */
.kv { display: grid; grid-template-columns: max-content 1fr; gap: 3px 14px; font-size: 12.5px; }
.kv dt { color: var(--text-2); } .kv dd { margin: 0; overflow-wrap: anywhere; }
.kv.wide { grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 4px 18px; }
.kv.wide > div { display: grid; grid-template-columns: 48% 52%; gap: 8px; border-bottom: 1px dashed var(--grid); padding: 2px 0; }
.kv.wide dt, .kv.wide dd { font-size: 12.5px; }
.group-title { margin: 18px 0 8px; font-size: 15px; color: var(--text-2); }
.vcard { scroll-margin-top: 12px; }
.vcard.hl { outline: 2px solid var(--accent); outline-offset: 2px; }
.vcard h3 { font-size: 14px; }
.badge.flag { background: transparent; border-color: var(--warn); color: var(--warn); font-weight: 600; }
.badge-hotspot { background: var(--warn-bg); color: var(--warn); border: 1px solid var(--warn);
                 padding: 1px 6px; border-radius: 4px; font-size: 11px; font-weight: 600;
                 display: inline-flex; align-items: center; gap: 4px; }
.badge-exon { background: var(--accent-soft); color: var(--text); border: 1px solid var(--ring);
              padding: 1px 6px; border-radius: 4px; font-size: 11px; }
details.seq summary { cursor: pointer; color: var(--accent); font-size: 12.5px; margin-top: 8px; }
details.seq pre { white-space: pre-wrap; word-break: break-all; margin-top: 6px; }

@media (max-width: 1000px) { .g2, .g3 { grid-template-columns: minmax(0, 1fr); } }
@media (max-width: 900px) { .detail-grid { grid-template-columns: minmax(0, 1fr); } }
@media (max-width: 720px) { main { padding: 12px; } }
@media print {
  .tabs, .toolbar, .btn, .qchips, .colmenu, .pager, .tbl-bar { display: none !important; }
  .panel[hidden] { display: block !important; }
  .card { break-inside: avoid; }
  .tbl-wrap { max-height: none; overflow: visible; }
  .js-plotly-plot .modebar { display: none !important; }
  main { padding: 0; }
}
@media (forced-colors: active) { .bar-fill, .swatch, .gdot { forced-color-adjust: none; } }
)-----"

.TALOS_REPORT_JS <- r"-----(
/* TALOS report -- client-side app (SPARKS style).
   Reads the JSON blob in #talos-data, builds one lazily-rendered panel per
   tab and draws figures with plotly.js. Tables and details work without
   plotly.js.

   Results-table system:
   - tiered columns: "Clinical view" shows the ~12 key columns; everything
     else sits behind the column chooser (per-group toggles + presets)
   - colour-coded column groups with header strips and a legend
   - inline row detail (expander) with key stats, status chips and warnings
   - sticky header + sticky leading columns, density switch, page size
   - quick filters (hotspot / flagged / ITD / PTD) + text filter + CSV export */
(function () {
  'use strict';

  // ------------------------------------------------------------------ data & utils
  var DATA = JSON.parse(document.getElementById('talos-data').textContent);
  var $ = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

  function h(tag, attrs) {
    var el = document.createElement(tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      var v = attrs[k];
      if (v === null || v === undefined || v === false) return;
      if (k === 'class') el.className = v;
      else if (k === 'text') el.textContent = v;
      else if (k.slice(0, 2) === 'on' && typeof v === 'function') el.addEventListener(k.slice(2), v);
      else if (k === 'style' && typeof v === 'object') Object.assign(el.style, v);
      else el.setAttribute(k, v === true ? '' : v);
    });
    for (var i = 2; i < arguments.length; i++) append(el, arguments[i]);
    return el;
  }
  function append(el, kid) {
    if (kid === null || kid === undefined || kid === false) return;
    if (Array.isArray(kid)) kid.forEach(function (k) { append(el, k); });
    else if (kid.nodeType) el.appendChild(kid);
    else el.appendChild(document.createTextNode(String(kid)));
  }
  function clear(el) { while (el.firstChild) el.removeChild(el.firstChild); return el; }
  function arr(x) { return x === null || x === undefined ? [] : (Array.isArray(x) ? x : [x]); }
  function range(n) { var a = new Array(n); for (var i = 0; i < n; i++) a[i] = i; return a; }
  function uniq(a) { return a.filter(function (v, i) { return a.indexOf(v) === i; }); }
  function fmt(v) {
    if (v === null || v === undefined || (typeof v === 'number' && !isFinite(v))) return '';
    if (typeof v !== 'number') return String(v);
    if (Number.isInteger(v)) return String(v);
    if (v !== 0 && Math.abs(v) < 1e-3) return v.toExponential(2);
    return String(Number(v.toPrecision(4)));
  }
  function fmtInt(v) {
    if (typeof v !== 'number' || !isFinite(v)) return fmt(v);
    var neg = v < 0 ? '-' : '';
    return neg + String(Math.round(Math.abs(v))).replace(/\B(?=(\d{3})+(?!\d))/g, ',');
  }
  function fmtPct(v, dig) {
    return typeof v === 'number' && isFinite(v) ? (v * 100).toFixed(dig == null ? 1 : dig) + '%' : 'N/A';
  }
  function median(a) {
    a = a.filter(function (x) { return typeof x === 'number' && isFinite(x); }).sort(function (p, q) { return p - q; });
    if (!a.length) return NaN;
    var m = a.length >> 1;
    return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
  }
  function download(name, text) {
    var blob = new Blob([text], { type: 'text/csv;charset=utf-8' });
    var a = h('a', { href: URL.createObjectURL(blob), download: name });
    document.body.appendChild(a); a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 0);
  }
  function csvField(v) {
    var s = v === null || v === undefined ? '' : String(v);
    return /[",\n\r\t]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  }
  function store(key, val) {
    try { if (val === undefined) return window.localStorage.getItem(key); window.localStorage.setItem(key, val); } catch (e) { /* storage may be unavailable */ }
    return null;
  }
  var T = DATA.table, NROW = T.n, HIDDEN = arr(DATA.hidden);
  var M = DATA.meta, VARIANTS = arr(DATA.variants), PANELS = arr(DATA.panels);
  var isAlu = M.mode === 'alu';
  function col(name) { return T.data[name] ? arr(T.data[name]) : []; }
  function hasCol(c) { return T.cols.indexOf(c) >= 0; }
  function cellRaw(c, i) { var a = col(c); return a.length ? a[i] : null; }
  function numv(c, i) { var v = cellRaw(c, i); return typeof v === 'number' && isFinite(v) ? v : null; }

  // ------------------------------------------------------------------ icons
  var ICONS = {
    chevron: '<polyline points="9 6 15 12 9 18"/>',
    warn: '<path d="M10.29 3.86L1.82 18a2 2 0 0 0 1.71 3h16.94a2 2 0 0 0 1.71-3L13.71 3.86a2 2 0 0 0-3.42 0z"/><line x1="12" y1="9" x2="12" y2="13"/><line x1="12" y1="17" x2="12.01" y2="17"/>',
    star: '<polygon points="12 2 15.09 8.26 22 9.27 17 14.14 18.18 21.02 12 17.77 5.82 21.02 7 14.14 2 9.27 8.91 8.26 12 2"/>',
    download: '<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" y1="15" x2="12" y2="3"/>',
    columns: '<path d="M12 3h7a2 2 0 0 1 2 2v14a2 2 0 0 1-2 2h-7"/><path d="M12 3H5a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h7"/><line x1="12" y1="3" x2="12" y2="21"/>',
    info: '<circle cx="12" cy="12" r="10"/><line x1="12" y1="16" x2="12" y2="12"/><line x1="12" y1="8" x2="12.01" y2="8"/>',
    jump: '<path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"/><polyline points="15 3 21 3 21 9"/><line x1="10" y1="14" x2="21" y2="3"/>',
    printer: '<polyline points="6 9 6 2 18 2 18 9"/><path d="M6 18H4a2 2 0 0 1-2-2v-5a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v5a2 2 0 0 1-2 2h-2"/><rect x="6" y="14" width="12" height="8"/>',
    eye: '<path d="M1 12s4-8 11-8 11 8 11 8-4 8-11 8-11-8-11-8z"/><circle cx="12" cy="12" r="3"/>',
    table: '<rect x="3" y="3" width="18" height="18" rx="2"/><line x1="3" y1="9" x2="21" y2="9"/><line x1="3" y1="15" x2="21" y2="15"/>',
    check: '<polyline points="20 6 9 17 4 12"/>'
  };
  function icon(name, size) {
    return '<svg viewBox="0 0 24 24" width="' + (size || 14) + '" height="' + (size || 14) +
      '" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false">' + ICONS[name] + '</svg>';
  }
  function iconEl(name, size) { var s = h('span', { class: 'ic' }); s.innerHTML = icon(name, size); return s; }

  // ------------------------------------------------------------------ theme & palette
  var CHROME = {
    light: { text: '#0b0b0b', secondary: '#52514e', muted: '#898781', grid: '#e1e0d9', axis: '#c3c2b7', surface: '#fcfcfb' },
    dark: { text: '#f0efec', secondary: '#c3c2b7', muted: '#898781', grid: '#2c2c2a', axis: '#383835', surface: '#1a1a19' }
  };
  var CAT = {
    light: ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948'],
    dark: ['#3987e5', '#d95926', '#199e70', '#c98500', '#d55181', '#008300', '#9085e9', '#e66767']
  };
  var EXON = { light: '#1baf7a', dark: '#199e70' };
  var FONT = "system-ui, -apple-system, 'Segoe UI', sans-serif";
  var themePref = store('talos-theme') || 'auto';
  function resolved() {
    if (themePref === 'light' || themePref === 'dark') return themePref;
    return window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
  }
  function applyThemeAttr() {
    if (themePref === 'auto') document.documentElement.removeAttribute('data-theme');
    else document.documentElement.setAttribute('data-theme', themePref);
  }
  function color(i, t) { return CAT[t][i % CAT[t].length]; }
  function alpha(hex, a) {
    var n = [1, 3, 5].map(function (i) { return parseInt(hex.slice(i, i + 2), 16); });
    return 'rgba(' + n.join(',') + ',' + a + ')';
  }

  // ------------------------------------------------------------------ plotly helpers
  var registry = [];
  var CONFIG = { displaylogo: false, responsive: true, modeBarButtonsToRemove: ['autoScale2d', 'lasso2d', 'select2d'] };
  function hasPlotly() { return typeof window.Plotly !== 'undefined'; }
  function axisStyle(t, extra) {
    var c = CHROME[t];
    return Object.assign({ gridcolor: c.grid, zerolinecolor: c.axis, linecolor: c.axis, tickfont: { color: c.muted },
      title: { font: { color: c.secondary } }, automargin: true }, extra || {});
  }
  function baseLayout(t, o) {
    var c = CHROME[t];
    return Object.assign({ paper_bgcolor: 'rgba(0,0,0,0)', plot_bgcolor: 'rgba(0,0,0,0)',
      font: { family: FONT, size: 12, color: c.secondary }, hoverlabel: { font: { family: FONT } },
      margin: { l: 64, r: 16, t: 16, b: 48 }, showlegend: false }, o || {});
  }
  // plotly.js is normally embedded; if it is missing or failed to load, the
  // CDN copy is fetched once and every pending figure is drawn when it lands.
  var PLOTLY_CDN = 'https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js';
  var NO_PLOTLY = 'Interactive plots need plotly.js: the embedded copy is missing and the CDN is unreachable. Tables and details still work.';
  var plotlyState = hasPlotly() ? 'ready' : 'missing';
  function ensurePlotly() {
    if (plotlyState !== 'missing') return;
    plotlyState = 'loading';
    var sc = document.createElement('script');
    sc.src = PLOTLY_CDN; sc.async = true;
    sc.onload = function () { plotlyState = hasPlotly() ? 'ready' : 'failed'; flushPending(); };
    sc.onerror = function () { plotlyState = 'failed'; flushPending(); };
    document.head.appendChild(sc);
  }
  function flushPending() {
    registry.forEach(function (d) {
      if (!d._render) return;
      if (plotlyState !== 'ready') { clear(d).appendChild(emptyNote(NO_PLOTLY)); return; }
      clear(d);
      if (d.offsetParent !== null) d._render(); else d._dirty = true;
    });
  }
  function plot(div, build) {
    div._render = function () {
      if (!hasPlotly()) return;
      var s = build(resolved());
      window.Plotly.react(div, s.data, s.layout, Object.assign({}, CONFIG, s.config || {}));
      div._dirty = false;
    };
    if (registry.indexOf(div) < 0) registry.push(div);
    if (hasPlotly()) { div._render(); return div; }
    clear(div).appendChild(emptyNote(plotlyState === 'failed' ? NO_PLOTLY : 'Loading plotly.js\u2026'));
    ensurePlotly();
    return div;
  }
  // Plotly attaches .on() at first render, which may come later (CDN fallback);
  // retries are capped so a missing library cannot spin forever.
  function onPlotClick(div, fn) {
    if (div._clickBound) return;
    if (div.on) { div._clickBound = true; div.on('plotly_click', fn); return; }
    if ((div._clickTries = (div._clickTries || 0) + 1) > 40) return;
    setTimeout(function () { onPlotClick(div, fn); }, 250);
  }
  function rerenderAll() {
    registry = registry.filter(function (d) { return document.body.contains(d) && d._render; });
    if (!hasPlotly()) return;
    registry.forEach(function (d) { if (d.offsetParent !== null) d._render(); else d._dirty = true; });
  }
  function refreshVisible(root) {
    $$('.plot', root).forEach(function (d) {
      if (!d._render) return;
      if (!hasPlotly()) return;
      if (d._dirty) d._render(); else if (d.data) window.Plotly.Plots.resize(d);
    });
    STICKY_TABLES.forEach(function (v) { v.applySticky(); });
  }
  function plotBox(cls) { return h('div', { class: 'plot ' + (cls || '') }); }
  function emptyNote(text) { return h('div', { class: 'empty', text: text }); }
  function card(title, hint, body, attrs) {
    return h('section', Object.assign({ class: 'card' }, attrs || {}),
      h('header', null, h('h3', { text: title }), hint ? h('span', { class: 'muted small hint', text: hint }) : null), body);
  }
  function kv(keys, values, wide) {
    var dl = h('dl', { class: 'kv' + (wide ? ' wide' : '') });
    arr(keys).forEach(function (k, i) {
      var v = arr(values)[i];
      if (wide) dl.appendChild(h('div', null, h('dt', { text: k }), h('dd', { text: v })));
      else { dl.appendChild(h('dt', { text: k })); dl.appendChild(h('dd', { text: v })); }
    });
    return dl;
  }

  // ------------------------------------------------------------------ column docs & groups
  var HELP = {
    Sample: 'Sample identifier (derived from the BAM file name unless set explicitly).',
    Gene: 'Target gene analysed.',
    Genome: 'Reference genome build of the run.',
    GenomeBuild: 'Reference genome build of the gene configuration.',
    TranscriptRef: 'Transcript used for HGVS notation ("none" when unavailable).',
    GenomicPosition: 'Leftmost genomic base of the duplicated segment (first copy).',
    InsertionSite: 'Median genomic position of the ALU insertion breakpoint cluster.',
    Length: 'Duplicated segment length (bp). 0 = PTD call from soft-clips only; values above ~1 kb are PTD-scale events.',
    DupExons: 'Duplicated exons in transcript numbering, from the exon-graph reconstruction.',
    ExonNumber: 'Exon containing the breakpoint, when exonic.',
    Region: 'exonic when the breakpoint falls inside a targeted exon, intronic otherwise.',
    HGVS_cDNA: 'HGVS c. notation of the duplication, relative to the configured transcript.',
    HGVS_Protein: 'HGVS p. notation (frameshifts are reported as p.XxxNNfsTer?).',
    Orientation: 'tandem = duplicated copy in the same orientation; inverted = reverse complemented.',
    CandidateTypes: 'How the call was found: bwd_jump / missing_kmer (k-mer walk), anchor_ext (soft-clip re-anchoring), cigar_ins, ptd_clip, ptd_graph (exon graph).',
    Hotspot: 'Overlaps a known hotspot from the annotation database.',
    HotspotName: 'Name of the matched hotspot.',
    AlleleFrequency: 'Size-bias-corrected variant allele fraction (supporting reads / depth at breakpoint).',
    SupportingReads: 'Supporting reads after size-bias correction.',
    WildtypeReads: 'Wild-type reads spanning the breakpoint (primary alignments, MAPQ-filtered).',
    DepthAtBreakpoint: 'Corrected support + wild-type depth at the breakpoint.',
    StrandBias: 'Fraction of supporting reads on the reverse strand; 0.25-0.75 is balanced.',
    MeanSupportMAPQ: 'Mean mapping quality of the supporting reads.',
    BreakpointSpread: 'Range of breakpoint positions within the supporting cluster (bp).',
    SoftclipFraction: 'Fraction of supporting reads detected through a soft-clip.',
    UniqueBreakpoints: 'Number of distinct breakpoint positions in the cluster.',
    CoverageDrop: 'Fold change of mean depth across the breakpoint vs flanks (window = coverage_window).',
    LocalCoverage: 'Minimum depth within +/-10 bp of the breakpoint.',
    SpanMinCoverage: 'Minimum depth over the whole duplicated span.',
    SpanMeanCoverage: 'Mean depth over the whole duplicated span.',
    SpanDepthFoldChange: 'Mean span depth / mean flanking depth. A true heterozygous duplication is roughly (1+VAF)x; very high values suggest a collapsed, multi-mapping region.',
    MedianMicrohomology: 'Median microhomology (bp) between soft-clip and reference at the junction.',
    DiscordantRatio: 'Fraction of everted / oversized pairs among pairs spanning the breakpoint.',
    RepeatEntropy: 'Dinucleotide Shannon entropy around the breakpoint (low = repetitive sequence).',
    SequenceImputed: 'TRUE when the ITD sequence was filled from the reference rather than observed in reads.',
    SequencePartial: 'TRUE when only part of the sequence was observed, the rest imputed.',
    SequenceSource: 'observed = read-derived; partial_read+ref = mixed; ref_imputed = reference-derived; merged_intervals = merged PTD.',
    RefMatch_Observed: '% identity between the observed (read-derived) sequence and the reference segment.',
    RefMatch_Total: '% identity between the full reconstructed sequence and the reference segment.',
    ITDReadCoverage: '% of the ITD sequence covered by supporting reads.',
    ITDCoverageRLE: 'Run-length encoding of per-base ITD read coverage (value:length,...).',
    SupportConsistency: '% of supporting reads aligning to the ITD sequence at >= 90% identity (informational).',
    AlignmentScore: 'Fraction of matching bases between the ITD sequence and the reference segment (0-1).',
    TotalSupportBases: 'Total soft-clip bases contributed by the supporting reads.',
    LengthPE: 'ITD length estimated from paired-end insert sizes (spanning pairs).',
    LengthPE_NSpanning: 'Number of spanning pairs used for the PE length estimate.',
    LengthExt: 'Anchor-extension length estimate (mode of per-read sizes).',
    PESoftclipSupport: 'Paired-end support among reads with a soft-clip.',
    PESoftclipEventPairs: 'Pairs whose insert size matches background + event length.',
    PESoftclipLongPairs: 'Pairs with an insert size >= 75% of the event length.',
    PEOrientationFR: 'FR-oriented pairs among soft-clip-supported pairs (normal).',
    PEOrientationRF: 'RF-oriented pairs (everted - tandem duplication signature).',
    PEOrientationFF: 'FF-oriented pairs (unusual).',
    PEOrientationRR: 'RR-oriented pairs (unusual).',
    PEOrientationOther: 'Pairs with another orientation.',
    PEOrientationDominant: 'Most frequent pair orientation.',
    LeftSoftclipCount: 'Supporting reads with a leading soft-clip.',
    RightSoftclipCount: 'Supporting reads with a trailing soft-clip.',
    BothSoftclipCount: 'Supporting reads soft-clipped on both sides.',
    LeftSoftclipPctSupport: '% of supporting reads with a leading soft-clip.',
    RightSoftclipPctSupport: '% of supporting reads with a trailing soft-clip.',
    LeftSoftclipPctWT: '% of wild-type reads with a leading soft-clip.',
    RightSoftclipPctWT: '% of wild-type reads with a trailing soft-clip.',
    ArtifactSuspect: 'Fully reference-imputed sequence with ~100% reference match - likely an alignment artifact.',
    MergedFromN: 'Number of duplication rows merged into this interval.',
    MergeFragmentWarning: 'Merged interval pools more fragments than plausibly come from one true PTD; review before reporting.',
    ITD_Sequence: 'Duplicated sequence (observed, partial or imputed - see SequenceSource).',
    ALUSubtype: 'Best-matching ALU consensus family (e.g. AluY, AluSx, AluJb).',
    ALUOrientation: 'sense = same orientation as the consensus; antisense = reverse complement.',
    TSD_Length: 'Target-site duplication length (bp); 5-20 bp supports a real insertion.',
    TSD_Sequence: 'Target-site duplication sequence.',
    PolyALength: 'Longest A/T run in the soft-clip (nt) - 3\u2032 poly-A tail of the insertion.',
    ALU5pTruncation: 'Bases missing from the 5\u2032 end of the ALU copy (truncated insertions are common).',
    ALU_RawScore: 'Normalised alignment score of the best soft-clip against the ALU consensus (0-1).',
    ALU_Sequence: 'Best-matching soft-clip sequence.'
  };

  var GROUPS_ITD = [
    { id: 'identity', label: 'Identity', cols: ['Sample', 'Gene', 'GenomeBuild', 'TranscriptRef'] },
    { id: 'call', label: 'Call', cols: ['GenomicPosition', 'Length', 'DupExons', 'ExonNumber', 'Region',
      'HGVS_cDNA', 'HGVS_Protein', 'Orientation', 'CandidateTypes', 'Hotspot', 'HotspotName'] },
    { id: 'support', label: 'Support', cols: ['AlleleFrequency', 'SupportingReads', 'WildtypeReads',
      'DepthAtBreakpoint', 'StrandBias', 'MeanSupportMAPQ'] },
    { id: 'evidence', label: 'Evidence', cols: ['LengthPE', 'LengthPE_NSpanning', 'LengthExt', 'PESoftclipSupport',
      'PESoftclipEventPairs', 'PESoftclipLongPairs', 'PEOrientationFR', 'PEOrientationRF', 'PEOrientationFF',
      'PEOrientationRR', 'PEOrientationOther', 'PEOrientationDominant', 'LeftSoftclipCount', 'RightSoftclipCount',
      'BothSoftclipCount', 'LeftSoftclipPctSupport', 'RightSoftclipPctSupport', 'LeftSoftclipPctWT',
      'RightSoftclipPctWT', 'TotalSupportBases', 'SoftclipFraction', 'BreakpointSpread', 'UniqueBreakpoints'] },
    { id: 'qc', label: 'Quality', cols: ['CoverageDrop', 'LocalCoverage', 'SpanMinCoverage', 'SpanMeanCoverage',
      'SpanDepthFoldChange', 'MedianMicrohomology', 'DiscordantRatio', 'RepeatEntropy', 'RefMatch_Observed',
      'RefMatch_Total', 'ITDReadCoverage', 'SupportConsistency', 'AlignmentScore', 'SequenceImputed',
      'SequencePartial', 'SequenceSource', 'ArtifactSuspect', 'MergeFragmentWarning', 'MergedFromN'] },
    { id: 'sequence', label: 'Sequence', cols: ['ITD_Sequence', 'ITDCoverageRLE'] }
  ];
  var KEY_ITD = ['Sample', 'Gene', 'GenomicPosition', 'Length', 'AlleleFrequency', 'HGVS_cDNA', 'HGVS_Protein',
    'SupportingReads', 'DepthAtBreakpoint', 'Hotspot', 'Region', 'GenomeBuild'];
  var GROUPS_ALU = [
    { id: 'identity', label: 'Identity', cols: ['Sample', 'Gene', 'GenomeBuild'] },
    { id: 'insertion', label: 'Insertion', cols: ['InsertionSite', 'ALUSubtype', 'ALUOrientation', 'TSD_Length',
      'TSD_Sequence', 'PolyALength', 'ALU5pTruncation', 'Region', 'ExonNumber', 'Hotspot', 'HotspotName'] },
    { id: 'support', label: 'Support', cols: ['AlleleFrequency', 'SupportingReads', 'WildtypeReads',
      'DepthAtBreakpoint', 'MeanSupportMAPQ', 'StrandBias'] },
    { id: 'evidence', label: 'Evidence', cols: ['ALU_RawScore', 'ALU_Sequence'] }
  ];
  var KEY_ALU = ['Sample', 'Gene', 'InsertionSite', 'ALUSubtype', 'ALUOrientation', 'AlleleFrequency',
    'SupportingReads', 'DepthAtBreakpoint', 'TSD_Length', 'PolyALength', 'Hotspot', 'Region'];
  var GROUPS = isAlu ? GROUPS_ALU : GROUPS_ITD;
  var KEY_COLS = isAlu ? KEY_ALU : KEY_ITD;

  function groupOf(c) {
    if (HIDDEN.indexOf(c) >= 0) return 'internal';
    for (var g = 0; g < GROUPS.length; g++) if (GROUPS[g].cols.indexOf(c) >= 0) return GROUPS[g].id;
    return 'other';
  }
  function groupLabel(id) {
    if (id === 'internal') return 'Internal';
    if (id === 'other') return 'Other';
    for (var g = 0; g < GROUPS.length; g++) if (GROUPS[g].id === id) return GROUPS[g].label;
    return 'Other';
  }
  // Column order: key columns first, then group order, then anything else.
  function orderCols(sel) {
    var out = [], seen = {};
    function push(list) {
      list.forEach(function (c) {
        if (!seen[c] && sel[c] && hasCol(c)) { seen[c] = 1; out.push(c); }
      });
    }
    push(KEY_COLS);
    GROUPS.forEach(function (g) { push(g.cols); });
    T.cols.forEach(function (c) { if (!seen[c] && sel[c]) { seen[c] = 1; out.push(c); } });
    return out;
  }

  // Column selection state (persisted per mode).
  var colStore = 'talos-cols-' + M.mode;
  function defaultSel() {
    var s = {};
    KEY_COLS.forEach(function (c) { if (hasCol(c)) s[c] = 1; });
    return s;
  }
  function readSel() {
    try { var v = store(colStore); if (v) return JSON.parse(v); } catch (e) { /* fall through */ }
    return defaultSel();
  }
  var sel = readSel();
  function saveSel() { store(colStore, JSON.stringify(sel)); }
  function selectPreset(kind) {
    sel = {};
    if (kind === 'all') T.cols.forEach(function (c) { sel[c] = 1; });
    else sel = defaultSel();
    saveSel();
  }

  // ------------------------------------------------------------------ table view
  var STICKY_TABLES = [];
  function colType(c) {
    var a = col(c);
    for (var i = 0; i < a.length; i++) { if (a[i] !== null && a[i] !== undefined) return typeof a[i]; }
    return 'string';
  }
  function tableView(host, opts) {
    opts = opts || {};
    var colsFn = opts.cols || function () { return T.cols; };
    var types = {}; T.cols.forEach(function (c) { types[c] = colType(c); });
    var pageSize = opts.pageSize || 25;
    var page = 0, filter = '', sortCol = null, sortDir = 'asc', pred = opts.pred || null;
    var expandedIds = {};
    var rowClassFn = opts.rowClass || function () { return ''; };

    var search = h('input', { type: 'search', placeholder: 'Filter rows\u2026', 'aria-label': 'Filter rows' });
    var count = h('span', { class: 'muted small' });
    var pager = h('span', { class: 'pager' });
    var thead = h('thead'), tbody = h('tbody');
    var tableEl = h('table', { class: 'tbl', 'aria-label': opts.label || 'Results' }, thead, tbody);
    var dl = h('button', { class: 'btn', type: 'button', title: 'Download the filtered rows (all columns) as CSV',
      onclick: function () {
        var lines = [T.cols.map(csvField).join(',')];
        idx.forEach(function (i) { lines.push(T.cols.map(function (c) { return csvField(col(c)[i]); }).join(',')); });
        download(opts.filename || 'TALOS_results.csv', lines.join('\n') + '\n');
      } }, iconEl('download', 13), 'CSV');
    var bar = h('div', { class: 'tbl-bar' }, search, count, dl);
    if (opts.toolbar) append(bar, opts.toolbar);
    if (opts.pageSizes) {
      var sizeSel = h('select', { 'aria-label': 'Rows per page', title: 'Rows per page' },
        opts.pageSizes.map(function (n) { return h('option', { value: n === Infinity ? 'all' : String(n), text: n === Infinity ? 'All' : String(n) }); }));
      var saved = store('talos-page-size');
      if (saved === 'all') pageSize = Infinity;
      else if (saved && opts.pageSizes.indexOf(+saved) >= 0) pageSize = +saved;
      sizeSel.value = pageSize === Infinity ? 'all' : String(pageSize);
      sizeSel.addEventListener('change', function () {
        pageSize = sizeSel.value === 'all' ? Infinity : +sizeSel.value;
        store('talos-page-size', sizeSel.value);
        page = 0; recompute();
      });
      bar.appendChild(sizeSel);
    }
    bar.appendChild(pager);
    host.appendChild(bar);
    host.appendChild(h('div', { class: 'tbl-wrap' }, tableEl));

    var view = {
      refresh: function () { draw(); },
      recompute: function () { recompute(); },
      setPredicate: function (p) { pred = p; page = 0; recompute(); },
      setDensity: function (d) {
        tableEl.classList.toggle('d-compact', d === 'compact');
        setTimeout(applySticky, 0);
      },
      applySticky: function () { applySticky(); }
    };
    STICKY_TABLES.push(view);

    function recompute() {
      var f = filter.toLowerCase();
      idx = range(NROW).filter(function (i) {
        if (pred && !pred(i)) return false;
        if (!f) return true;
        return T.cols.some(function (c) { var v = col(c)[i]; return v !== null && v !== undefined && String(v).toLowerCase().indexOf(f) >= 0; });
      });
      if (sortCol) {
        var a = col(sortCol), dir = sortDir === 'asc' ? 1 : -1, num = types[sortCol] === 'number';
        idx.sort(function (p, q) {
          var x = a[p], y = a[q];
          if (x === null || x === undefined) return 1;
          if (y === null || y === undefined) return -1;
          return (num ? x - y : String(x).localeCompare(String(y))) * dir;
        });
      }
      page = Math.min(page, Math.max(0, Math.ceil(idx.length / pageSize) - 1));
      draw();
    }

    function draw() {
      var cs = colsFn(), hasDetail = !!opts.detail;
      var tr = h('tr');
      if (hasDetail) tr.appendChild(h('th', { class: 'exp-h no-sort', scope: 'col', 'aria-label': 'Expand row' }));
      cs.forEach(function (c) {
        tr.appendChild(h('th', {
          scope: 'col', tabindex: '0', text: c, title: HELP[c] || '',
          class: 'g-' + groupOf(c) + (types[c] === 'number' ? ' num' : ''),
          'aria-sort': sortCol === c ? (sortDir === 'asc' ? 'ascending' : 'descending') : 'none',
          onclick: function () { sortDir = sortCol === c && sortDir === 'asc' ? 'desc' : 'asc'; sortCol = c; recompute(); },
          onkeydown: function (e) { if (e.key === 'Enter') this.click(); }
        }));
      });
      clear(thead).appendChild(tr);
      clear(tbody);
      var start = pageSize === Infinity ? 0 : page * pageSize;
      var end = Math.min(idx.length, start + pageSize);
      for (var k = start; k < end; k++) (function (i) {
        var cls = [opts.onRow ? 'clickable' : '', rowClassFn(i)].filter(Boolean).join(' ');
        var row = h('tr', { class: cls, title: opts.onRow ? 'Open full details' : null });
        if (hasDetail) {
          row.appendChild(h('td', { class: 'exp-c' },
            h('button', { class: 'exp', type: 'button', 'aria-expanded': expandedIds[i] ? 'true' : 'false',
              'aria-label': (expandedIds[i] ? 'Collapse' : 'Expand') + ' row details',
              onclick: function (e) { e.stopPropagation(); toggleDetail(i); } },
              iconEl('chevron', 12))));
        }
        cs.forEach(function (c) {
          var v = col(c)[i], r = opts.render && opts.render[c];
          if (r) { row.appendChild(h('td', { class: types[c] === 'number' ? 'num' : '' }, r(v, i))); return; }
          var td = h('td', { class: types[c] === 'number' ? 'num' : '' }, fmt(v));
          if (typeof v === 'string' && v.length > 48) { td.title = v; td.textContent = v.slice(0, 45) + '\u2026'; }
          row.appendChild(td);
        });
        tbody.appendChild(row);
        if (hasDetail && expandedIds[i]) tbody.appendChild(buildDetail(i, cs.length + 1));
      })(idx[k]);
      count.textContent = idx.length + (idx.length === NROW ? '' : ' of ' + NROW) + ' rows';
      drawPager();
      applySticky();
    }

    function toggleDetail(i) {
      if (expandedIds[i]) delete expandedIds[i];
      else expandedIds[i] = 1;
      draw();
    }
    function buildDetail(i, colspan) {
      var body = opts.detail(i);
      return h('tr', { class: 'detail-row' },
        h('td', { colspan: String(colspan) }, body == null ? emptyNote('No details available.') : body));
    }

    function drawPager() {
      var pages = Math.max(1, Math.ceil(idx.length / pageSize));
      clear(pager);
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '\u2039', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
      pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '\u203a', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
    }

    var STICKY_N = opts.sticky == null ? 2 : opts.sticky;
    function applySticky() {
      if (tableEl.offsetParent === null) return;
      var rows = tableEl.rows;
      if (!rows.length) return;
      var ncells = rows[0].cells.length;
      var want = (opts.detail ? 1 : 0) + STICKY_N;
      if (ncells < want) return;
      $$('th.sticky, td.sticky', tableEl).forEach(function (c) {
        c.classList.remove('sticky', 'sticky-end'); c.style.left = '';
      });
      var x = 0, last = want - 1;
      for (var c = 0; c < want; c++) {
        var w = rows[0].cells[c].getBoundingClientRect().width;
        for (var r = 0; r < rows.length; r++) {
          if (rows[r].cells.length !== ncells) continue;
          var cell = rows[r].cells[c];
          cell.classList.add('sticky');
          if (c === last) cell.classList.add('sticky-end');
          cell.style.left = x + 'px';
        }
        x += w;
      }
    }

    var idx = range(NROW);
    var timer = null;
    search.addEventListener('input', function () { clearTimeout(timer); timer = setTimeout(function () { filter = search.value; page = 0; recompute(); }, 120); });
    recompute();
    return view;
  }
  var resizeTimer = null;
  window.addEventListener('resize', function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () { STICKY_TABLES.forEach(function (v) { v.applySticky(); }); }, 150);
  });

  // ------------------------------------------------------------------ tabs
  var TABS = [];
  function tab(id, label, available, build) { TABS.push({ id: id, label: label, available: available, build: build, built: false }); }
  function showVariant(id) {
    activate('details', true);
    var el = document.getElementById('v-' + id);
    if (!el) return;
    $$('.vcard.hl').forEach(function (x) { x.classList.remove('hl'); });
    el.classList.add('hl');
    el.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }
  var resultsAPI = null, pendingQF = null;
  function gotoResults(qf) {
    pendingQF = qf;
    activate('results', true);
    if (resultsAPI && pendingQF) { resultsAPI.setQuick(pendingQF); pendingQF = null; }
  }
  function panelIndexOfVariant(vid) {
    for (var p = 0; p < PANELS.length; p++) {
      var vs = arr(PANELS[p].variants);
      for (var k = 0; k < vs.length; k++) if (vs[k].vid === vid) return p;
    }
    return -1;
  }
  var genomeAPI = null;
  function showInGenome(vid) {
    var p = panelIndexOfVariant(vid);
    if (p < 0) return false;
    activate('genome', true);
    if (genomeAPI) genomeAPI.focusPanel(p, vid);
    return true;
  }

  // ---- Overview
  tab('overview', 'Overview', function () { return true; }, function (root) {
    var samples = uniq(col('Sample')), genes = uniq(col('Gene'));
    var vaf = col('AlleleFrequency').map(function (v) { return typeof v === 'number' ? v * 100 : NaN; });
    var hot = col('Hotspot').filter(function (v) { return v === true; }).length;
    var flagged = col('MergeFragmentWarning').concat(col('ArtifactSuspect'))
                  .filter(function (v) { return v === true; }).length;
    var imputed = col('SequenceSource').filter(function (v) { return v === 'ref_imputed'; }).length;
    var ptd = 0;
    if (!isAlu && hasCol('Length')) {
      ptd = col('Length').filter(function (l) { return typeof l === 'number' && (l === 0 || l > 1000); }).length;
    }

    if (NROW === 0) {
      root.appendChild(h('div', { class: 'banner neutral' },
        h('span', { class: 'btitle' }, iconEl('info', 15), 'No variants detected'),
        h('span', { class: 'small', style: { color: 'var(--text-2)' },
          text: 'No ' + M.events + ' passed the filters in this run. The Results and Genome view tabs are hidden; ' +
                'check the Configuration and Run tabs for the thresholds that were applied.' })));
    } else {
      var banner = h('div', { class: 'banner' + (flagged ? '' : ' neutral') },
        h('span', { class: 'btitle' }, iconEl(flagged ? 'warn' : 'check', 15),
          flagged ? 'Review required' : 'Run summary'));
      var chips = h('span', { class: 'bchips' });
      function bchip(label, qf, title) {
        chips.appendChild(h('button', { class: 'bchip', type: 'button', title: title || '', text: label,
          onclick: function () { gotoResults(qf); } }));
      }
      bchip(NROW + ' ' + M.events + ' called', null, 'Open the results table');
      if (hot) bchip('\u2605 ' + hot + ' in known hotspots', 'hotspot');
      if (flagged) bchip('\u26A0 ' + flagged + ' flagged for review', 'flagged');
      if (imputed) bchip(imputed + ' with reference-imputed sequence', null);
      if (!isAlu && ptd) bchip(ptd + ' PTD-scale calls', 'ptd');
      if (!flagged && !hot) chips.appendChild(h('span', { class: 'small', style: { color: 'var(--text-2)' }, text: 'Nothing flagged for review.' }));
      banner.appendChild(chips);
      root.appendChild(banner);
    }

    var tiles = [[NROW, M.events.charAt(0).toUpperCase() + M.events.slice(1) + ' called'], [samples.length, 'Samples'],
      [genes.length, 'Genes'], [hot, 'In hotspots']];
    if (flagged > 0) tiles.push([flagged, '\u26A0 Flagged calls', 'warn']);
    var mv = median(vaf);
    if (isFinite(mv)) tiles.push([mv.toFixed(1) + '%', 'Median VAF']);
    if (!isAlu && hasCol('Length')) {
      var big = col('Length').filter(function (l) { return typeof l === 'number' && l > 1000; }).length;
      tiles.push([big, 'Calls > 1 kb (PTD scale)']);
    }
    root.appendChild(h('div', { class: 'tiles' }, tiles.map(function (t) {
      return h('div', { class: 'tile' + (t[2] ? ' ' + t[2] : '') }, h('div', { class: 'v', text: String(t[0]) }), h('div', { class: 'l', text: t[1] }));
    })));
    if (!NROW) { root.appendChild(h('div', { style: { marginTop: '12px' } }, card('Results', null, emptyNote('No variants detected.')))); return; }

    var grid = h('div', { class: 'grid g2', style: { marginTop: '12px' } });
    var tb = h('tbody'), groups = {};
    VARIANTS.forEach(function (v) { var k = v.sample + '\u0001' + v.gene; (groups[k] = groups[k] || []).push(v); });
    Object.keys(groups).forEach(function (k) {
      var g = groups[k], vs = g.map(function (v) { var x = col('AlleleFrequency')[v.id - 1]; return typeof x === 'number' ? x * 100 : NaN; });
      tb.appendChild(h('tr', null, h('td', { text: g[0].sample }), h('td', { text: g[0].gene }), h('td', { class: 'num', text: String(g.length) }),
        h('td', { class: 'num', text: isFinite(Math.max.apply(null, vs)) ? Math.max.apply(null, vs).toFixed(1) + '%' : '' })));
    });
    grid.appendChild(card('Calls per sample and gene', null, h('div', { class: 'tbl-wrap', style: { maxHeight: '240px' } }, h('table', { class: 'tbl' },
      h('thead', null, h('tr', null, ['Sample', 'Gene', 'Calls', 'Max VAF'].map(function (x) { return h('th', { text: x }); }))), tb))));
    var box = plotBox();
    grid.appendChild(card('Allele frequency per call', 'click a bar for details', box));
    root.appendChild(grid);
    var labels = VARIANTS.map(function (v) { return v.sample + ' \u00b7 ' + v.gene + ' #' + v.n_in_group; });
    plot(box, function (t) {
      var byGene = uniq(VARIANTS.map(function (v) { return v.gene; }));
      // One trace per gene so the bar colours are explained by a legend.
      var data = byGene.map(function (g, gi) {
        var vs = VARIANTS.filter(function (v) { return v.gene === g; });
        return { type: 'bar', name: g, x: vs.map(function (v) { return v.sample + ' \u00b7 ' + v.gene + ' #' + v.n_in_group; }),
          y: vs.map(function (v) { var x = col('AlleleFrequency')[v.id - 1]; return typeof x === 'number' ? x * 100 : null; }),
          customdata: vs.map(function (v) { return v.id; }), text: vs.map(function (v) { return v.title; }),
          marker: { color: color(gi, t) }, hovertemplate: '%{text}<br>VAF %{y:.1f}%<extra></extra>' };
      });
      return { data: data,
        layout: baseLayout(t, { showlegend: true,
          legend: { orientation: 'h', x: 0, y: 1.12, font: { size: 11 } },
          xaxis: axisStyle(t, { type: 'category', tickangle: -30 }), yaxis: axisStyle(t, { title: { text: 'VAF (%)' }, rangemode: 'tozero' }), margin: { l: 56, r: 16, t: 46, b: 110 } }) };
    });
    onPlotClick(box, function (ev) { if (ev.points && ev.points[0]) showVariant(ev.points[0].customdata); });
  });

  // ---- Shared cell renderers for the results table
  function dash() { return h('span', { class: 'cell-muted', text: '\u2014' }); }
  function tagEl(text, cls, title, ic) {
    return h('span', { class: 'tag ' + (cls || ''), title: title || '' }, ic ? iconEl(ic, 11) : null, text);
  }
  function warnIcon(v, title) {
    return v === true ? h('span', { class: 'cell-warn', title: title || '', text: '\u26A0' }) : dash();
  }
  function pctCell(v, warnBelow) {
    if (typeof v !== 'number' || !isFinite(v)) return dash();
    var bad = warnBelow != null && v < warnBelow;
    return h('span', { class: bad ? 'cell-warn' : '', text: v.toFixed(1) + '%' });
  }
  function monoCell(v) {
    if (v === null || v === undefined || v === '') return dash();
    var s = String(v);
    var el = h('span', { class: 'mono', text: s.length > 16 ? s.slice(0, 14) + '\u2026' : s });
    if (s.length > 16) el.title = s;
    return el;
  }
  function makeResultRenderers() {
    var maxSupport = 0;
    col('SupportingReads').forEach(function (v) { if (typeof v === 'number' && v > maxSupport) maxSupport = v; });
    var r = {
      AlleleFrequency: function (v) {
        if (typeof v !== 'number' || !isFinite(v)) return dash();
        return tagEl(fmtPct(v), v >= 0.10 ? 'ok' : v >= 0.01 ? 'info' : 'mut', 'Variant allele frequency (size-bias corrected)');
      },
      GenomicPosition: function (v) { return h('span', { text: fmtInt(v) }); },
      InsertionSite: function (v) { return h('span', { text: fmtInt(v) }); },
      DepthAtBreakpoint: function (v) { return h('span', { text: fmtInt(v) }); },
      WildtypeReads: function (v) { return h('span', { text: fmtInt(v) }); },
      SupportingReads: function (v) {
        if (typeof v !== 'number') return dash();
        var w = maxSupport > 0 ? Math.max(2, Math.round(v / maxSupport * 40)) : 2;
        return h('span', { class: 'cell-support' },
                 h('span', { class: 'bar', style: { width: w + 'px' } }),
                 h('span', { text: fmtInt(v) }));
      },
      Length: function (v) {
        if (typeof v !== 'number' || !isFinite(v)) return dash();
        if (v === 0) return tagEl('PTD', 'ptd', 'Zero-length call from soft-clips only');
        if (v > 1000) return tagEl(fmtInt(v) + ' bp \u00b7 PTD', 'ptd', 'PTD-scale duplication');
        return h('span', { text: fmtInt(v) + ' bp' });
      },
      Region: function (v) {
        if (!v) return dash();
        return v === 'exonic' ? tagEl('exonic', 'info') : tagEl(String(v), 'mut');
      },
      Orientation: function (v) {
        if (!v) return dash();
        return tagEl(String(v), v === 'tandem' ? 'ok' : v === 'inverted' ? 'danger' : 'mut');
      },
      SequenceSource: function (v) {
        if (!v) return dash();
        return tagEl(String(v), v === 'observed' ? 'ok' : v === 'ref_imputed' ? 'mut' : v === 'merged_intervals' ? 'mut' : 'warn');
      },
      SequenceImputed: function (v) { return v === true ? tagEl('imputed', 'mut') : dash(); },
      Hotspot: function (v, i) {
        return v === true ? h('span', { class: 'cell-hotspot', title: cellRaw('HotspotName', i) || 'Known hotspot', text: '\u2605' }) : dash();
      },
      MergeFragmentWarning: function (v) { return warnIcon(v, 'Merged interval pools more fragments than expected - review before reporting'); },
      ArtifactSuspect: function (v) { return warnIcon(v, 'Fully imputed sequence with ~100% reference match - likely artifact'); },
      StrandBias: function (v) {
        if (typeof v !== 'number' || !isFinite(v)) return dash();
        var bad = v < 0.2 || v > 0.8;
        return h('span', { class: bad ? 'cell-warn' : '', text: v.toFixed(3), title: bad ? 'Strong strand bias' : '' });
      },
      MeanSupportMAPQ: function (v) {
        if (typeof v !== 'number' || !isFinite(v)) return dash();
        return h('span', { class: v < 20 ? 'cell-warn' : '', text: v.toFixed(1), title: v < 20 ? 'Low mean MAPQ' : '' });
      },
      ITDReadCoverage: function (v) { return pctCell(v, 50); },
      RefMatch_Observed: function (v) { return pctCell(v); },
      RefMatch_Total: function (v) { return pctCell(v); },
      ITD_Sequence: monoCell, ITDCoverageRLE: monoCell, TSD_Sequence: monoCell,
      HotspotName: monoCell, CandidateTypes: monoCell, DupExons: monoCell
    };
    return r;
  }

  // ---- Per-call warnings for the inline detail row
  function callWarnings(i) {
    var out = [];
    function add(sev, ic, text) { out.push({ sev: sev, ic: ic, text: text }); }
    if (cellRaw('ArtifactSuspect', i) === true) add('danger', 'warn', 'Artifact suspect: sequence fully imputed from the reference with ~100% identity match.');
    if (cellRaw('MergeFragmentWarning', i) === true) add('warn', 'warn', 'Merged from ' + fmt(cellRaw('MergedFromN', i)) + ' fragments - more than expected for a single true PTD; review before reporting.');
    var src = cellRaw('SequenceSource', i);
    if (src === 'ref_imputed') add('warn', 'info', 'ITD sequence imputed from the reference (not read-derived).');
    else if (src === 'partial_read+ref') add('info', 'info', 'ITD sequence partially imputed from the reference.');
    var sb = numv('StrandBias', i);
    if (sb != null && (sb < 0.2 || sb > 0.8)) add('warn', 'warn', 'Strong strand bias (' + sb.toFixed(3) + ' of supporting reads on one strand).');
    var mq = numv('MeanSupportMAPQ', i);
    if (mq != null && mq < 20) add('warn', 'warn', 'Low mean MAPQ of supporting reads (' + mq.toFixed(1) + ').');
    var cov = numv('ITDReadCoverage', i), len = numv('Length', i);
    if (cov != null && cov < 50 && !(len != null && len > 1000)) add('warn', 'warn', 'Only ' + cov.toFixed(0) + '% of the ITD sequence is covered by reads.');
    if (cellRaw('Region', i) === 'intronic') add('info', 'info', 'Breakpoint is intronic.');
    if (cellRaw('Orientation', i) === 'inverted') add('warn', 'warn', 'Duplicated segment appears inverted.');
    if (cellRaw('Hotspot', i) === true) add('good', 'star', 'Known hotspot' + (cellRaw('HotspotName', i) ? ': ' + cellRaw('HotspotName', i) : '') + '.');
    var af = numv('AlleleFrequency', i);
    if (af != null && af < 0.02) add('info', 'info', 'Low allele frequency (' + fmtPct(af, 2) + ').');
    return out;
  }

  // ---- Inline detail row body
  function resultDetail(i) {
    var grid = h('div', { class: 'detail-grid' });
    var af = numv('AlleleFrequency', i), len = numv('Length', i);
    var stats;
    if (isAlu) {
      stats = [
        ['VAF', af == null ? 'N/A' : fmtPct(af), af != null && af >= 0.1 ? 'good' : ''],
        ['Subtype', fmt(cellRaw('ALUSubtype', i)) || '\u2014', ''],
        ['Orientation', fmt(cellRaw('ALUOrientation', i)) || '\u2014', ''],
        ['TSD', numv('TSD_Length', i) == null ? 'N/A' : fmtInt(numv('TSD_Length', i)) + ' bp', ''],
        ['Poly-A', numv('PolyALength', i) == null ? 'N/A' : fmtInt(numv('PolyALength', i)) + ' nt', ''],
        ['Support', fmtInt(numv('SupportingReads', i)) || 'N/A', ''],
        ['Depth', fmtInt(numv('DepthAtBreakpoint', i)) || 'N/A', '']
      ];
    } else {
      var lenTxt = len == null ? 'N/A' : len === 0 ? 'PTD (clip)' : fmtInt(len) + ' bp';
      var mq = numv('MeanSupportMAPQ', i), sb = numv('StrandBias', i);
      stats = [
        ['VAF', af == null ? 'N/A' : fmtPct(af), af != null && af >= 0.1 ? 'good' : ''],
        ['Length', lenTxt, ''],
        ['Support', fmtInt(numv('SupportingReads', i)) || 'N/A', ''],
        ['Depth', fmtInt(numv('DepthAtBreakpoint', i)) || 'N/A', ''],
        ['MAPQ', mq == null ? 'N/A' : mq.toFixed(1), mq != null && mq < 20 ? 'warn' : ''],
        ['Strand bias', sb == null ? 'N/A' : sb.toFixed(3), sb != null && (sb < 0.2 || sb > 0.8) ? 'warn' : ''],
        ['ITD read cov.', numv('ITDReadCoverage', i) == null ? 'N/A' : numv('ITDReadCoverage', i).toFixed(0) + '%',
          numv('ITDReadCoverage', i) != null && numv('ITDReadCoverage', i) < 50 && !(len != null && len > 1000) ? 'warn' : '']
      ];
    }
    grid.appendChild(h('div', { class: 'mini-stats' }, stats.map(function (s) {
      return h('div', { class: 'mini-stat' },
        h('div', { class: 'v ' + (s[2] || ''), text: String(s[1]) }),
        h('div', { class: 'l', text: s[0] }));
    })));

    var mid = h('div');
    var chips = h('div', { class: 'chips' });
    if (cellRaw('Region', i)) chips.appendChild(tagEl(String(cellRaw('Region', i)), cellRaw('Region', i) === 'exonic' ? 'info' : 'mut', 'Region'));
    if (cellRaw('Orientation', i)) chips.appendChild(tagEl(String(cellRaw('Orientation', i)), cellRaw('Orientation', i) === 'tandem' ? 'ok' : 'danger', 'Orientation'));
    var src = cellRaw('SequenceSource', i);
    if (src) chips.appendChild(tagEl('seq: ' + src, src === 'observed' ? 'ok' : src === 'ref_imputed' ? 'mut' : 'warn', 'Sequence source'));
    if (cellRaw('Hotspot', i) === true) chips.appendChild(tagEl(cellRaw('HotspotName', i) || 'hotspot', 'warn', 'Known hotspot', 'star'));
    var hg = cellRaw('HGVS_cDNA', i);
    if (hg) chips.appendChild(tagEl(String(hg), 'mut', 'HGVS cDNA'));
    var hp = cellRaw('HGVS_Protein', i);
    if (hp) chips.appendChild(tagEl(String(hp), 'mut', 'HGVS protein'));
    if (cellRaw('DupExons', i)) chips.appendChild(tagEl('exons ' + cellRaw('DupExons', i), 'ptd', 'Duplicated exons'));
    if (!chips.childNodes.length) chips.appendChild(h('span', { class: 'small muted', text: 'No status chips for this call.' }));
    mid.appendChild(chips);
    var warns = callWarnings(i);
    if (warns.length) {
      mid.appendChild(h('div', { class: 'warn-list' }, warns.map(function (w) {
        return h('div', { class: 'warn-item sev-' + w.sev, title: '' }, iconEl(w.ic, 13),
          h('span', { text: w.text }));
      })));
    } else {
      mid.appendChild(h('div', { class: 'warn-list' }, h('div', { class: 'warn-item sev-good' }, iconEl('check', 13), h('span', { text: 'No warnings for this call.' }))));
    }
    grid.appendChild(mid);

    var right = h('div', { class: 'stack', style: { gap: '8px' } });
    right.appendChild(h('button', { class: 'btn primary', type: 'button', onclick: function () { showVariant(i + 1); } }, iconEl('eye', 12), 'Full details'));
    if (panelIndexOfVariant(i + 1) >= 0) {
      right.appendChild(h('button', { class: 'btn', type: 'button', onclick: function () { showInGenome(i + 1); } }, iconEl('jump', 12), 'Genome view'));
    }
    var seq = cellRaw(isAlu ? 'ALU_Sequence' : 'ITD_Sequence', i);
    if (seq) right.appendChild(h('details', { class: 'seq' }, h('summary', { text: (isAlu ? 'ALU' : 'ITD') + ' sequence (' + seq.length + ' bp)' }), h('pre', { text: seq })));
    grid.appendChild(right);
    return grid;
  }

  // ---- Results table
  tab('results', 'Results table', function () { return NROW > 0; }, function (root) {
    var host = h('div');
    var render = makeResultRenderers();
    var rowClass = function (i) {
      var cls = [];
      if (hasCol('Hotspot') && cellRaw('Hotspot', i) === true) cls.push('row-hotspot');
      if (cellRaw('MergeFragmentWarning', i) === true || cellRaw('ArtifactSuspect', i) === true) cls.push('row-flagged');
      var len = numv('Length', i);
      if (len != null && (len === 0 || len > 1000)) cls.push('row-ptd');
      return cls.join(' ');
    };
    function visible() { return orderCols(sel); }

    // -- column chooser (Illumina variant-grid style)
    var colBoxes = [];
    var shownLabel = h('span', { class: 'small muted' });
    var menu = h('details', { class: 'colmenu' });
    var sum = h('summary', { class: 'btn', role: 'button', tabindex: '0', title: 'Choose which columns are shown' }, iconEl('columns', 13), 'Columns');
    menu.appendChild(sum);
    var panel = h('div', { class: 'colmenu-panel' });
    panel.appendChild(shownLabel);
    panel.appendChild(h('div', { class: 'chips', style: { margin: '6px 0 4px' } },
      h('button', { class: 'btn', type: 'button', text: 'Clinical view', title: 'Show only the key columns', onclick: function () { selectPreset('key'); syncMenu(); view.refresh(); } }),
      h('button', { class: 'btn', type: 'button', text: 'All columns', title: 'Show every column', onclick: function () { selectPreset('all'); syncMenu(); view.refresh(); } }),
      h('button', { class: 'btn', type: 'button', text: 'Done', title: 'Close the column chooser', onclick: function () { menu.open = false; } })));
    panel.appendChild(h('div', { class: 'hintline', text: 'Changes apply immediately and are remembered for future reports.' }));
    var menuGroups = GROUPS.slice();
    var internalCols = HIDDEN.filter(hasCol);
    if (internalCols.length) menuGroups.push({ id: 'internal', label: 'Internal', cols: internalCols });
    var otherCols = T.cols.filter(function (c) { return groupOf(c) === 'other' && HIDDEN.indexOf(c) < 0; });
    if (otherCols.length) menuGroups.push({ id: 'other', label: 'Other', cols: otherCols });
    menuGroups.forEach(function (g) {
      var fs = h('fieldset');
      // Hidden columns are listed under Internal only, so no column appears twice.
      var present = g.cols.filter(function (c) { return hasCol(c) && (g.id === 'internal' || HIDDEN.indexOf(c) < 0); });
      if (!present.length) return;
      fs.appendChild(h('legend', null, h('span', { class: 'gdot', style: { background: 'var(--gc-' + g.id + ')' } }), g.label + ' (' + present.length + ')'));
      var gridEl = h('div', { class: 'cols' });
      present.forEach(function (c) {
        var cb = h('input', { type: 'checkbox' });
        cb.checked = !!sel[c];
        cb.addEventListener('change', function () { if (cb.checked) sel[c] = 1; else delete sel[c]; saveSel(); syncMenu(); view.refresh(); });
        colBoxes.push({ cb: cb, c: c });
        gridEl.appendChild(h('label', { title: HELP[c] || '' }, cb, c));
      });
      fs.appendChild(gridEl);
      panel.appendChild(fs);
    });
    menu.appendChild(panel);
    function syncMenu() {
      shownLabel.textContent = orderCols(sel).length + ' of ' + T.cols.length + ' columns shown';
      colBoxes.forEach(function (b) { b.cb.checked = !!sel[b.c]; });
    }

    // -- quick filters
    var qmode = null;
    function qpred(i) {
      if (qmode === 'hotspot') return cellRaw('Hotspot', i) === true;
      if (qmode === 'flagged') return cellRaw('MergeFragmentWarning', i) === true || cellRaw('ArtifactSuspect', i) === true;
      if (qmode === 'itd') { var l = numv('Length', i); return l != null && l > 0 && l <= 1000; }
      if (qmode === 'ptd') { var l2 = numv('Length', i); return l2 != null && (l2 === 0 || l2 > 1000); }
      return true;
    }
    var qchips = h('div', { class: 'qchips', role: 'group', 'aria-label': 'Quick filters' });
    var qdefs = [['all', 'All'], ['hotspot', '\u2605 Hotspot'], ['flagged', '\u26A0 Flagged']];
    if (!isAlu) { qdefs.push(['itd', 'ITD (\u22641 kb)']); qdefs.push(['ptd', 'PTD']); }
    var qbtns = [];
    qdefs.forEach(function (d) {
      var b = h('button', { class: 'qchip', type: 'button', 'aria-pressed': d[0] === 'all' ? 'true' : 'false', text: d[1],
        onclick: function () { setQuick(d[0] === 'all' ? null : d[0]); } });
      qbtns.push({ b: b, id: d[0] }); qchips.appendChild(b);
    });
    function setQuick(mode) {
      qmode = mode;
      qbtns.forEach(function (q) { q.b.setAttribute('aria-pressed', (mode || 'all') === q.id ? 'true' : 'false'); });
      view.setPredicate(mode ? qpred : null);
    }

    // -- density
    var density = store('talos-density') || 'comfortable';
    var dseg = h('div', { class: 'seg', role: 'group', 'aria-label': 'Table density' });
    var dBtns = {};
    [['compact', 'Compact'], ['comfortable', 'Comfortable']].forEach(function (d) {
      var b = h('button', { class: 'seg-btn', type: 'button', 'aria-pressed': density === d[0] ? 'true' : 'false', text: d[1],
        onclick: function () { density = d[0]; store('talos-density', density); syncDensity(); } });
      dBtns[d[0]] = b; dseg.appendChild(b);
    });
    function syncDensity() {
      dBtns.compact.setAttribute('aria-pressed', density === 'compact' ? 'true' : 'false');
      dBtns.comfortable.setAttribute('aria-pressed', density === 'comfortable' ? 'true' : 'false');
      view.setDensity(density);
    }

    var view = tableView(host, {
      cols: visible, toolbar: h('span', { style: { display: 'inline-flex', gap: '10px', alignItems: 'center', flexWrap: 'wrap' } }, qchips, menu, dseg),
      pageSizes: [10, 25, 50, 100, Infinity], filename: 'TALOS_results.csv',
      onRow: function (i) { showVariant(i + 1); },
      render: render, rowClass: rowClass, detail: resultDetail, sticky: 2, label: 'Results'
    });
    // initial quick filter / density / column state
    syncMenu();
    syncDensity();
    if (pendingQF) { setQuick(pendingQF); pendingQF = null; }
    resultsAPI = { setQuick: setQuick };

    // group legend above the table
    function legendFor() {
      var vis = visible(), ids = uniq(vis.map(groupOf));
      return h('div', { class: 'glegend' }, ids.map(function (id) {
        return h('span', { class: 'gk' }, h('span', { class: 'gdot', style: { background: 'var(--gc-' + id + ')' } }), groupLabel(id));
      }));
    }
    var legendHost = h('div');
    var refresh0 = view.refresh;
    view.refresh = function () { refresh0(); clear(legendHost).appendChild(legendFor()); };

    var subtitle = 'Rows: \u2605 hotspot \u00b7 \u26A0 flagged \u00b7 purple edge = PTD \u00b7 bar = support \u00b7 hover a header for its meaning \u00b7 click \u25B8 to expand a row';
    root.appendChild(card(
      M.events.charAt(0).toUpperCase() + M.events.slice(1) + ' detected: ' + NROW +
      ' across ' + uniq(col('Sample')).length + ' sample(s)',
      subtitle, h('div', null, legendHost, host)));
    clear(legendHost).appendChild(legendFor());
  });

  // ---- Genome view
  tab('genome', 'Genome view', function () { return PANELS.length > 0; }, function (root) {
    var cur = 0, box = plotBox('tall'), legend = h('div', { class: 'legend' }), hidden = {};
    var focusHost = h('span', { class: 'field' });
    var sel = h('select', { 'aria-label': 'Sample and gene' },
      PANELS.map(function (p, i) { return h('option', { value: i, text: p.sample + ' \u00b7 ' + p.gene + ' (' + p.chrom + ', ' + p.build + ')' }); }));
    sel.addEventListener('change', function () { cur = +sel.value; hidden = {}; draw(); });
    var resetBtn = h('button', { class: 'btn', type: 'button', text: 'Reset zoom',
      onclick: function () {
        if (!hasPlotly()) return;
        var P = PANELS[cur];
        window.Plotly.relayout(box, { 'xaxis.range': [P.from, P.to] });
      } });
    var bar = h('div', { class: 'toolbar' });
    if (PANELS.length > 1) bar.appendChild(h('label', { class: 'field' }, 'Sample \u00b7 gene', sel));
    bar.appendChild(focusHost);
    bar.appendChild(resetBtn);
    bar.appendChild(h('span', { class: 'note',
      text: 'Hover for details; click a call to open it. Exons are drawn to scale, introns compressed.' }));
    root.appendChild(bar);
    root.appendChild(card('Coverage, exons and calls', 'click a call for details', h('div', null, box, legend)));

    function zoomTo(vid) {
      var vs = arr(PANELS[cur].variants);
      for (var k = 0; k < vs.length; k++) {
        if (vs[k].vid !== vid) continue;
        var v = vs[k], pad = Math.max(30, (v.x1 - v.x0) * 0.3);
        if (hasPlotly()) window.Plotly.relayout(box, { 'xaxis.range': [v.x0 - pad, v.x1 + pad] });
        return;
      }
    }
    genomeAPI = {
      focusPanel: function (p, vid) {
        cur = p; sel.value = String(p); hidden = {}; draw();
        if (vid != null) setTimeout(function () { zoomTo(vid); }, 80);
      }
    };

    function draw() {
      var P = PANELS[cur];
      clear(legend); clear(focusHost);
      var vs = arr(P.variants);
      if (vs.length > 1) {
        var focusSel = h('select', { 'aria-label': 'Focus on call' },
          h('option', { value: '', text: 'All ' + vs.length + ' call(s)' }),
          vs.map(function (v, k) { return h('option', { value: String(k), text: v.label }); }));
        focusSel.addEventListener('change', function () {
          if (!hasPlotly()) return;
          var k = +focusSel.value;
          if (focusSel.value === '' || isNaN(k)) {
            window.Plotly.relayout(box, { 'xaxis.range': [P.from, P.to] });
          } else {
            var v = vs[k], pad = Math.max(30, (v.x1 - v.x0) * 0.3);
            window.Plotly.relayout(box, { 'xaxis.range': [v.x0 - pad, v.x1 + pad] });
            showVariant(v.vid);
          }
        });
        focusHost.appendChild(h('span', { class: 'small muted', text: 'Focus: ' }));
        focusHost.appendChild(focusSel);
      }
      vs.forEach(function (v, k) {
        var b = h('button', { type: 'button',
          'aria-pressed': hidden[v.vid] ? 'false' : 'true',
          onclick: function () { hidden[v.vid] = !hidden[v.vid]; draw(); } },
          h('span', { class: 'swatch', style: { background: color(k, resolved()) } }),
          v.label);
        legend.appendChild(b);
      });
      plot(box, function (t) { return genomeFigure(P, t, hidden); });
      onPlotClick(box, function (ev) {
        var p = ev.points && ev.points[0];
        if (p && p.data && p.data.meta && p.data.meta.vid) showVariant(p.data.meta.vid);
      });
    }
    draw();
  });

  function genomeFigure(P, t, hidden) {
    var c = CHROME[t], data = [], shapes = [], ann = [];
    var vs = arr(P.variants), n = Math.max(1, vs.length);
    if (P.cov) data.push({ type: 'scatter', mode: 'lines', x: arr(P.cov.x), y: arr(P.cov.y), customdata: arr(P.cov.pos), fill: 'tozeroy', name: 'Coverage',
      line: { color: CAT[t][0], width: 1.2 }, fillcolor: alpha(CAT[t][0], 0.22), yaxis: 'y', hovertemplate: P.chrom + ':%{customdata}<br>depth %{y}<extra></extra>' });
    arr(P.exons).forEach(function (e) {
      shapes.push({ type: 'rect', xref: 'x', yref: 'y2', x0: e.x0, x1: e.x1, y0: 0.1, y1: 0.9, fillcolor: alpha(EXON[t], 0.75), line: { color: EXON[t], width: 1 } });
      ann.push({ xref: 'x', yref: 'y2', x: (e.x0 + e.x1) / 2, y: 0.5, text: e.label, showarrow: false, font: { size: 10, color: c.text } });
    });
    if (arr(P.exons).length) data.push({ type: 'scatter', mode: 'markers', x: arr(P.exons).map(function (e) { return (e.x0 + e.x1) / 2; }), y: arr(P.exons).map(function () { return 0.5; }),
      text: arr(P.exons).map(function (e) { return '<b>' + e.label + '</b><br>' + P.chrom + ':' + e.start + '-' + e.end; }), marker: { size: 1, opacity: 0 }, yaxis: 'y2', hoverinfo: 'text', name: 'Exons' });
    arr(P.itd_cov).forEach(function (d) {
      var k = indexOfVid(vs, d.vid);
      if (hidden[d.vid]) return;
      data.push({ type: 'scatter', mode: 'lines', x: arr(d.x), y: arr(d.y), customdata: arr(d.pos), fill: 'tozeroy', yaxis: 'y3', name: 'ITD cov',
        line: { color: color(k, t), width: 2 }, fillcolor: alpha(color(k, t), 0.18), hovertemplate: 'ITD pos %{customdata}<br>reads %{y}<extra></extra>' });
    });
    var step = Math.max(0.1, Math.min(0.18, 0.85 / n));
    vs.forEach(function (v, k) {
      if (hidden[v.vid]) return;
      var y = 0.08 + k * step, col_ = color(k, t);
      var point = v.len <= 0 || Math.abs(v.x1 - v.x0) < 0.5;
      data.push(point
        ? { type: 'scatter', mode: 'markers', x: [v.x0], y: [y], yaxis: 'y4', marker: { symbol: 'line-ns-open', size: 18, color: col_, line: { width: 3, color: col_ } }, text: [v.hover], hoverinfo: 'text', meta: { vid: v.vid }, name: v.label }
        : { type: 'scatter', mode: 'lines', x: [v.x0, v.x1], y: [y, y], yaxis: 'y4', line: { color: col_, width: 12 }, text: [v.hover, v.hover], hoverinfo: 'text', meta: { vid: v.vid }, name: v.label });
    });
    var hide = function (extra) { return axisStyle(t, Object.assign({ showticklabels: false, fixedrange: true }, extra)); };
    return {
      data: data,
      layout: baseLayout(t, {
        shapes: shapes, annotations: ann, hovermode: 'closest', margin: { l: 64, r: 16, t: 16, b: 56 },
        xaxis: axisStyle(t, { anchor: 'y4', range: [P.from, P.to], tickvals: arr(P.ticks), ticktext: arr(P.tick_labels), title: { text: P.chrom + ' (compressed, ' + P.build + (P.strand ? ', strand ' + P.strand : '') + ')' } }),
        yaxis: axisStyle(t, { domain: [0.56, 1], title: { text: 'Depth' }, rangemode: 'tozero' }),
        yaxis2: hide({ domain: [0.46, 0.53], range: [0, 1], anchor: 'x', title: { text: 'Exons' } }),
        yaxis3: axisStyle(t, { domain: [0.32, 0.43], anchor: 'x', title: { text: isAlu ? '' : 'ITD cov' }, rangemode: 'tozero', fixedrange: true }),
        yaxis4: hide({ domain: [0, 0.29], anchor: 'x', range: [0, Math.max(1, 0.16 + n * step)], title: { text: isAlu ? 'Insertions' : 'Calls' } })
      })
    };
  }
  function indexOfVid(vs, vid) { for (var i = 0; i < vs.length; i++) if (vs[i].vid === vid) return i; return 0; }

  // ---- Variant details
  tab('details', 'Variant details', function () { return VARIANTS.length > 0; }, function (root) {
    var q = h('input', { type: 'search', placeholder: 'Search variants\u2026', 'aria-label': 'Search variants' });
    root.appendChild(h('div', { class: 'toolbar' }, h('label', { class: 'field' }, 'Search', q)));
    var list = h('div', { class: 'stack' });
    root.appendChild(list);
    var lastGroup = null, cards = [];
    VARIANTS.forEach(function (v) {
      var g = v.sample + ' \u00b7 ' + v.gene;
      if (g !== lastGroup) { list.appendChild(h('h2', { class: 'group-title', text: 'Sample: ' + v.sample + ' \u00b7 ' + v.gene })); lastGroup = g; }
      var body = h('div', null, kv(v.keys, v.values, true));
      if (v.sequence) body.appendChild(h('details', { class: 'seq' }, h('summary', { text: (isAlu ? 'Matched ALU sequence' : 'ITD sequence') + ' (' + v.sequence.length + ' bp)' }), h('pre', { text: v.sequence })));
      var c = card(v.title, v.flag ? 'review flagged' : null, body, { id: 'v-' + v.id, class: 'card vcard' });
      var hdr = $('header', c);
      if (v.hotspot) hdr.appendChild(h('span', { class: 'badge badge-hotspot', text: '\u2605 ' + (v.hotspot_name || 'hotspot') }));
      if (v.flag) hdr.appendChild(h('span', { class: 'badge flag', text: '\u26A0 flagged' }));
      if (panelIndexOfVariant(v.id) >= 0) {
        hdr.appendChild(h('span', { class: 'hint', style: { marginLeft: 'auto', display: 'inline-flex', gap: '6px' } },
          h('button', { class: 'btn', type: 'button', title: 'Show this call in the genome view', onclick: function () { showInGenome(v.id); } }, iconEl('jump', 12), 'Genome view')));
      }
      c._text = (v.title + ' ' + arr(v.values).join(' ')).toLowerCase();
      cards.push(c); list.appendChild(c);
    });
    q.addEventListener('input', function () { var f = q.value.toLowerCase(); cards.forEach(function (c) { c.hidden = f && c._text.indexOf(f) < 0; }); });
  });

  // ---- Configuration
  tab('config', 'Configuration', function () { return arr(DATA.configs).length > 0; }, function (root) {
    root.appendChild(h('div', { class: 'grid g2' }, arr(DATA.configs).map(function (c) { return card('Analysis configuration: ' + c.gene, null, kv(c.keys, c.values)); })));
  });

  // ---- Run
  tab('run', 'Run', function () { return true; }, function (root) {
    var R = DATA.run;
    root.appendChild(card('Run information', M.plotly_inline ? 'plotly.js embedded (works offline)' : 'plotly.js loaded from CDN', kv(R.keys, R.values)));
  });

  // ------------------------------------------------------------------ shell
  var panels = {}, tabButtons = {}, active = null;
  function activate(id, push) {
    var t = TABS.filter(function (x) { return x.id === id && x.available(); })[0] || TABS.filter(function (x) { return x.available(); })[0];
    if (!t) return;
    id = t.id;
    Object.keys(panels).forEach(function (k) {
      var on = k === id; panels[k].hidden = !on;
      tabButtons[k].setAttribute('aria-selected', on ? 'true' : 'false'); tabButtons[k].tabIndex = on ? 0 : -1;
    });
    if (!t.built) {
      t.built = true;
      try { t.build(panels[id]); } catch (e) { panels[id].appendChild(emptyNote('This tab could not be built: ' + e.message)); if (window.console) console.error(e); }
    } else refreshVisible(panels[id]);
    setTimeout(function () { refreshVisible(panels[id]); }, 60);
    active = id;
    if (push) { try { history.replaceState(null, '', '#tab=' + id); } catch (e) { /* file:// may refuse */ } }
  }
  function init() {
    applyThemeAttr();
    var tabsEl = $('#tabs'), main = $('#panels');
    var avail = TABS.filter(function (t) { return t.available(); });
    avail.forEach(function (t, i) {
      var b = h('button', { class: 'tab', role: 'tab', id: 'tab-' + t.id, 'aria-controls': 'panel-' + t.id, 'aria-selected': 'false', tabindex: '-1', type: 'button', text: t.label,
        onclick: function () { activate(t.id, true); },
        onkeydown: function (e) {
          var k = e.key === 'ArrowRight' ? 1 : e.key === 'ArrowLeft' ? -1 : 0;
          if (!k) return; e.preventDefault();
          var nx = avail[(i + k + avail.length) % avail.length]; activate(nx.id, true); tabButtons[nx.id].focus();
        } });
      tabButtons[t.id] = b; tabsEl.appendChild(b);
      panels[t.id] = h('section', { class: 'panel', role: 'tabpanel', id: 'panel-' + t.id, 'aria-labelledby': 'tab-' + t.id, hidden: true });
      main.appendChild(panels[t.id]);
    });
    var btn = $('#theme-btn');
    var label = function () { btn.textContent = 'Theme: ' + themePref; };
    label();
    btn.addEventListener('click', function () {
      themePref = themePref === 'auto' ? 'light' : themePref === 'light' ? 'dark' : 'auto';
      store('talos-theme', themePref); applyThemeAttr(); label(); rerenderAll();
    });
    var printBtn = $('#print-btn');
    if (printBtn) printBtn.addEventListener('click', function () { window.print(); });
    if (window.matchMedia) {
      var mq = window.matchMedia('(prefers-color-scheme: dark)');
      if (mq.addEventListener) mq.addEventListener('change', function () { if (themePref === 'auto') rerenderAll(); });
    }
    // Print: build every tab first so the printed page shows all sections.
    window.addEventListener('beforeprint', function () {
      TABS.forEach(function (t) {
        if (t.available() && !t.built) {
          t.built = true;
          try { t.build(panels[t.id]); } catch (e) { if (window.console) console.error(e); }
        }
      });
    });
    var m = /tab=([\w-]+)/.exec(location.hash);
    activate(m ? m[1] : 'overview', false);
    window.addEventListener('hashchange', function () { var mm = /tab=([\w-]+)/.exec(location.hash); if (mm && mm[1] !== active) activate(mm[1], false); });
    document.documentElement.setAttribute('data-ready', '1');
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
)-----"

# ============================================================================
# Report assets embedded in R: stylesheet and client-side app of the
# self-contained HTML report (no inst/ files, no R Markdown). Both strings are
# written verbatim into every report. Raw strings need R >= 4.0; the text is
# ASCII-only (non-ASCII characters are escaped) so the file is locale-safe.
# ============================================================================

.TALOS_REPORT_CSS <- r"-----(
/* TALOS report -- design tokens, layout and components (SPARKS style).
   Colours follow the validated palette: categorical hues are applied by JS in
   fixed order; chrome (surface / ink / grid) swaps between light and dark. */
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
  --accent: #2a78d6;
  --accent-ink: #ffffff;
  --accent-soft: #e7f0fc;
  --good: #006300;
  --warn: #8a5a00;
  --radius: 10px;
  --font: system-ui, -apple-system, "Segoe UI", sans-serif;
  --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) {
  :root:where(:not([data-theme="light"])) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
    --muted: #898781; --grid: #2c2c2a; --axis: #383835;
    --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
    --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --text: #f0efec; --text-2: #c3c2b7;
  --muted: #898781; --grid: #2c2c2a; --axis: #383835;
  --ring: rgba(255, 255, 255, 0.10); --accent: #3987e5; --accent-ink: #0b0b0b;
  --accent-soft: #1c2a3d; --good: #0ca30c; --warn: #fab219;
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

/* ---- Header & tabs ---- */
.top { background: var(--surface); border-bottom: 1px solid var(--ring); }
.top-inner { max-width: 1440px; margin: 0 auto; padding: 14px 16px 0; }
.title-row { display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }
.title-row h1 { font-size: 20px; }
.badge { display: inline-flex; align-items: center; gap: 6px; padding: 2px 10px;
         border-radius: 999px; background: var(--accent-soft); color: var(--text-2);
         font-size: 12px; border: 1px solid var(--ring); }
.badge.demo { background: transparent; border-color: var(--warn); color: var(--warn); font-weight: 600; }
.spacer { flex: 1; }
.tabs { display: flex; gap: 2px; overflow-x: auto; scrollbar-width: thin; margin-top: 10px; }
.tab { appearance: none; border: 0; background: none; color: var(--text-2); font: inherit;
       padding: 10px 14px; cursor: pointer; border-bottom: 2px solid transparent;
       white-space: nowrap; border-radius: 6px 6px 0 0; }
.tab:hover { background: var(--accent-soft); color: var(--text); }
.tab[aria-selected="true"] { color: var(--text); border-bottom-color: var(--accent); font-weight: 600; }
.tab:focus-visible, button:focus-visible, select:focus-visible, input:focus-visible {
  outline: 2px solid var(--accent); outline-offset: 2px; }
main { max-width: 1440px; margin: 0 auto; padding: 16px; }
.panel[hidden] { display: none; }

/* ---- Cards, grids, controls ---- */
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius);
        padding: 14px; min-width: 0; }
.card > header { display: flex; align-items: baseline; gap: 8px; margin-bottom: 8px; flex-wrap: wrap; }
.card > header .hint { margin-left: auto; }
.grid { display: grid; gap: 12px; }
.g2 { grid-template-columns: repeat(2, minmax(0, 1fr)); }
.g3 { grid-template-columns: repeat(3, minmax(0, 1fr)); }
.g-main { grid-template-columns: minmax(0, 1fr) 360px; align-items: start; }
.stack { display: flex; flex-direction: column; gap: 12px; }
.toolbar { display: flex; gap: 12px; flex-wrap: wrap; align-items: flex-end; margin-bottom: 12px; }
.field { display: flex; flex-direction: column; gap: 3px; font-size: 12px; color: var(--text-2); min-width: 0; max-width: 100%; }
.field.inline { flex-direction: row; align-items: center; gap: 6px; }
select, input[type="text"], input[type="number"], input[type="search"] {
  font: inherit; color: var(--text); background: var(--surface); border: 1px solid var(--axis);
  border-radius: 6px; padding: 5px 8px; min-width: 120px; max-width: 280px; }
input[type="range"] { accent-color: var(--accent); }
button.btn { font: inherit; padding: 5px 12px; border-radius: 6px; cursor: pointer;
             border: 1px solid var(--axis); background: var(--surface); color: var(--text); }
button.btn:hover { background: var(--accent-soft); }
button.btn.primary { background: var(--accent); border-color: var(--accent); color: var(--accent-ink); }
button.btn[disabled] { opacity: 0.45; cursor: not-allowed; }
.plot { width: 100%; min-height: 360px; }
.plot.tall { min-height: 560px; } .plot.short { min-height: 280px; }
.empty { padding: 28px 12px; text-align: center; color: var(--muted); }

/* ---- Stat tiles ---- */
.tiles { display: grid; gap: 12px; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); }
.tile { background: var(--surface); border: 1px solid var(--ring); border-radius: var(--radius); padding: 12px 14px; }
.tile .v { font-size: 26px; font-weight: 600; letter-spacing: -0.02em; }
.tile .l { color: var(--text-2); font-size: 12px; }

/* ---- Legend (HTML so categories can be toggled) ---- */
.legend { display: flex; flex-wrap: wrap; gap: 4px 12px; margin-top: 8px; font-size: 12px; color: var(--text-2); }
.legend button { appearance: none; border: 0; background: none; font: inherit; color: inherit;
                 display: inline-flex; gap: 6px; align-items: center; cursor: pointer; padding: 2px 4px; border-radius: 4px; }
.legend button:hover { background: var(--accent-soft); }
.legend button[aria-pressed="false"] { opacity: 0.4; text-decoration: line-through; }
.swatch { width: 10px; height: 10px; border-radius: 3px; display: inline-block; }

/* ---- Tables ---- */
.tbl-wrap { overflow: auto; border: 1px solid var(--ring); border-radius: 8px; max-height: 520px; }
table.tbl { border-collapse: collapse; width: 100%; font-size: 12.5px; }
.tbl th, .tbl td { padding: 5px 10px; text-align: left; white-space: nowrap; border-bottom: 1px solid var(--grid); }
.tbl th { position: sticky; top: 0; background: var(--surface); cursor: pointer; user-select: none;
          font-weight: 600; color: var(--text-2); border-bottom: 1px solid var(--axis); z-index: 1; }
.tbl th[aria-sort="ascending"]::after { content: " \25B2"; font-size: 9px; }
.tbl th[aria-sort="descending"]::after { content: " \25BC"; font-size: 9px; }
.tbl td.num { text-align: right; font-variant-numeric: tabular-nums; }
.tbl tbody tr:hover { background: var(--accent-soft); }
.tbl-bar { display: flex; gap: 10px; align-items: center; margin: 8px 0; flex-wrap: wrap; }
.pager { display: inline-flex; gap: 6px; align-items: center; margin-left: auto; }

/* ---- Selection panel ---- */
.bars { display: flex; flex-direction: column; gap: 5px; }
.bar-row { display: grid; grid-template-columns: 110px 1fr 70px; gap: 8px; align-items: center; font-size: 12px; }
.bar-row .lab { overflow: hidden; text-overflow: ellipsis; white-space: nowrap; color: var(--text-2); }
.bar-track { background: var(--grid); border-radius: 4px; height: 10px; overflow: hidden; }
.bar-fill { height: 100%; border-radius: 4px 4px 4px 4px; }
.note { font-size: 12px; color: var(--muted); }
.kv { display: grid; grid-template-columns: max-content 1fr; gap: 3px 14px; font-size: 12.5px; }
.kv dt { color: var(--text-2); } .kv dd { margin: 0; overflow-wrap: anywhere; }

/* ---- Gallery ---- */
.gallery { display: grid; gap: 12px; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); }
.shot { background: var(--surface); border: 1px solid var(--ring); border-radius: 8px; overflow: hidden; margin: 0; }
.shot button { all: unset; display: block; cursor: zoom-in; width: 100%; }
.shot img { display: block; width: 100%; height: 190px; object-fit: contain; background: #fff; }
.shot figcaption { padding: 6px 10px; font-size: 12px; color: var(--text-2); overflow-wrap: anywhere; }
.shot .nopreview { height: 190px; display: grid; place-items: center; color: var(--muted); background: var(--page); }
dialog.lightbox { border: 0; padding: 0; background: transparent; max-width: 96vw; max-height: 96vh; }
dialog.lightbox::backdrop { background: rgba(0, 0, 0, 0.7); }
dialog.lightbox img { max-width: 96vw; max-height: 92vh; background: #fff; border-radius: 8px; display: block; }

@media (max-width: 1000px) { .g-main, .g2, .g3 { grid-template-columns: minmax(0, 1fr); } }
@media (max-width: 720px) { main { padding: 12px; } .bar-row { grid-template-columns: 80px 1fr 56px; } }
@media print {
  .tabs, .toolbar, .btn { display: none !important; }
  .panel[hidden] { display: block !important; }
  .card { break-inside: avoid; }
}
@media (forced-colors: active) { .bar-fill, .swatch { forced-color-adjust: none; } }

/* ---- TALOS additions ---- */
.badge.flag { background: transparent; border-color: var(--warn); color: var(--warn); font-weight: 600; }
.vcard { scroll-margin-top: 12px; }
.vcard.hl { outline: 2px solid var(--accent); outline-offset: 2px; }
.vcard h3 { font-size: 14px; }
.group-title { margin: 18px 0 8px; font-size: 15px; color: var(--text-2); }
details.seq summary { cursor: pointer; color: var(--accent); font-size: 12.5px; margin-top: 8px; }
details.seq pre { white-space: pre-wrap; word-break: break-all; margin-top: 6px; }
.tbl tbody tr.clickable { cursor: pointer; }
.chips { display: flex; flex-wrap: wrap; gap: 6px; }
.chip { font-size: 12px; padding: 2px 8px; border-radius: 999px; border: 1px solid var(--ring); background: var(--page); color: var(--text-2); }
.kv.wide { grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 4px 18px; }
.kv.wide > div { display: grid; grid-template-columns: 48% 52%; gap: 8px; border-bottom: 1px dashed var(--grid); padding: 2px 0; }
.kv.wide dt, .kv.wide dd { font-size: 12.5px; }

/* ---- Highlight important results (TALOS) ---- */
.tbl tbody tr.row-hotspot             { background: rgba(250,178,25,0.08); }
:root[data-theme="dark"] .tbl tbody tr.row-hotspot { background: rgba(250,178,25,0.12); }
.tbl tbody tr.row-hotspot td:first-child { box-shadow: inset 3px 0 0 var(--warn); }
.tbl tbody tr.row-flagged             { background: rgba(138,90,0,0.05); }
.tbl tbody tr.row-flagged td:first-child { box-shadow: inset 3px 0 0 var(--warn); }
.tbl tbody tr.row-hotspot:hover,
.tbl tbody tr.row-flagged:hover       { background: var(--accent-soft); }

.cell-flag, .cell-hotspot             { color: var(--warn); font-weight: 600; }
.cell-muted                            { color: var(--muted); }
.cell-hi                               { font-weight: 600; color: var(--good); }
.cell-lo                               { color: var(--muted); }
.cell-ptd                              { color: var(--accent); font-style: italic; font-weight: 600; }
.cell-support                          { display: inline-flex; align-items: center; gap: 6px; }
.cell-support .bar                     { display: inline-block; height: 8px; background: var(--accent);
                                         border-radius: 2px; min-width: 2px; }
.badge-hotspot { background: rgba(250,178,25,0.15); color: var(--warn); border: 1px solid var(--warn);
                 padding: 1px 6px; border-radius: 4px; font-size: 11px; font-weight: 600; }
.badge-exon    { background: var(--accent-soft); color: var(--text); border: 1px solid var(--ring);
                 padding: 1px 6px; border-radius: 4px; font-size: 11px; }
)-----"

.TALOS_REPORT_JS <- r"-----(
/* TALOS report -- client-side app (SPARKS style).
   Reads the JSON blob in #talos-data, builds one lazily-rendered panel per tab
   and draws figures with plotly.js. Tables and details work without plotly. */
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
  var T = DATA.table, NROW = T.n;
  function col(name) { return T.data[name] ? arr(T.data[name]) : []; }

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
  // Plotly attaches .on() at first render, which may come later (CDN fallback).
  function onPlotClick(div, fn) {
    if (div._clickBound) return;
    if (div.on) { div._clickBound = true; div.on('plotly_click', fn); return; }
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

  // ------------------------------------------------------------------ table view
  function colType(c) {
    var a = col(c);
    for (var i = 0; i < a.length; i++) { if (a[i] !== null && a[i] !== undefined) return typeof a[i]; }
    return 'string';
  }
  function tableView(host, cols, opts) {
    opts = opts || {};
    var pageSize = opts.pageSize || 25, page = 0, filter = '', sortCol = null, sortDir = 'asc', idx = range(NROW);
    var types = {}; T.cols.forEach(function (c) { types[c] = colType(c); });
    var search = h('input', { type: 'search', placeholder: 'Filter rows\u2026', 'aria-label': 'Filter rows' });
    var count = h('span', { class: 'muted small' });
    var pager = h('span', { class: 'pager' });
    var thead = h('thead'), tbody = h('tbody');
    var rowClassFn = opts.rowClass || function () { return ''; };
    var dl = h('button', { class: 'btn', type: 'button', text: 'Download CSV', onclick: function () {
      var lines = [T.cols.map(csvField).join(',')];
      idx.forEach(function (i) { lines.push(T.cols.map(function (c) { return csvField(col(c)[i]); }).join(',')); });
      download(opts.filename || 'TALOS_results.csv', lines.join('\n') + '\n');
    } });
    var bar = h('div', { class: 'tbl-bar' }, search, count, dl);
    if (opts.toolbar) append(bar, opts.toolbar);
    bar.appendChild(pager);
    host.appendChild(bar);
    host.appendChild(h('div', { class: 'tbl-wrap' }, h('table', { class: 'tbl' }, thead, tbody)));
    function current() { return typeof cols === 'function' ? cols() : cols; }
    function recompute() {
      var f = filter.toLowerCase();
      idx = range(NROW).filter(function (i) {
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
      var cs = current(), tr = h('tr');
      cs.forEach(function (c) {
        tr.appendChild(h('th', { scope: 'col', tabindex: '0', text: c,
          'aria-sort': sortCol === c ? (sortDir === 'asc' ? 'ascending' : 'descending') : 'none',
          onclick: function () { sortDir = sortCol === c && sortDir === 'asc' ? 'desc' : 'asc'; sortCol = c; recompute(); },
          onkeydown: function (e) { if (e.key === 'Enter') this.click(); } }));
      });
      clear(thead).appendChild(tr);
      clear(tbody);
      var start = page * pageSize, end = Math.min(idx.length, start + pageSize);
      for (var k = start; k < end; k++) {
        (function (i) {
          var rc = rowClassFn(i);
          var cls = [opts.onRow ? 'clickable' : '', rc].filter(Boolean).join(' ');
          var row = h('tr', { class: cls, title: opts.onRow ? 'Show details' : null,
            onclick: opts.onRow ? function () { opts.onRow(i); } : null });
          cs.forEach(function (c) {
            var v = col(c)[i], r = opts.render && opts.render[c];
            row.appendChild(h('td', { class: types[c] === 'number' ? 'num' : '' }, r ? r(v, i) : fmt(v)));
          });
          tbody.appendChild(row);
        })(idx[k]);
      }
      count.textContent = idx.length + (idx.length === NROW ? '' : ' of ' + NROW) + ' rows';
      var pages = Math.max(1, Math.ceil(idx.length / pageSize));
      clear(pager);
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '\u2039', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
      pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '\u203a', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
    }
    var timer = null;
    search.addEventListener('input', function () { clearTimeout(timer); timer = setTimeout(function () { filter = search.value; page = 0; recompute(); }, 120); });
    recompute();
    return { refresh: draw };
  }

  // ------------------------------------------------------------------ tabs
  var TABS = [];
  function tab(id, label, available, build) { TABS.push({ id: id, label: label, available: available, build: build, built: false }); }
  var M = DATA.meta, VARIANTS = arr(DATA.variants), PANELS = arr(DATA.panels), HIDDEN = arr(DATA.hidden);
  var isAlu = M.mode === 'alu';
  function variantIndex(id) { for (var i = 0; i < VARIANTS.length; i++) if (VARIANTS[i].id === id) return i; return -1; }
  function showVariant(id) {
    activate('details', true);
    var el = document.getElementById('v-' + id);
    if (!el) return;
    $$('.vcard.hl').forEach(function (x) { x.classList.remove('hl'); });
    el.classList.add('hl');
    el.scrollIntoView({ behavior: 'smooth', block: 'start' });
  }

  // ---- Overview
  tab('overview', 'Overview', function () { return true; }, function (root) {
    var samples = uniq(col('Sample')), genes = uniq(col('Gene'));
    var vaf = col('AlleleFrequency').map(function (v) { return typeof v === 'number' ? v * 100 : NaN; });
    var hot = col('Hotspot').filter(function (v) { return v === true; }).length;
    var flagged = col('MergeFragmentWarning').concat(col('ArtifactSuspect'))
                  .filter(function (v) { return v === true; }).length;
    var tiles = [[NROW, M.events.charAt(0).toUpperCase() + M.events.slice(1) + ' called'], [samples.length, 'Samples'],
      [genes.length, 'Genes'], [hot, 'In hotspots']];
    if (flagged > 0) tiles.push([flagged, '\u26A0 Flagged calls']);
    var mv = median(vaf);
    if (isFinite(mv)) tiles.push([mv.toFixed(1) + '%', 'Median VAF']);
    if (!isAlu) {
      var ptd = col('Length').filter(function (l) { return typeof l === 'number' && l > 1000; }).length;
      tiles.push([ptd, 'Calls > 1 kb (PTD scale)']);
    }
    root.appendChild(h('div', { class: 'tiles' }, tiles.map(function (t) { return h('div', { class: 'tile' }, h('div', { class: 'v', text: String(t[0]) }), h('div', { class: 'l', text: t[1] })); })));
    if (!NROW) { root.appendChild(h('div', { style: { marginTop: '12px' } }, card('Results', null, emptyNote('No variants detected.')))); return; }

    var grid = h('div', { class: 'grid g2', style: { marginTop: '12px' } });
    // calls per sample / gene
    var tb = h('tbody'), groups = {};
    VARIANTS.forEach(function (v) { var k = v.sample + '\u0001' + v.gene; (groups[k] = groups[k] || []).push(v); });
    Object.keys(groups).forEach(function (k) {
      var g = groups[k], vs = g.map(function (v) { var x = col('AlleleFrequency')[v.id - 1]; return typeof x === 'number' ? x * 100 : NaN; });
      tb.appendChild(h('tr', null, h('td', { text: g[0].sample }), h('td', { text: g[0].gene }), h('td', { class: 'num', text: String(g.length) }),
        h('td', { class: 'num', text: isFinite(Math.max.apply(null, vs)) ? Math.max.apply(null, vs).toFixed(1) + '%' : '' })));
    });
    grid.appendChild(card('Calls per sample and gene', null, h('div', { class: 'tbl-wrap' }, h('table', { class: 'tbl' },
      h('thead', null, h('tr', null, ['Sample', 'Gene', 'Calls', 'Max VAF'].map(function (x) { return h('th', { text: x }); }))), tb))));
    // VAF per call
    var box = plotBox();
    grid.appendChild(card('Allele frequency per call', 'click a bar for details', box));
    root.appendChild(grid);
    var labels = VARIANTS.map(function (v) { return v.sample + ' \u00b7 ' + v.gene + ' #' + v.n_in_group; });
    plot(box, function (t) {
      var byGene = uniq(VARIANTS.map(function (v) { return v.gene; }));
      return { data: [{ type: 'bar', x: labels, y: VARIANTS.map(function (v) { var x = col('AlleleFrequency')[v.id - 1]; return typeof x === 'number' ? x * 100 : null; }),
        customdata: VARIANTS.map(function (v) { return v.id; }), text: VARIANTS.map(function (v) { return v.title; }),
        marker: { color: VARIANTS.map(function (v) { return color(byGene.indexOf(v.gene), t); }) },
        hovertemplate: '%{text}<br>VAF %{y:.1f}%<extra></extra>' }],
        layout: baseLayout(t, { xaxis: axisStyle(t, { type: 'category', tickangle: -30 }), yaxis: axisStyle(t, { title: { text: 'VAF (%)' }, rangemode: 'tozero' }), margin: { l: 56, r: 16, t: 16, b: 110 } }) };
    });
    onPlotClick(box, function (ev) { if (ev.points && ev.points[0]) showVariant(ev.points[0].customdata); });
  });

  // ---- Shared cell renderers for the results table
  function makeResultRenderers() {
    var maxSupport = 0;
    col('SupportingReads').forEach(function (v) { if (typeof v === 'number' && v > maxSupport) maxSupport = v; });
    return {
      AlleleFrequency: function (v) {
        if (typeof v !== 'number') return h('span', { class: 'cell-muted', text: 'N/A' });
        var pct = (v * 100).toFixed(1) + '%';
        var cls = v >= 0.1 ? 'cell-hi' : v < 0.01 ? 'cell-lo' : '';
        return h('span', { class: cls, text: pct });
      },
      Length: function (v) {
        if (typeof v !== 'number') return '';
        if (v === 0)  return h('span', { class: 'cell-ptd', text: 'PTD (clip)' });
        if (v > 1000) return h('span', { class: 'cell-ptd', text: fmt(v) + ' bp (PTD)' });
        return h('span', { text: fmt(v) + ' bp' });
      },
      Region: function (v) {
        if (!v) return h('span', { class: 'cell-muted', text: '\u2014' });
        return h('span', { class: v === 'exonic' ? 'badge-exon' : 'cell-muted', text: v });
      },
      SupportingReads: function (v) {
        if (typeof v !== 'number') return '';
        var w = maxSupport > 0 ? Math.max(2, Math.round(v / maxSupport * 40)) : 2;
        return h('span', { class: 'cell-support' },
                 h('span', { class: 'bar', style: { width: w + 'px' } }),
                 h('span', { text: fmt(v) }));
      },
      Hotspot:              function (v) { return v === true ? h('span', { class: 'cell-hotspot', text: '\u2605' }) : h('span', { class: 'cell-muted', text: '\u2014' }); },
      MergeFragmentWarning: function (v) { return v === true ? h('span', { class: 'cell-flag', text: '\u26A0' }) : ''; },
      ArtifactSuspect:      function (v) { return v === true ? h('span', { class: 'cell-flag', text: '\u26A0' }) : ''; },
      SequenceImputed:      function (v) { return v === true ? h('span', { class: 'cell-muted', text: 'imputed' }) : ''; }
    };
  }

  // ---- Results table
  tab('results', 'Results table', function () { return NROW > 0; }, function (root) {
    var showAll = false, host = h('div');
    var visible = function () { return showAll ? T.cols : T.cols.filter(function (c) { return HIDDEN.indexOf(c) < 0; }); };
    var toggle = h('label', { class: 'field inline small' },
      h('input', { type: 'checkbox', onchange: function (e) { showAll = e.target.checked; view.refresh(); } }),
      'All columns');
    var hasCol = function (c) { return T.cols.indexOf(c) >= 0; };
    var render = makeResultRenderers();
    var rowClass = function (i) {
      var cls = [];
      if (hasCol('Hotspot') && col('Hotspot')[i] === true) cls.push('row-hotspot');
      if ((hasCol('MergeFragmentWarning') && col('MergeFragmentWarning')[i] === true) ||
          (hasCol('ArtifactSuspect')      && col('ArtifactSuspect')[i]      === true)) cls.push('row-flagged');
      if (hasCol('Length')) {
        var len = col('Length')[i];
        if (typeof len === 'number' && (len === 0 || len > 1000)) cls.push('row-ptd');
      }
      return cls.join(' ');
    };
    var view = tableView(host, visible, {
      toolbar: toggle, filename: 'TALOS_results.csv',
      onRow: function (i) { showVariant(i + 1); },
      render: render, rowClass: rowClass
    });
    var subtitle = 'Rows are colour-coded: \u2605 hotspot  \u26A0 flagged  italic = PTD  \u25CF bar = support';
    root.appendChild(card(
      M.events.charAt(0).toUpperCase() + M.events.slice(1) + ' detected: ' + NROW +
      ' across ' + uniq(col('Sample')).length + ' sample(s)',
      subtitle, host));
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
      if (v.flag)    hdr.appendChild(h('span', { class: 'badge flag', text: '\u26A0 flagged' }));
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
    if (window.matchMedia) {
      var mq = window.matchMedia('(prefers-color-scheme: dark)');
      if (mq.addEventListener) mq.addEventListener('change', function () { if (themePref === 'auto') rerenderAll(); });
    }
    var m = /tab=([\w-]+)/.exec(location.hash);
    activate(m ? m[1] : 'overview', false);
    window.addEventListener('hashchange', function () { var mm = /tab=([\w-]+)/.exec(location.hash); if (mm && mm[1] !== active) activate(mm[1], false); });
    document.documentElement.setAttribute('data-ready', '1');
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
)-----"
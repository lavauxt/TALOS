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
  function plot(div, build) {
    if (!hasPlotly()) { clear(div).appendChild(emptyNote('Interactive plots need plotly.js (offline and the CDN was unreachable).')); return div; }
    div._render = function () { var s = build(resolved()); window.Plotly.react(div, s.data, s.layout, Object.assign({}, CONFIG, s.config || {})); div._dirty = false; };
    if (registry.indexOf(div) < 0) registry.push(div);
    div._render();
    return div;
  }
  function rerenderAll() {
    registry = registry.filter(function (d) { return document.body.contains(d) && d._render; });
    registry.forEach(function (d) { if (d.offsetParent !== null) d._render(); else d._dirty = true; });
  }
  function refreshVisible(root) {
    $$('.plot', root).forEach(function (d) {
      if (!d._render) return;
      if (d._dirty) d._render(); else if (hasPlotly() && d.data) window.Plotly.Plots.resize(d);
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
    var search = h('input', { type: 'search', placeholder: 'Filter rows…', 'aria-label': 'Filter rows' });
    var count = h('span', { class: 'muted small' });
    var pager = h('span', { class: 'pager' });
    var thead = h('thead'), tbody = h('tbody');
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
          var row = h('tr', { class: opts.onRow ? 'clickable' : '', title: opts.onRow ? 'Show details' : null,
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
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '‹', 'aria-label': 'Previous page', disabled: page === 0, onclick: function () { page--; draw(); } }));
      pager.appendChild(h('span', { class: 'small', text: (page + 1) + ' / ' + pages }));
      pager.appendChild(h('button', { class: 'btn', type: 'button', text: '›', 'aria-label': 'Next page', disabled: page >= pages - 1, onclick: function () { page++; draw(); } }));
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
    var tiles = [[NROW, M.events.charAt(0).toUpperCase() + M.events.slice(1) + ' called'], [samples.length, 'Samples'],
      [genes.length, 'Genes'], [hot, 'In hotspots']];
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
    var labels = VARIANTS.map(function (v) { return v.sample + ' · ' + v.gene + ' #' + v.n_in_group; });
    plot(box, function (t) {
      var byGene = uniq(VARIANTS.map(function (v) { return v.gene; }));
      return { data: [{ type: 'bar', x: labels, y: VARIANTS.map(function (v) { var x = col('AlleleFrequency')[v.id - 1]; return typeof x === 'number' ? x * 100 : null; }),
        customdata: VARIANTS.map(function (v) { return v.id; }), text: VARIANTS.map(function (v) { return v.title; }),
        marker: { color: VARIANTS.map(function (v) { return color(byGene.indexOf(v.gene), t); }) },
        hovertemplate: '%{text}<br>VAF %{y:.1f}%<extra></extra>' }],
        layout: baseLayout(t, { xaxis: axisStyle(t, { type: 'category', tickangle: -30 }), yaxis: axisStyle(t, { title: { text: 'VAF (%)' }, rangemode: 'tozero' }), margin: { l: 56, r: 16, t: 16, b: 110 } }) };
    });
    if (box.on) box.on('plotly_click', function (ev) { if (ev.points && ev.points[0]) showVariant(ev.points[0].customdata); });
  });

  // ---- Results table
  tab('results', 'Results table', function () { return NROW > 0; }, function (root) {
    var showAll = false, host = h('div');
    var visible = function () { return showAll ? T.cols : T.cols.filter(function (c) { return HIDDEN.indexOf(c) < 0; }); };
    var toggle = h('label', { class: 'field inline small' }, h('input', { type: 'checkbox', onchange: function (e) { showAll = e.target.checked; view.refresh(); } }), 'All columns');
    var pct = function (v) { return typeof v === 'number' ? (v * 100).toFixed(1) + '%' : fmt(v); };
    var view = tableView(host, visible, { toolbar: toggle, filename: 'TALOS_results.csv', onRow: function (i) { showVariant(i + 1); },
      render: { AlleleFrequency: pct } });
    root.appendChild(card(M.events.charAt(0).toUpperCase() + M.events.slice(1) + ' detected: ' + NROW + ' across ' + uniq(col('Sample')).length + ' sample(s)',
      'click a row for details; CSV includes every column', host));
  });

  // ---- Genome view
  tab('genome', 'Genome view', function () { return PANELS.length > 0; }, function (root) {
    var cur = 0, box = plotBox('tall'), legend = h('div', { class: 'legend' }), hidden = {};
    var sel = h('select', { 'aria-label': 'Sample and gene' }, PANELS.map(function (p, i) { return h('option', { value: i, text: p.sample + ' · ' + p.gene + ' (' + p.chrom + ', ' + p.build + ')' }); }));
    sel.addEventListener('change', function () { cur = +sel.value; hidden = {}; draw(); });
    var bar = h('div', { class: 'toolbar' });
    if (PANELS.length > 1) bar.appendChild(h('label', { class: 'field' }, 'Sample · gene', sel));
    bar.appendChild(h('span', { class: 'note', text: 'Introns are compressed; exons are drawn to scale. Zero-length clip calls are shown as ticks.' }));
    root.appendChild(bar);
    root.appendChild(card('Coverage, exons and calls', 'click a call for details', h('div', null, box, legend)));
    function draw() {
      var P = PANELS[cur];
      clear(legend);
      arr(P.variants).forEach(function (v, k) {
        var b = h('button', { type: 'button', 'aria-pressed': hidden[v.vid] ? 'false' : 'true', onclick: function () { hidden[v.vid] = !hidden[v.vid]; draw(); } },
          h('span', { class: 'swatch', style: { background: color(k, resolved()) } }), v.label);
        legend.appendChild(b);
      });
      plot(box, function (t) { return genomeFigure(P, t, hidden); });
      if (box.on && !box._clickBound) { box._clickBound = true; box.on('plotly_click', function (ev) { var p = ev.points && ev.points[0]; if (p && p.data && p.data.meta && p.data.meta.vid) showVariant(p.data.meta.vid); }); }
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
    var q = h('input', { type: 'search', placeholder: 'Search variants…', 'aria-label': 'Search variants' });
    root.appendChild(h('div', { class: 'toolbar' }, h('label', { class: 'field' }, 'Search', q)));
    var list = h('div', { class: 'stack' });
    root.appendChild(list);
    var lastGroup = null, cards = [];
    VARIANTS.forEach(function (v) {
      var g = v.sample + ' · ' + v.gene;
      if (g !== lastGroup) { list.appendChild(h('h2', { class: 'group-title', text: 'Sample: ' + v.sample + ' · ' + v.gene })); lastGroup = g; }
      var body = h('div', null, kv(v.keys, v.values, true));
      if (v.sequence) body.appendChild(h('details', { class: 'seq' }, h('summary', { text: (isAlu ? 'Matched ALU sequence' : 'ITD sequence') + ' (' + v.sequence.length + ' bp)' }), h('pre', { text: v.sequence })));
      var c = card(v.title, v.flag ? 'review flagged' : null, body, { id: 'v-' + v.id, class: 'card vcard' });
      if (v.flag) $('header', c).appendChild(h('span', { class: 'badge flag', text: 'flagged' }));
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

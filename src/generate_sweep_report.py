#!/usr/bin/env python3
"""
Builds a single self-contained HTML report for a parameter sweep produced by
run_param_sweep.sh - charts + full data tables, instead of just summary.csv.

Usage:
    python3 generate_sweep_report.py SWEEP_DIR [--out FILE]

SWEEP_DIR is the sweep_results/PARAM_NAME directory (containing one
subdirectory per tested value). Reuses aggregate_sweep_results.py's data
collection so the report always matches what summary.csv/death_causes_summary.csv
would produce.
"""
import argparse
import json
import os
import sys

import aggregate_sweep_results as agg

CELL_TYPE_LABELS = {
    "zona_pellucida": "Zona pellucida",
    "sox2_high": "Sox2 high",
    "sox17_high": "Sox17 high",
    "loser": "Loser",
    "undifferentiated": "Undifferentiated",
    "differentiated": "Differentiated",
}


def build_report_data(sweep_dir):
    param = os.path.basename(sweep_dir.rstrip("/"))
    results = agg.collect(sweep_dir)
    if not results:
        sys.exit(f"No data found in {sweep_dir}")

    values = list(results.keys())

    cell_types = []
    for cell_type in agg.CELL_TYPES:
        points = []
        for value in values:
            stats = results[value].get(cell_type)
            if stats is None:
                continue
            points.append({
                "value": value,
                "n_orgs": stats["n_orgs"],
                "mean_n_killed": stats["mean_n_killed"],
                "frac_extinct": stats["frac_extinct"],
                "mean_extinction_time": stats["mean_extinction_time"] if stats["mean_extinction_time"] != "" else None,
            })
        if points:
            cell_types.append({
                "key": cell_type,
                "label": CELL_TYPE_LABELS.get(cell_type, cell_type),
                "points": points,
            })

    n_orgs_by_value = {}
    for value in values:
        for cell_type in agg.CELL_TYPES:
            stats = results[value].get(cell_type)
            if stats is not None:
                n_orgs_by_value[value] = stats["n_orgs"]
                break

    death_causes_raw = agg.collect_death_causes(sweep_dir)
    death_causes = None
    if death_causes_raw:
        death_causes = []
        for value, stats in death_causes_raw.items():
            total = stats["mean_lonely"] + stats["mean_signal"]
            death_causes.append({
                "value": value,
                "n_orgs": stats["n_orgs"],
                "mean_lonely": stats["mean_lonely"],
                "mean_signal": stats["mean_signal"],
                "frac_signal": (stats["mean_signal"] / total) if total else None,
            })

    return {
        "param": param,
        "values": values,
        "nOrgsByValue": n_orgs_by_value,
        "cellTypes": cell_types,
        "deathCauses": death_causes,
    }


PAGE_TEMPLATE = r"""<title>Sweep report: {param}</title>
<div class="viz-root">
<style>
.viz-root {{
  color-scheme: light;
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  --page:           #f9f9f7;
  --surface-1:      #fcfcfb;
  --text-primary:   #0b0b0b;
  --text-secondary: #52514e;
  --text-muted:     #898781;
  --gridline:       #e1e0d9;
  --baseline:       #c3c2b7;
  --border:         rgba(11,11,11,0.10);
  --series-lonely:  #2a78d6;
  --series-signal:  #1baf7a;
  --series-killed:  #2a78d6;
  background: var(--page);
  color: var(--text-primary);
  padding: 24px;
}}
@media (prefers-color-scheme: dark) {{
  :root:where(:not([data-theme="light"])) .viz-root {{
    color-scheme: dark;
    --page:           #0d0d0d;
    --surface-1:      #1a1a19;
    --text-primary:   #ffffff;
    --text-secondary: #c3c2b7;
    --text-muted:     #898781;
    --gridline:       #2c2c2a;
    --baseline:       #383835;
    --border:         rgba(255,255,255,0.10);
    --series-lonely:  #3987e5;
    --series-signal:  #199e70;
    --series-killed:  #3987e5;
  }}
}}
:root[data-theme="dark"] .viz-root {{
  color-scheme: dark;
  --page:           #0d0d0d;
  --surface-1:      #1a1a19;
  --text-primary:   #ffffff;
  --text-secondary: #c3c2b7;
  --text-muted:     #898781;
  --gridline:       #2c2c2a;
  --baseline:       #383835;
  --border:         rgba(255,255,255,0.10);
  --series-lonely:  #3987e5;
  --series-signal:  #199e70;
  --series-killed:  #3987e5;
}}

.viz-root * {{ box-sizing: border-box; }}
.report-header h1 {{ font-size: 1.5rem; margin: 0 0 4px; }}
.report-header p {{ color: var(--text-secondary); margin: 0 0 24px; }}
section {{ margin-bottom: 40px; }}
section > h2 {{ font-size: 1.1rem; margin: 0 0 4px; }}
section > .section-note {{ color: var(--text-secondary); font-size: 0.85rem; margin: 0 0 16px; }}

.grid {{
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
  gap: 16px;
}}
.card {{
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 14px 16px 8px;
  position: relative;
  overflow: visible;
}}
.card h3 {{ font-size: 0.9rem; margin: 0 0 8px; font-weight: 600; }}
.card.wide {{ grid-column: 1 / -1; }}

.legend {{ display: flex; gap: 16px; flex-wrap: wrap; margin: 0 0 8px; font-size: 0.8rem; color: var(--text-secondary); }}
.legend-item {{ display: flex; align-items: center; gap: 6px; }}
.legend-key {{ width: 14px; height: 2px; border-radius: 1px; display: inline-block; }}
.legend-key.rect {{ width: 10px; height: 10px; border-radius: 2px; }}

svg.chart {{ width: 100%; height: auto; overflow: visible; display: block; }}
.gridline {{ stroke: var(--gridline); stroke-width: 1; }}
.baseline {{ stroke: var(--baseline); stroke-width: 1; }}
.axis-label {{ fill: var(--text-muted); font-size: 10px; }}
.crosshair {{ stroke: var(--text-muted); stroke-width: 1; pointer-events: none; opacity: 0; }}

.tooltip {{
  position: fixed;
  pointer-events: none;
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 8px 10px;
  font-size: 0.78rem;
  box-shadow: 0 4px 16px rgba(0,0,0,0.15);
  opacity: 0;
  transition: opacity 0.08s ease;
  z-index: 1000;
  max-width: 220px;
}}
.tooltip .tt-x {{ color: var(--text-secondary); margin-bottom: 4px; }}
.tooltip .tt-row {{ display: flex; align-items: center; gap: 6px; }}
.tooltip .tt-key {{ width: 10px; height: 2px; flex: none; border-radius: 1px; }}
.tooltip .tt-name {{ color: var(--text-secondary); flex: 1; }}
.tooltip .tt-val {{ color: var(--text-primary); font-weight: 600; }}

table.data-table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 0.82rem;
  background: var(--surface-1);
  border: 1px solid var(--border);
  border-radius: 8px;
  overflow: hidden;
}}
table.data-table th, table.data-table td {{
  text-align: right;
  padding: 6px 12px;
  border-bottom: 1px solid var(--gridline);
  font-variant-numeric: tabular-nums;
  white-space: nowrap;
}}
table.data-table th:first-child, table.data-table td:first-child {{ text-align: left; font-variant-numeric: normal; }}
table.data-table th {{ color: var(--text-secondary); font-weight: 600; background: var(--page); }}
table.data-table tr:last-child td {{ border-bottom: none; }}
.table-wrap {{ overflow-x: auto; }}
</style>

<div class="report-header">
  <h1>Sweep report: {param}</h1>
  <p>{subtitle}</p>
</div>

<section>
  <h2>Cells killed, by type</h2>
  <p class="section-note">Mean number of cells killed per organism, across the swept values of {param}.</p>
  <div class="grid" id="killed-grid"></div>
</section>

<section>
  <h2>Extinction rate, by type</h2>
  <p class="section-note">Fraction of organisms where that cell type went extinct.</p>
  <div class="grid" id="extinct-grid"></div>
</section>

<section id="cause-section" style="display:none">
  <h2>Cause of death</h2>
  <p class="section-note">Signal-based apoptosis (NeighbourBasedApoptosis) vs. lonely-cell/blastocoel extrusion (ToxictoLonelyCells).</p>
  <div class="grid">
    <div class="card wide"><h3>Mean cells killed, by cause</h3><div id="cause-stackedbar"></div></div>
    <div class="card wide"><h3>Share of deaths from signal apoptosis</h3><div id="cause-fracline"></div></div>
  </div>
</section>

<section>
  <h2>Full data</h2>
  <p class="section-note">Underlying numbers for every chart above.</p>
  <div class="table-wrap" id="killed-table"></div>
  <div style="height:16px"></div>
  <div class="table-wrap" id="cause-table"></div>
</section>

<div class="tooltip" id="tooltip"></div>

</div>

<script type="application/json" id="report-data">{data_json}</script>
<script>
{js_engine}
</script>
"""

JS_ENGINE = r"""
(function () {
  const DATA = JSON.parse(document.getElementById('report-data').textContent);
  const tooltip = document.getElementById('tooltip');
  const root = document.querySelector('.viz-root');
  const cssColor = (name) => getComputedStyle(root).getPropertyValue(name).trim();

  function isNumeric(values) {
    return values.every(v => v !== '' && v !== null && !isNaN(parseFloat(v)) && isFinite(v));
  }

  function niceMax(max) {
    if (max <= 0) return 1;
    const pow = Math.pow(10, Math.floor(Math.log10(max)));
    const n = max / pow;
    let step;
    if (n <= 1) step = 1;
    else if (n <= 2) step = 2;
    else if (n <= 5) step = 5;
    else step = 10;
    return step * pow;
  }

  function fmt(v, digits) {
    if (v === null || v === undefined) return '–';
    if (digits === undefined) digits = 2;
    let s = Number(v).toFixed(digits);
    if (s.indexOf('.') !== -1) s = s.replace(/0+$/, '').replace(/\.$/, '');
    return s;
  }

  function svgEl(tag, attrs) {
    const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const k in attrs) el.setAttribute(k, attrs[k]);
    return el;
  }

  function setText(el, text) {
    el.textContent = text;
    return el;
  }

  // X scale: numeric (linear) if every label parses as a number, else evenly-spaced/ordinal.
  function makeXScale(labels, plotW, padL) {
    const numeric = isNumeric(labels);
    if (numeric) {
      const nums = labels.map(Number);
      const min = Math.min(...nums), max = Math.max(...nums);
      const span = (max - min) || 1;
      return { x: (lbl) => padL + ((Number(lbl) - min) / span) * plotW };
    }
    const n = labels.length;
    const step = n > 1 ? plotW / (n - 1) : 0;
    return { x: (lbl) => padL + labels.indexOf(lbl) * step };
  }

  function showTooltip(evt, xLabel, rows) {
    tooltip.innerHTML = '';
    tooltip.appendChild(setText(document.createElement('div'), xLabel)).className = 'tt-x';
    rows.forEach(r => {
      const row = document.createElement('div');
      row.className = 'tt-row';
      const key = document.createElement('span');
      key.className = 'tt-key';
      key.style.background = r.color;
      const name = document.createElement('span');
      name.className = 'tt-name';
      name.textContent = r.name;
      const val = document.createElement('span');
      val.className = 'tt-val';
      val.textContent = r.value;
      row.appendChild(key);
      row.appendChild(name);
      row.appendChild(val);
      tooltip.appendChild(row);
    });
    const pad = 14;
    tooltip.style.opacity = '1';
    tooltip.style.left = (evt.clientX + pad) + 'px';
    tooltip.style.top = (evt.clientY + pad) + 'px';
    const rect = tooltip.getBoundingClientRect();
    let left = evt.clientX + pad, top = evt.clientY + pad;
    if (rect.right > window.innerWidth) left = evt.clientX - rect.width - pad;
    if (rect.bottom > window.innerHeight) top = evt.clientY - rect.height - pad;
    tooltip.style.left = Math.max(4, left) + 'px';
    tooltip.style.top = Math.max(4, top) + 'px';
  }

  function hideTooltip() { tooltip.style.opacity = '0'; }

  function renderLegend(container, series) {
    if (series.length < 2) return;
    const legend = document.createElement('div');
    legend.className = 'legend';
    series.forEach(s => {
      const item = document.createElement('div');
      item.className = 'legend-item';
      const key = document.createElement('span');
      key.className = 'legend-key' + (s.rect ? ' rect' : '');
      key.style.background = s.color;
      item.appendChild(key);
      const label = document.createElement('span');
      label.textContent = s.name;
      item.appendChild(label);
      legend.appendChild(item);
    });
    container.appendChild(legend);
  }

  // opts: {xLabels, xTitle, series:[{name,color,values}], yFormat(v), yDigits, showLegend}
  function renderLineChart(container, opts) {
    const W = 400, H = 190;
    const padL = 34, padR = 12, padT = 10, padB = 22;
    const plotW = W - padL - padR, plotH = H - padT - padB;

    const allVals = opts.series.flatMap(s => s.values.filter(v => v !== null && v !== undefined));
    const yMax = niceMax(Math.max(0.0001, ...allVals));
    const yScale = (v) => padT + plotH - (v / yMax) * plotH;
    const xs = makeXScale(opts.xLabels, plotW, padL);

    if (opts.showLegend !== false) renderLegend(container, opts.series);

    const svg = svgEl('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}` });

    const steps = 4;
    for (let i = 0; i <= steps; i++) {
      const v = (yMax / steps) * i;
      const y = yScale(v);
      svg.appendChild(svgEl('line', { class: 'gridline', x1: padL, x2: W - padR, y1: y, y2: y }));
      const label = svgEl('text', { class: 'axis-label', x: padL - 6, y: y + 3, 'text-anchor': 'end' });
      label.textContent = fmt(v, opts.yDigits !== undefined ? opts.yDigits : 1);
      svg.appendChild(label);
    }
    svg.appendChild(svgEl('line', { class: 'baseline', x1: padL, x2: W - padR, y1: padT + plotH, y2: padT + plotH }));

    opts.xLabels.forEach((lbl) => {
      const x = xs.x(lbl);
      const label = svgEl('text', { class: 'axis-label', x: x, y: H - 6, 'text-anchor': 'middle' });
      label.textContent = lbl;
      svg.appendChild(label);
    });

    const crosshair = svgEl('line', { class: 'crosshair', x1: padL, x2: padL, y1: padT, y2: padT + plotH });
    svg.appendChild(crosshair);

    const markers = [];
    opts.series.forEach(s => {
      const pts = opts.xLabels.map((lbl, i) => {
        const v = s.values[i];
        return (v === null || v === undefined) ? null : [xs.x(lbl), yScale(v)];
      });
      const validPts = pts.filter(p => p !== null);
      if (validPts.length) {
        const d = validPts.map((p, i) => (i === 0 ? 'M' : 'L') + p[0].toFixed(2) + ' ' + p[1].toFixed(2)).join(' ');
        svg.appendChild(svgEl('path', { d, fill: 'none', stroke: s.color, 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }));
      }
      pts.forEach((p, i) => {
        if (!p) return;
        svg.appendChild(svgEl('circle', { cx: p[0], cy: p[1], r: 4, fill: s.color, stroke: cssColor('--surface-1'), 'stroke-width': 2 }));
        markers.push({ x: p[0], xLabel: opts.xLabels[i], name: s.name, color: s.color, value: s.values[i] });
      });
    });

    svg.addEventListener('pointermove', (evt) => {
      const rect = svg.getBoundingClientRect();
      const localX = (evt.clientX - rect.left) * (W / rect.width);
      let nearestLabel = opts.xLabels[0], nearestDist = Infinity, nearestX = padL;
      opts.xLabels.forEach(lbl => {
        const x = xs.x(lbl);
        const dist = Math.abs(x - localX);
        if (dist < nearestDist) { nearestDist = dist; nearestLabel = lbl; nearestX = x; }
      });
      crosshair.setAttribute('x1', nearestX);
      crosshair.setAttribute('x2', nearestX);
      crosshair.style.opacity = '1';
      const rows = markers.filter(m => m.xLabel === nearestLabel).map(m => ({
        name: m.name, color: m.color, value: opts.yFormat ? opts.yFormat(m.value) : fmt(m.value)
      }));
      const xLabel = (opts.xTitle ? opts.xTitle + ' = ' : '') + nearestLabel;
      showTooltip(evt, xLabel, rows);
    });
    svg.addEventListener('pointerleave', () => { crosshair.style.opacity = '0'; hideTooltip(); });

    container.appendChild(svg);
  }

  // opts: {xLabels, xTitle, series:[{name,color,values}], yFormat(v)} - series stacked bottom-up in given order.
  function renderStackedBarChart(container, opts) {
    const W = 640, H = 220;
    const padL = 40, padR = 12, padT = 10, padB = 22;
    const plotW = W - padL - padR, plotH = H - padT - padB;

    const totals = opts.xLabels.map((_, i) => opts.series.reduce((sum, s) => sum + (s.values[i] || 0), 0));
    const yMax = niceMax(Math.max(0.0001, ...totals));
    const yScale = (v) => (v / yMax) * plotH;

    renderLegend(container, opts.series);

    const svg = svgEl('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}` });

    const steps = 4;
    for (let i = 0; i <= steps; i++) {
      const v = (yMax / steps) * i;
      const y = padT + plotH - yScale(v);
      svg.appendChild(svgEl('line', { class: 'gridline', x1: padL, x2: W - padR, y1: y, y2: y }));
      const label = svgEl('text', { class: 'axis-label', x: padL - 6, y: y + 3, 'text-anchor': 'end' });
      label.textContent = fmt(v, 1);
      svg.appendChild(label);
    }
    svg.appendChild(svgEl('line', { class: 'baseline', x1: padL, x2: W - padR, y1: padT + plotH, y2: padT + plotH }));

    const n = opts.xLabels.length;
    const slot = plotW / n;
    const barW = Math.min(28, slot * 0.55);
    const gap = 2;
    const bars = []; // for hover lookup

    opts.xLabels.forEach((lbl, i) => {
      const cx = padL + slot * (i + 0.5);
      const label = svgEl('text', { class: 'axis-label', x: cx, y: H - 6, 'text-anchor': 'middle' });
      label.textContent = lbl;
      svg.appendChild(label);

      let yCursor = padT + plotH;
      const segRows = [];
      // The topmost drawn segment gets the rounded data-end, not necessarily
      // the last series in the array - a zero-valued top series must not
      // leave the visible top of the bar square.
      let lastVisibleIdx = -1;
      opts.series.forEach((s, si) => { if ((s.values[i] || 0) > 0) lastVisibleIdx = si; });
      opts.series.forEach((s, si) => {
        const v = s.values[i] || 0;
        const h = Math.max(0, yScale(v) - (si < lastVisibleIdx ? gap / 2 : 0));
        if (h <= 0) return;
        const yTop = yCursor - h;
        const isLast = si === lastVisibleIdx;
        const r = isLast ? 4 : 0;
        const path = isLast
          ? `M ${cx - barW / 2} ${yCursor} L ${cx - barW / 2} ${yTop + r} Q ${cx - barW / 2} ${yTop} ${cx - barW / 2 + r} ${yTop} `
            + `L ${cx + barW / 2 - r} ${yTop} Q ${cx + barW / 2} ${yTop} ${cx + barW / 2} ${yTop + r} L ${cx + barW / 2} ${yCursor} Z`
          : `M ${cx - barW / 2} ${yCursor} L ${cx - barW / 2} ${yTop} L ${cx + barW / 2} ${yTop} L ${cx + barW / 2} ${yCursor} Z`;
        svg.appendChild(svgEl('path', { d: path, fill: s.color }));
        segRows.push({ name: s.name, color: s.color, value: opts.yFormat ? opts.yFormat(v) : fmt(v) });
        yCursor = yTop - gap;
      });
      bars.push({ cx, halfW: slot / 2, xLabel: lbl, rows: segRows.reverse() });
    });

    svg.addEventListener('pointermove', (evt) => {
      const rect = svg.getBoundingClientRect();
      const localX = (evt.clientX - rect.left) * (W / rect.width);
      const bar = bars.find(b => Math.abs(b.cx - localX) <= b.halfW) || bars[bars.length - 1];
      const xLabel = (opts.xTitle ? opts.xTitle + ' = ' : '') + bar.xLabel;
      showTooltip(evt, xLabel, bar.rows);
    });
    svg.addEventListener('pointerleave', hideTooltip);

    container.appendChild(svg);
  }

  function renderTable(container, columns, rows) {
    const table = document.createElement('table');
    table.className = 'data-table';
    const thead = document.createElement('thead');
    const headRow = document.createElement('tr');
    columns.forEach(c => {
      const th = document.createElement('th');
      th.textContent = c.label;
      headRow.appendChild(th);
    });
    thead.appendChild(headRow);
    table.appendChild(thead);
    const tbody = document.createElement('tbody');
    rows.forEach(r => {
      const tr = document.createElement('tr');
      columns.forEach(c => {
        const td = document.createElement('td');
        td.textContent = c.get(r);
        tr.appendChild(td);
      });
      tbody.appendChild(tr);
    });
    table.appendChild(tbody);
    container.appendChild(table);
  }

  // --- Header subtitle already rendered server-side ---

  // --- Cells killed / extinction rate small multiples ---
  const killedGrid = document.getElementById('killed-grid');
  const extinctGrid = document.getElementById('extinct-grid');
  const xLabels = DATA.values;

  DATA.cellTypes.forEach(ct => {
    const killedCard = document.createElement('div');
    killedCard.className = 'card';
    killedCard.appendChild(setText(document.createElement('h3'), ct.label));
    renderLineChart(killedCard, {
      xLabels, xTitle: DATA.param,
      series: [{ name: ct.label, color: cssColor('--series-killed'), values: ct.points.map(p => p.mean_n_killed) }],
      showLegend: false,
      yFormat: (v) => fmt(v, 2),
    });
    killedGrid.appendChild(killedCard);

    const extinctCard = document.createElement('div');
    extinctCard.className = 'card';
    extinctCard.appendChild(setText(document.createElement('h3'), ct.label));
    renderLineChart(extinctCard, {
      xLabels, xTitle: DATA.param,
      series: [{ name: ct.label, color: cssColor('--series-killed'), values: ct.points.map(p => p.frac_extinct * 100) }],
      showLegend: false,
      yFormat: (v) => fmt(v, 0) + '%',
      yDigits: 0,
    });
    extinctGrid.appendChild(extinctCard);
  });

  // --- Cause of death ---
  if (DATA.deathCauses) {
    document.getElementById('cause-section').style.display = '';
    const dcLabels = DATA.deathCauses.map(d => d.value);
    const lonelyColor = cssColor('--series-lonely');
    const signalColor = cssColor('--series-signal');

    renderStackedBarChart(document.getElementById('cause-stackedbar'), {
      xLabels: dcLabels, xTitle: DATA.param,
      series: [
        { name: 'Lonely / blastocoel', color: lonelyColor, rect: true, values: DATA.deathCauses.map(d => d.mean_lonely) },
        { name: 'Signal apoptosis', color: signalColor, rect: true, values: DATA.deathCauses.map(d => d.mean_signal) },
      ],
      yFormat: (v) => fmt(v, 2),
    });

    renderLineChart(document.getElementById('cause-fracline'), {
      xLabels: dcLabels, xTitle: DATA.param,
      series: [{ name: 'Share from signal apoptosis', color: signalColor, values: DATA.deathCauses.map(d => d.frac_signal === null ? null : d.frac_signal * 100) }],
      showLegend: false,
      yFormat: (v) => fmt(v, 0) + '%',
      yDigits: 0,
    });
  }

  // --- Full data tables ---
  const killedRows = [];
  DATA.cellTypes.forEach(ct => {
    ct.points.forEach(p => killedRows.push({ label: ct.label, ...p }));
  });
  renderTable(document.getElementById('killed-table'), [
    { label: 'Cell type', get: r => r.label },
    { label: DATA.param, get: r => r.value },
    { label: 'n orgs', get: r => r.n_orgs },
    { label: 'Mean killed', get: r => fmt(r.mean_n_killed, 2) },
    { label: 'Frac. extinct', get: r => fmt(r.frac_extinct * 100, 0) + '%' },
    { label: 'Mean extinction t', get: r => r.mean_extinction_time === null ? '–' : fmt(r.mean_extinction_time, 0) },
  ], killedRows);

  if (DATA.deathCauses) {
    renderTable(document.getElementById('cause-table'), [
      { label: DATA.param, get: r => r.value },
      { label: 'n orgs', get: r => r.n_orgs },
      { label: 'Mean lonely killed', get: r => fmt(r.mean_lonely, 2) },
      { label: 'Mean signal killed', get: r => fmt(r.mean_signal, 2) },
      { label: 'Share signal', get: r => r.frac_signal === null ? '–' : fmt(r.frac_signal * 100, 0) + '%' },
    ], DATA.deathCauses);
  }
})();
"""


def render_html(data):
    param = data["param"]
    n_orgs_values = sorted(set(data["nOrgsByValue"].values()))
    n_orgs_str = str(n_orgs_values[0]) if len(n_orgs_values) == 1 else "/".join(str(n) for n in n_orgs_values)
    subtitle = f"{len(data['values'])} values tested, {n_orgs_str} organisms per value."
    return PAGE_TEMPLATE.format(
        param=param,
        subtitle=subtitle,
        data_json=json.dumps(data).replace("</", "<\\/"),
        js_engine=JS_ENGINE,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sweep_dir", help="sweep_results/PARAM_NAME directory produced by run_param_sweep.sh")
    parser.add_argument("--out", default=None, help="HTML file to write (default: SWEEP_DIR/report.html)")
    args = parser.parse_args()

    sweep_dir = args.sweep_dir.rstrip("/")
    out_path = args.out or os.path.join(sweep_dir, "report.html")

    data = build_report_data(sweep_dir)
    html = render_html(data)
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Report written to {out_path}")


if __name__ == "__main__":
    main()

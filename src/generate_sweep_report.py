#!/usr/bin/env python3
"""
Builds a single self-contained HTML report for a parameter sweep produced by
run_param_sweep.py.

The report shows one quantity: the proportion of a cell type's cells that
were killed by the end of the simulation, averaged across the organisms
simulated at each swept value (or parameter combination). Two controls pick
the representation:
  - which cell type(s) to show (one small-multiple panel per type)
  - which death cause to show: combined (every death), cell sorting
    (lonely-cell/blastocoel extrusion), or apoptosis (neighbour-based
    signalling)
Each panel is a line chart when one parameter was swept, or a phase-diagram
heatmap when exactly two were. A full data table underneath always lists
every cell type/value combination and both causes, regardless of the
controls above.

Usage:
    python3 generate_sweep_report.py SWEEP_DIR [--out FILE]

SWEEP_DIR is the sweep_results/PARAM_NAME (or PARAM1+PARAM2) directory
produced by run_param_sweep.py. Reuses aggregate_sweep_results.py's
collect_death_causes() so the report always matches death_causes_summary.csv.
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
    "undifferentiated": "Undifferentiated",
    "differentiated": "Differentiated",
    "total": "Total (whole organism)",
}


def combined_frac(stats):
    """Proportion killed by either cause combined, or None if neither
    fraction is available (sweep run before per-type cause tracking)."""
    lonely, signal = stats["mean_frac_lonely_of_born"], stats["mean_frac_signal_of_born"]
    if lonely == "" and signal == "":
        return None
    return (lonely or 0) + (signal or 0)


def build_grid(results, parsed, x_values, y_values, param_x, param_y, get_value):
    """Pivots {combo_key: by_type} into a 2D grid (rows = y values, columns =
    x values) of get_value(by_type), or None if no combination has a value."""
    x_index = {v: i for i, v in enumerate(x_values)}
    y_index = {v: i for i, v in enumerate(y_values)}
    grid = [[None] * len(x_values) for _ in y_values]
    any_value = False
    for key, by_type in results.items():
        pairs = parsed.get(key)
        if pairs is None:
            continue
        value = get_value(by_type)
        if value is None:
            continue
        grid[y_index[pairs[param_y]]][x_index[pairs[param_x]]] = value
        any_value = True
    return grid if any_value else None


def build_report_data(sweep_dir):
    param = os.path.basename(sweep_dir.rstrip("/"))
    results = agg.collect_death_causes(sweep_dir)
    if not results:
        sys.exit(
            f"No per-cell-type death-cause data found in {sweep_dir}.\n"
            "This report needs death_causes-org-*.dat files with per-cell-type "
            "columns - re-run the sweep with the current embryo_multi."
        )

    values = list(results.keys())
    axes = agg.phase_diagram_axes(results)

    cell_types = []
    for cell_type in agg.DEATH_CAUSE_TYPES:
        points = []
        for value in values:
            stats = results[value].get(cell_type)
            if stats is None:
                continue
            lonely = stats["mean_frac_lonely_of_born"]
            signal = stats["mean_frac_signal_of_born"]
            points.append({
                "value": value,
                "n_orgs": stats["n_orgs"],
                "fracSorting": lonely if lonely != "" else None,
                "fracApoptosis": signal if signal != "" else None,
                "fracCombined": combined_frac(stats),
            })
        if points:
            cell_types.append({"key": cell_type, "label": CELL_TYPE_LABELS.get(cell_type, cell_type), "points": points})

    phase_diagram = None
    if axes:
        param_x, param_y = axes
        parsed = {key: dict(agg.parse_combo_key(key)) for key in results}
        x_values = sorted({p[param_x] for p in parsed.values()}, key=agg.sort_key)
        y_values = sorted({p[param_y] for p in parsed.values()}, key=agg.sort_key)

        grids_by_type = {}
        for ct in cell_types:
            key = ct["key"]

            def get(field, by_type):
                stats = by_type.get(key)
                if stats is None or stats[field] == "":
                    return None
                return stats[field]

            grids = {}
            sorting_grid = build_grid(results, parsed, x_values, y_values, param_x, param_y,
                                       lambda by_type: get("mean_frac_lonely_of_born", by_type))
            apoptosis_grid = build_grid(results, parsed, x_values, y_values, param_x, param_y,
                                         lambda by_type: get("mean_frac_signal_of_born", by_type))
            combined_grid = build_grid(results, parsed, x_values, y_values, param_x, param_y,
                                        lambda by_type, key=key: (
                                            combined_frac(by_type[key]) if by_type.get(key) is not None else None
                                        ))
            if sorting_grid is not None:
                grids["sorting"] = sorting_grid
            if apoptosis_grid is not None:
                grids["apoptosis"] = apoptosis_grid
            if combined_grid is not None:
                grids["combined"] = combined_grid
            if grids:
                grids_by_type[key] = grids

        phase_diagram = {
            "paramX": param_x, "paramY": param_y,
            "xValues": x_values, "yValues": y_values,
            "grids": grids_by_type,
        }

    return {
        "param": param,
        "values": values,
        "cellTypes": cell_types,
        "phaseDiagram": phase_diagram,
    }


PAGE_TEMPLATE = r"""<title>Sweep report: {param}</title>
<div class="viz-root">
<style>
.viz-root {{
  color-scheme: light;
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  --page:              #f9f9f7;
  --surface-1:         #fcfcfb;
  --text-primary:      #0b0b0b;
  --text-secondary:    #52514e;
  --text-muted:        #898781;
  --gridline:          #e1e0d9;
  --baseline:          #c3c2b7;
  --border:            rgba(11,11,11,0.10);
  --cause-combined:    #4a3aa7;
  --cause-sorting:     #2a78d6;
  --cause-apoptosis:   #1baf7a;
  --heat-low:          #cde2fb;
  --heat-high:         #0d366b;
  background: var(--page);
  color: var(--text-primary);
  padding: 24px;
}}
@media (prefers-color-scheme: dark) {{
  :root:where(:not([data-theme="light"])) .viz-root {{
    color-scheme: dark;
    --page:              #0d0d0d;
    --surface-1:         #1a1a19;
    --text-primary:      #ffffff;
    --text-secondary:    #c3c2b7;
    --text-muted:        #898781;
    --gridline:          #2c2c2a;
    --baseline:          #383835;
    --border:            rgba(255,255,255,0.10);
    --cause-combined:    #9085e9;
    --cause-sorting:     #3987e5;
    --cause-apoptosis:   #199e70;
    --heat-low:          #184f95;
    --heat-high:         #b7d3f6;
  }}
}}
:root[data-theme="dark"] .viz-root {{
  color-scheme: dark;
  --page:              #0d0d0d;
  --surface-1:         #1a1a19;
  --text-primary:      #ffffff;
  --text-secondary:    #c3c2b7;
  --text-muted:        #898781;
  --gridline:          #2c2c2a;
  --baseline:          #383835;
  --border:            rgba(255,255,255,0.10);
  --cause-combined:    #9085e9;
  --cause-sorting:     #3987e5;
  --cause-apoptosis:   #199e70;
  --heat-low:          #184f95;
  --heat-high:         #b7d3f6;
}}

.viz-root * {{ box-sizing: border-box; }}
.report-header h1 {{ font-size: 1.5rem; margin: 0 0 4px; }}
.report-header p {{ color: var(--text-secondary); margin: 0 0 24px; }}
section {{ margin-bottom: 32px; }}
section > h2 {{ font-size: 1.1rem; margin: 0 0 4px; }}
section > .section-note {{ color: var(--text-secondary); font-size: 0.85rem; margin: 0 0 16px; }}

.controls {{ display: flex; flex-direction: column; gap: 10px; margin-bottom: 28px; }}
.control-row {{ display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }}
.control-row .control-label {{ color: var(--text-secondary); font-size: 0.85rem; margin-right: 2px; }}
.chip {{
  display: inline-flex; align-items: center; gap: 6px;
  padding: 5px 12px; border: 1px solid var(--border); border-radius: 999px;
  font-size: 0.82rem; cursor: pointer; user-select: none; background: var(--surface-1);
}}
.chip input {{ margin: 0; }}
.chip.unchecked {{ color: var(--text-muted); }}
.chip-link {{
  font-size: 0.78rem; color: var(--text-secondary); cursor: pointer;
  text-decoration: underline; background: none; border: none; padding: 0 4px;
}}

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
.card .card-subnote {{ font-size: 0.76rem; color: var(--text-muted); margin: -6px 0 8px; }}

.heat-scale {{ display: flex; align-items: center; gap: 6px; font-size: 0.72rem; color: var(--text-secondary); margin-top: 4px; }}
.heat-scale-bar {{ width: 80px; height: 8px; border-radius: 4px; background: linear-gradient(to right, var(--heat-low), var(--heat-high)); }}

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
.tooltip .tt-key {{ width: 10px; height: 10px; flex: none; border-radius: 2px; }}
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

<div class="controls">
  <div class="control-row" id="type-filter">
    <span class="control-label">Cell type(s) to show:</span>
    <button class="chip-link" id="filter-all" type="button">All</button>
    <button class="chip-link" id="filter-none" type="button">None</button>
  </div>
  <div class="control-row" id="cause-filter">
    <span class="control-label">Death cause to show:</span>
  </div>
</div>

<section>
  <h2>Proportion of cells killed by the end of the simulation</h2>
  <p class="section-note" id="chart-section-note"></p>
  <div class="grid" id="chart-grid"></div>
</section>

<section>
  <h2>Full data</h2>
  <p class="section-note">Every cell type/value combination, both causes - independent of the controls above.</p>
  <div class="table-wrap" id="data-table"></div>
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

  const CAUSES = [
    { key: 'combined', label: 'Combined (all deaths)', field: 'fracCombined', color: '--cause-combined' },
    { key: 'sorting', label: 'Cell sorting (lonely / blastocoel)', field: 'fracSorting', color: '--cause-sorting' },
    { key: 'apoptosis', label: 'Apoptosis (signal)', field: 'fracApoptosis', color: '--cause-apoptosis' },
  ];
  let selectedCause = CAUSES[0];

  function fmtPct(v) {
    if (v === null || v === undefined) return '–';
    return (v * 100).toFixed(1).replace(/\.0$/, '') + '%';
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

  function svgEl(tag, attrs) {
    const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const k in attrs) el.setAttribute(k, attrs[k]);
    return el;
  }

  function isNumeric(values) {
    return values.every(v => v !== '' && v !== null && !isNaN(parseFloat(v)) && isFinite(v));
  }

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
    tooltip.appendChild(Object.assign(document.createElement('div'), { textContent: xLabel, className: 'tt-x' }));
    rows.forEach(r => {
      const row = document.createElement('div');
      row.className = 'tt-row';
      const key = document.createElement('span');
      key.className = 'tt-key';
      key.style.background = r.color;
      const name = Object.assign(document.createElement('span'), { className: 'tt-name', textContent: r.name });
      const val = Object.assign(document.createElement('span'), { className: 'tt-val', textContent: r.value });
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

  // opts: {xLabels, xTitle, values, color, name}
  function renderLineChart(container, opts) {
    const W = 400, H = 190;
    const padL = 34, padR = 12, padT = 10, padB = 22;
    const plotW = W - padL - padR, plotH = H - padT - padB;

    const present = opts.values.filter(v => v !== null && v !== undefined);
    const yMax = niceMax(Math.max(0.0001, ...present, 0.05));
    const yScale = (v) => padT + plotH - (v / yMax) * plotH;
    const xs = makeXScale(opts.xLabels, plotW, padL);

    const svg = svgEl('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}` });

    const steps = 4;
    for (let i = 0; i <= steps; i++) {
      const v = (yMax / steps) * i;
      const y = yScale(v);
      svg.appendChild(svgEl('line', { class: 'gridline', x1: padL, x2: W - padR, y1: y, y2: y }));
      const label = svgEl('text', { class: 'axis-label', x: padL - 6, y: y + 3, 'text-anchor': 'end' });
      label.textContent = fmtPct(v);
      svg.appendChild(label);
    }
    svg.appendChild(svgEl('line', { class: 'baseline', x1: padL, x2: W - padR, y1: padT + plotH, y2: padT + plotH }));

    const maxTicks = 6;
    const tickStep = opts.xLabels.length > maxTicks ? Math.ceil(opts.xLabels.length / maxTicks) : 1;
    opts.xLabels.forEach((lbl, i) => {
      if (i % tickStep !== 0 && i !== opts.xLabels.length - 1) return;
      const x = xs.x(lbl);
      const label = svgEl('text', { class: 'axis-label', x: x, y: H - 6, 'text-anchor': 'middle' });
      label.textContent = lbl;
      svg.appendChild(label);
    });

    const crosshair = svgEl('line', { class: 'crosshair', x1: padL, x2: padL, y1: padT, y2: padT + plotH });
    svg.appendChild(crosshair);

    const pts = opts.xLabels.map((lbl, i) => {
      const v = opts.values[i];
      return (v === null || v === undefined) ? null : [xs.x(lbl), yScale(v)];
    });
    const validPts = pts.filter(p => p !== null);
    if (validPts.length) {
      const d = validPts.map((p, i) => (i === 0 ? 'M' : 'L') + p[0].toFixed(2) + ' ' + p[1].toFixed(2)).join(' ');
      svg.appendChild(svgEl('path', { d, fill: 'none', stroke: opts.color, 'stroke-width': 2, 'stroke-linejoin': 'round', 'stroke-linecap': 'round' }));
    }
    pts.forEach((p, i) => {
      if (!p) return;
      svg.appendChild(svgEl('circle', { cx: p[0], cy: p[1], r: 4, fill: opts.color, stroke: cssColor('--surface-1'), 'stroke-width': 2 }));
    });

    svg.addEventListener('pointermove', (evt) => {
      const rect = svg.getBoundingClientRect();
      const localX = (evt.clientX - rect.left) * (W / rect.width);
      let nearestLabel = opts.xLabels[0], nearestDist = Infinity, nearestX = padL, nearestI = 0;
      opts.xLabels.forEach((lbl, i) => {
        const x = xs.x(lbl);
        const dist = Math.abs(x - localX);
        if (dist < nearestDist) { nearestDist = dist; nearestLabel = lbl; nearestX = x; nearestI = i; }
      });
      const v = opts.values[nearestI];
      if (v === null || v === undefined) { hideTooltip(); return; }
      crosshair.setAttribute('x1', nearestX);
      crosshair.setAttribute('x2', nearestX);
      crosshair.style.opacity = '1';
      const xLabel = (opts.xTitle ? opts.xTitle + ' = ' : '') + nearestLabel;
      showTooltip(evt, xLabel, [{ name: opts.name, color: opts.color, value: fmtPct(v) }]);
    });
    svg.addEventListener('pointerleave', () => { crosshair.style.opacity = '0'; hideTooltip(); });

    container.appendChild(svg);
  }

  function lerpColor(hexLow, hexHigh, t) {
    t = Math.max(0, Math.min(1, t));
    const parse = (hex) => [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
    const [r1, g1, b1] = parse(hexLow);
    const [r2, g2, b2] = parse(hexHigh);
    const r = Math.round(r1 + (r2 - r1) * t);
    const g = Math.round(g1 + (g2 - g1) * t);
    const b = Math.round(b1 + (b2 - b1) * t);
    return `rgb(${r}, ${g}, ${b})`;
  }

  // opts: {xLabels, yLabels, xTitle, yTitle, grid (grid[yi][xi] = number|null)}
  function renderHeatmap(container, opts) {
    const cellW = 56, cellH = 32;
    const padL = 70, padR = 12, padT = 10, padB = 34;
    const nX = opts.xLabels.length, nY = opts.yLabels.length;
    const W = padL + cellW * nX + padR;
    const H = padT + cellH * nY + padB;
    const gap = 2;

    const values = opts.grid.flat().filter(v => v !== null && v !== undefined);
    const min = values.length ? Math.min(...values) : 0;
    const max = values.length ? Math.max(...values) : 1;
    const span = (max - min) || 1;
    const heatLow = cssColor('--heat-low');
    const heatHigh = cssColor('--heat-high');
    const noDataColor = cssColor('--gridline');

    const svg = svgEl('svg', { class: 'chart', viewBox: `0 0 ${W} ${H}` });

    opts.yLabels.forEach((lbl, yi) => {
      const label = svgEl('text', { class: 'axis-label', x: padL - 8, y: padT + yi * cellH + cellH / 2 + 3, 'text-anchor': 'end' });
      label.textContent = lbl;
      svg.appendChild(label);
    });
    opts.xLabels.forEach((lbl, xi) => {
      const label = svgEl('text', { class: 'axis-label', x: padL + xi * cellW + cellW / 2, y: H - padB + 16, 'text-anchor': 'middle' });
      label.textContent = lbl;
      svg.appendChild(label);
    });

    const cells = [];
    opts.yLabels.forEach((yLbl, yi) => {
      opts.xLabels.forEach((xLbl, xi) => {
        const v = opts.grid[yi][xi];
        const x = padL + xi * cellW + gap / 2, y = padT + yi * cellH + gap / 2;
        const w = cellW - gap, h = cellH - gap;
        const fill = (v === null || v === undefined) ? noDataColor : lerpColor(heatLow, heatHigh, (v - min) / span);
        svg.appendChild(svgEl('rect', { x, y, width: w, height: h, rx: 3, fill }));
        cells.push({ x, y, w, h, xLbl, yLbl, v });
      });
    });

    const highlight = svgEl('rect', {
      class: 'crosshair', x: 0, y: 0, width: cellW - gap, height: cellH - gap, rx: 3,
      fill: 'none', stroke: cssColor('--text-primary'), 'stroke-width': 2,
    });
    svg.appendChild(highlight);

    svg.addEventListener('pointermove', (evt) => {
      const rect = svg.getBoundingClientRect();
      const localX = (evt.clientX - rect.left) * (W / rect.width);
      const localY = (evt.clientY - rect.top) * (H / rect.height);
      const cell = cells.find(c => localX >= c.x && localX <= c.x + c.w && localY >= c.y && localY <= c.y + c.h);
      if (!cell) { highlight.style.opacity = '0'; hideTooltip(); return; }
      highlight.setAttribute('x', cell.x);
      highlight.setAttribute('y', cell.y);
      highlight.style.opacity = '1';
      const xLabel = `${opts.xTitle} = ${cell.xLbl}, ${opts.yTitle} = ${cell.yLbl}`;
      const value = cell.v === null || cell.v === undefined ? 'no data' : fmtPct(cell.v);
      showTooltip(evt, xLabel, [{ name: 'Value', color: cell.v === null ? noDataColor : heatHigh, value }]);
    });
    svg.addEventListener('pointerleave', () => { highlight.style.opacity = '0'; hideTooltip(); });

    container.appendChild(svg);

    if (values.length) {
      const scale = document.createElement('div');
      scale.className = 'heat-scale';
      const barEl = document.createElement('div');
      barEl.className = 'heat-scale-bar';
      scale.appendChild(Object.assign(document.createElement('span'), { textContent: fmtPct(min) }));
      scale.appendChild(barEl);
      scale.appendChild(Object.assign(document.createElement('span'), { textContent: fmtPct(max) }));
      container.appendChild(scale);
    }
  }

  function renderTable(container, cellTypes) {
    const table = document.createElement('table');
    table.className = 'data-table';
    const thead = document.createElement('thead');
    const headRow = document.createElement('tr');
    ['Cell type', DATA.param, 'n orgs', 'Killed - cell sorting', 'Killed - apoptosis', 'Killed - combined'].forEach(text => {
      const th = document.createElement('th');
      th.textContent = text;
      headRow.appendChild(th);
    });
    thead.appendChild(headRow);
    table.appendChild(thead);
    const tbody = document.createElement('tbody');
    cellTypes.forEach(ct => {
      ct.points.forEach(p => {
        const tr = document.createElement('tr');
        [ct.label, p.value, p.n_orgs, fmtPct(p.fracSorting), fmtPct(p.fracApoptosis), fmtPct(p.fracCombined)].forEach(text => {
          const td = document.createElement('td');
          td.textContent = text;
          tr.appendChild(td);
        });
        tbody.appendChild(tr);
      });
    });
    table.appendChild(tbody);
    container.innerHTML = '';
    container.appendChild(table);
  }

  // --- Controls: cell-type chips + cause radio ---
  const pd = DATA.phaseDiagram;
  const selectedTypes = new Set(DATA.cellTypes.map(ct => ct.key));

  const typeFilterRow = document.getElementById('type-filter');
  const allBtn = document.getElementById('filter-all');
  const noneBtn = document.getElementById('filter-none');
  const typeCheckboxes = [];
  DATA.cellTypes.forEach(ct => {
    const chip = document.createElement('label');
    chip.className = 'chip';
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.checked = true;
    box.addEventListener('change', () => {
      if (box.checked) selectedTypes.add(ct.key); else selectedTypes.delete(ct.key);
      chip.classList.toggle('unchecked', !box.checked);
      renderAll();
    });
    chip.appendChild(box);
    chip.appendChild(Object.assign(document.createElement('span'), { textContent: ct.label }));
    typeFilterRow.insertBefore(chip, allBtn);
    typeCheckboxes.push(box);
  });
  allBtn.addEventListener('click', () => typeCheckboxes.forEach(b => { b.checked = true; b.dispatchEvent(new Event('change')); }));
  noneBtn.addEventListener('click', () => typeCheckboxes.forEach(b => { b.checked = false; b.dispatchEvent(new Event('change')); }));

  const causeFilterRow = document.getElementById('cause-filter');
  CAUSES.forEach((cause, i) => {
    const chip = document.createElement('label');
    chip.className = 'chip';
    const radio = document.createElement('input');
    radio.type = 'radio';
    radio.name = 'cause';
    radio.checked = i === 0;
    radio.addEventListener('change', () => { selectedCause = cause; renderAll(); });
    chip.appendChild(radio);
    chip.appendChild(Object.assign(document.createElement('span'), { textContent: cause.label }));
    causeFilterRow.appendChild(chip);
  });

  function renderAll() {
    const note = document.getElementById('chart-section-note');
    note.textContent = pd
      ? `Mean across organisms, ${pd.xValues.length} × ${pd.yValues.length} combinations (${pd.paramX} × ${pd.paramY}).`
      : `Mean across organisms, ${DATA.values.length} values of ${DATA.param} tested.`;

    const chartGrid = document.getElementById('chart-grid');
    chartGrid.innerHTML = '';
    const color = cssColor(selectedCause.color);

    DATA.cellTypes.forEach(ct => {
      if (!selectedTypes.has(ct.key)) return;
      const card = document.createElement('div');
      card.className = 'card';
      card.appendChild(Object.assign(document.createElement('h3'), { textContent: ct.label }));

      if (pd) {
        const grid = (pd.grids[ct.key] || {})[selectedCause.key];
        if (grid) {
          renderHeatmap(card, { xLabels: pd.xValues, yLabels: pd.yValues, xTitle: pd.paramX, yTitle: pd.paramY, grid });
        } else {
          card.appendChild(Object.assign(document.createElement('p'), { className: 'card-subnote', textContent: 'No data for this cause.' }));
        }
      } else {
        const values = ct.points.map(p => p[selectedCause.field]);
        renderLineChart(card, {
          xLabels: ct.points.map(p => p.value), xTitle: DATA.param,
          values, color, name: selectedCause.label,
        });
      }
      chartGrid.appendChild(card);
    });

    renderTable(document.getElementById('data-table'), DATA.cellTypes.filter(ct => selectedTypes.has(ct.key)));
  }

  renderAll();
})();
"""


def render_html(data):
    param = data["param"]
    n_orgs_values = sorted({p["n_orgs"] for ct in data["cellTypes"] for p in ct["points"]})
    n_orgs_str = str(n_orgs_values[0]) if len(n_orgs_values) == 1 else "/".join(str(n) for n in n_orgs_values)
    if data["phaseDiagram"]:
        pd = data["phaseDiagram"]
        subtitle = (
            f"{len(pd['xValues'])} × {len(pd['yValues'])} combinations tested "
            f"({pd['paramX']} × {pd['paramY']}), {n_orgs_str} organisms per combination."
        )
    else:
        subtitle = f"{len(data['values'])} values tested, {n_orgs_str} organisms per value."
    return PAGE_TEMPLATE.format(
        param=param,
        subtitle=subtitle,
        data_json=json.dumps(data).replace("</", "<\\/"),
        js_engine=JS_ENGINE,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sweep_dir", help="sweep_results/PARAM_NAME directory produced by run_param_sweep.py")
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

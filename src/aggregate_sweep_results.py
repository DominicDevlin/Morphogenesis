#!/usr/bin/env python3
"""
Aggregates the per-value death_summary.csv files produced by
run_param_sweep.sh into one table showing, for each swept value and cell
type, the average outcome across the organisms simulated at that value.

This is the "what happens on average" view: run_param_sweep.sh already
tests each value on par.n_orgs organisms, but each value's results are
spread over one death_summary.csv per value. This script averages over the
organisms within each value, and lines the values up so the trend across
the sweep is visible directly.

It also reports the cause-of-death split (sweep_results/PARAM/death_causes_summary.csv):
for each swept value and cell type, the mean number of cells killed by
lonely/blastocoel extrusion (ToxictoLonelyCells) vs. neighbour-competition
signalling (NeighbourBasedApoptosis), read from each organism's
org-data/death_causes-org-*.dat. Sweeps run before this per-cell-type
breakdown existed are skipped for this part.

The "total" pseudo cell type (the whole organism) also reports
mean_frac_eliminated: the proportion of the starting population (before
which no apoptosis mechanism runs) no longer present in the final sample -
how much cell sorting culled overall.

Works the same way for sweeps over two parameters at once (run_param_sweep.sh
PARAM1 v1,v2 PARAM2 v1,v2): each combination is just another "value", named
"PARAM1=v1,PARAM2=v2" by run_param_sweep.sh. generate_sweep_report.py detects
that naming and renders a phase-diagram heatmap instead of line charts.

Usage:
    python3 aggregate_sweep_results.py SWEEP_DIR [--out FILE]

SWEEP_DIR is the sweep_results/PARAM_NAME (or PARAM1+PARAM2) directory
created by run_param_sweep.sh (containing one subdirectory per tested
value/combination, each with its own death_summary.csv).
"""
import argparse
import csv
import glob
import os
import re
import sys

CELL_TYPES = ["zona_pellucida", "sox2_high", "sox17_high", "undifferentiated", "differentiated", "total"]
RAW_DEATH_CAUSE_TYPES = ["zona_pellucida", "sox2_high", "sox17_high", "undifferentiated"]
DEATH_CAUSE_TYPES = RAW_DEATH_CAUSE_TYPES + ["differentiated", "total"]


def sort_key(value):
    try:
        return (0, float(value))
    except ValueError:
        return (1, value)


def parse_combo_key(key):
    """Parses a sweep_results leaf directory name into an ordered list of
    (param, value) pairs, for multi-parameter sweeps (run_param_sweep.sh
    names those "param1=v1,param2=v2,..."). Returns None for a
    single-parameter sweep, where the leaf name is just the bare value."""
    if "=" not in key:
        return None
    pairs = []
    for part in key.split(","):
        name, sep, val = part.partition("=")
        if not sep:
            return None
        pairs.append((name, val))
    return pairs


def phase_diagram_axes(results):
    """If every key in `results` is a multi-parameter combo over the same two
    parameter names (in the same order), returns (param_x, param_y); else
    None. Used to decide whether a phase-diagram heatmap can be built - only
    supported for exactly two swept parameters."""
    keys = list(results.keys())
    if not keys:
        return None
    parsed = [parse_combo_key(k) for k in keys]
    if any(p is None for p in parsed):
        return None
    param_names = [name for name, _ in parsed[0]]
    if len(param_names) != 2:
        return None
    if any([name for name, _ in p] != param_names for p in parsed):
        return None
    return tuple(param_names)


def load_value_dir(path):
    """Reads death_summary.csv in `path` and returns per-cell-type stats."""
    csv_path = os.path.join(path, "death_summary.csv")
    rows_by_type = {ct: [] for ct in CELL_TYPES}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            # Ignore cell types no longer tracked (e.g. death_summary.csv files
            # from a sweep run before "loser" was retired) rather than crash.
            if row["cell_type"] in rows_by_type:
                rows_by_type[row["cell_type"]].append(row)

    stats = {}
    for cell_type, rows in rows_by_type.items():
        n_orgs = len(rows)
        if n_orgs == 0:
            continue
        mean_n_killed = sum(int(r["n_killed"]) for r in rows) / n_orgs
        extinct_rows = [r for r in rows if r["extinction_time"] != ""]
        frac_extinct = len(extinct_rows) / n_orgs
        mean_extinction_time = (
            sum(int(r["extinction_time"]) for r in extinct_rows) / len(extinct_rows)
            if extinct_rows else ""
        )
        eliminated_rows = [r for r in rows if r.get("frac_eliminated", "")]
        mean_frac_eliminated = (
            sum(float(r["frac_eliminated"]) for r in eliminated_rows) / len(eliminated_rows)
            if eliminated_rows else ""
        )
        stats[cell_type] = {
            "n_orgs": n_orgs,
            "mean_n_killed": mean_n_killed,
            "frac_extinct": frac_extinct,
            "mean_extinction_time": mean_extinction_time,
            "mean_frac_eliminated": mean_frac_eliminated,
        }
    return stats


def collect(sweep_dir):
    """Returns {value: {cell_type: stats}}, values sorted numerically when possible."""
    value_dirs = [d for d in sorted(glob.glob(os.path.join(sweep_dir, "*"))) if os.path.isdir(d)]
    if not value_dirs:
        sys.exit(f"No value subdirectories found in {sweep_dir}")

    results = {}
    for value_dir in value_dirs:
        value = os.path.basename(value_dir.rstrip("/"))
        csv_path = os.path.join(value_dir, "death_summary.csv")
        if not os.path.isfile(csv_path):
            print(f"Skipping {value_dir}: no death_summary.csv", file=sys.stderr)
            continue
        results[value] = load_value_dir(value_dir)

    return dict(sorted(results.items(), key=lambda kv: sort_key(kv[0])))


def _org_index_from_filename(path):
    m = re.search(r"-org-(\d+)\.dat$", path)
    return int(m.group(1)) if m else None


def load_total_ever_born(value_dir):
    """Reads each org's celltypes-org-*.dat and returns {org_index:
    total_ever_born} - initial_count plus every division along the way (same
    net-increase convention as n_killed uses net-decreases), i.e. everyone
    who will ever have been born by the end of the run, not just the
    starting population. {} for sweeps run before initial_count was
    tracked."""
    org_files = sorted(glob.glob(os.path.join(value_dir, "org-data", "celltypes-org-*.dat")))
    result = {}
    for path in org_files:
        org_index = _org_index_from_filename(path)
        if org_index is None:
            continue
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows or not rows[0].get("initial_count"):
            continue
        initial_count = int(rows[0]["initial_count"])
        if not initial_count:
            continue
        cumulative_born = 0
        prev_total = None
        for r in rows:
            cur_total = int(r["total"])
            if prev_total is not None and cur_total > prev_total:
                cumulative_born += cur_total - prev_total
            prev_total = cur_total
        result[org_index] = initial_count + cumulative_born
    return result


def load_final_type_counts(value_dir):
    """Reads each org's celltypes-org-*.dat and returns {org_index:
    {cell_type: final_alive_count}} from the last sampled row, for whichever
    raw types (sox2_high, sox17_high, undifferentiated, zona_pellucida) this
    model version's output includes, plus a derived "differentiated"
    (sox2_high + sox17_high) when both are present."""
    org_files = sorted(glob.glob(os.path.join(value_dir, "org-data", "celltypes-org-*.dat")))
    result = {}
    for path in org_files:
        org_index = _org_index_from_filename(path)
        if org_index is None:
            continue
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows:
            continue
        last = rows[-1]
        counts = {ct: int(last[ct]) for ct in RAW_DEATH_CAUSE_TYPES if ct in last}
        if "sox2_high" in counts and "sox17_high" in counts:
            counts["differentiated"] = counts["sox2_high"] + counts["sox17_high"]
        result[org_index] = counts
    return result


def load_death_causes(value_dir):
    """Reads each org's death_causes-org-*.dat (time series of cumulative
    lonely/signal kill counts, per cell type) and returns the final (total)
    per-cell-type counts per org, plus each count's share of that cell
    type's own cohort (frac_lonely/frac_signal_of_born: None if
    unavailable). For a specific cell type, that cohort is "everyone who was
    ever seen as this type" - cells still alive as this type at the end of
    the run, plus cells that died as this type along the way (killed_of_type
    / (final_alive_of_type + killed_of_type)) - not the whole organism's
    total_ever_born, so e.g. undifferentiated (which nearly always empties
    out as cells commit to a lineage or die) can legitimately reach 100%.
    "total" (the whole organism) is the one exception: its cohort is
    everyone ever born (initial population plus every division), since
    there's no "different type" for the whole organism to have transitioned
    into. Returns [] for sweeps run before the per-cell-type cause breakdown
    existed (old files only had a single overall lonely/signal total, no
    per-type columns)."""
    total_ever_born_by_org = load_total_ever_born(value_dir)
    final_counts_by_org = load_final_type_counts(value_dir)
    org_files = sorted(glob.glob(os.path.join(value_dir, "org-data", "death_causes-org-*.dat")))
    per_org = []
    for path in org_files:
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows or "sox2_high_lonely" not in rows[-1]:
            continue  # empty, or pre-per-cell-type format
        last = rows[-1]
        org_index = _org_index_from_filename(path)
        total_ever_born = total_ever_born_by_org.get(org_index)
        final_counts = final_counts_by_org.get(org_index, {})

        # Only build an entry for cell types whose columns actually exist in
        # this file - e.g. zona_pellucida has been dropped from some model
        # versions' output (it never dies here), and old files might lack it
        # while newer ones do or vice versa.
        counts = {}
        for cell_type in RAW_DEATH_CAUSE_TYPES:
            if f"{cell_type}_lonely" not in last:
                continue
            lonely = int(last[f"{cell_type}_lonely"])
            signal = int(last[f"{cell_type}_signal"])
            final_alive = final_counts.get(cell_type)
            type_cohort = (final_alive + lonely + signal) if final_alive is not None else None
            counts[cell_type] = {
                "lonely": lonely,
                "signal": signal,
                "frac_lonely_of_born": (lonely / type_cohort) if type_cohort else None,
                "frac_signal_of_born": (signal / type_cohort) if type_cohort else None,
            }
        if "sox2_high" in counts and "sox17_high" in counts:
            d2, d17 = counts["sox2_high"], counts["sox17_high"]
            diff_lonely, diff_signal = d2["lonely"] + d17["lonely"], d2["signal"] + d17["signal"]
            diff_final = final_counts.get("differentiated")
            diff_cohort = (diff_final + diff_lonely + diff_signal) if diff_final is not None else None
            counts["differentiated"] = {
                "lonely": diff_lonely,
                "signal": diff_signal,
                "frac_lonely_of_born": (diff_lonely / diff_cohort) if diff_cohort else None,
                "frac_signal_of_born": (diff_signal / diff_cohort) if diff_cohort else None,
            }
        # "total" reads the file's own whole-organism cumulative columns
        # directly, rather than summing per-type ones - always present
        # regardless of which per-type columns this model version writes -
        # and normalizes against total_ever_born (see docstring: the whole
        # organism has no "other type" to have transitioned into).
        if "total_lonely" in last:
            total_lonely = int(last["total_lonely"])
            total_signal = int(last["total_signal"])
            counts["total"] = {
                "lonely": total_lonely,
                "signal": total_signal,
                "frac_lonely_of_born": (total_lonely / total_ever_born) if total_ever_born else None,
                "frac_signal_of_born": (total_signal / total_ever_born) if total_ever_born else None,
            }
        per_org.append(counts)
    return per_org


def collect_death_causes(sweep_dir):
    """Returns {value: {cell_type: {n_orgs, mean_lonely, mean_signal,
    mean_frac_lonely_of_born, mean_frac_signal_of_born}}}; skips values with
    no (or pre-per-cell-type) death_causes-org-*.dat. The mean_frac_*_of_born
    fields average, across organisms, each cause's share of that organism's
    total_ever_born - "" if initial_count wasn't tracked for that sweep."""
    value_dirs = [d for d in sorted(glob.glob(os.path.join(sweep_dir, "*"))) if os.path.isdir(d)]
    results = {}
    for value_dir in value_dirs:
        value = os.path.basename(value_dir.rstrip("/"))
        per_org = load_death_causes(value_dir)
        if not per_org:
            continue
        by_type = {}
        for cell_type in DEATH_CAUSE_TYPES:
            # Not every org necessarily has this cell type's columns (e.g.
            # zona_pellucida dropped from some model versions' output).
            present = [o[cell_type] for o in per_org if cell_type in o]
            if not present:
                continue
            m = len(present)
            lonely_fracs = [o["frac_lonely_of_born"] for o in present if o["frac_lonely_of_born"] is not None]
            signal_fracs = [o["frac_signal_of_born"] for o in present if o["frac_signal_of_born"] is not None]
            by_type[cell_type] = {
                "n_orgs": m,
                "mean_lonely": sum(o["lonely"] for o in present) / m,
                "mean_signal": sum(o["signal"] for o in present) / m,
                "mean_frac_lonely_of_born": (sum(lonely_fracs) / len(lonely_fracs)) if lonely_fracs else "",
                "mean_frac_signal_of_born": (sum(signal_fracs) / len(signal_fracs)) if signal_fracs else "",
            }
        if not by_type:
            continue
        results[value] = by_type
    return dict(sorted(results.items(), key=lambda kv: sort_key(kv[0])))


def load_elimination_timeseries(value_dir):
    """Reads each org's celltypes-org-*.dat and returns, for each sampled
    time common to every organism, the mean cumulative fraction of cells
    killed so far relative to *everyone who will ever have been born by the
    end of the run* (initial_count plus every division along the way, not
    just the starting population) - averaged across organisms. [] for
    sweeps run before initial_count was tracked.

    Cumulative kills and cumulative births both use the same net-change
    convention as analyse_differentiated_deaths.py's n_killed: a drop in
    total population between two consecutive samples counts as that many
    deaths, a rise counts as that many births. The denominator
    (initial_count + total births) is fixed at the run's final value, so
    the reported curve is monotonically increasing - killed so far divided
    by the total cohort size once the run ends, not by the (shrinking or
    growing) population at time t."""
    total_ever_born_by_org = load_total_ever_born(value_dir)
    org_files = sorted(glob.glob(os.path.join(value_dir, "org-data", "celltypes-org-*.dat")))
    per_org_series = []
    for path in org_files:
        total_ever_born = total_ever_born_by_org.get(_org_index_from_filename(path))
        if not total_ever_born:
            continue
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

        killed_by_time = {}
        cumulative_killed = 0
        prev_total = None
        for r in rows:
            cur_total = int(r["total"])
            if prev_total is not None and cur_total < prev_total:
                cumulative_killed += prev_total - cur_total
            prev_total = cur_total
            killed_by_time[int(r["time"])] = cumulative_killed

        per_org_series.append({t: k / total_ever_born for t, k in killed_by_time.items()})

    if not per_org_series:
        return []
    common_times = sorted(set.intersection(*(set(s.keys()) for s in per_org_series)))
    return [
        {"time": t, "mean_frac_killed": sum(s[t] for s in per_org_series) / len(per_org_series)}
        for t in common_times
    ]


def collect_elimination_timeseries(sweep_dir):
    """Returns {value: [{time, mean_frac_killed}, ...]}; skips values with no
    initial_count tracking (sweeps run before it existed)."""
    value_dirs = [d for d in sorted(glob.glob(os.path.join(sweep_dir, "*"))) if os.path.isdir(d)]
    results = {}
    for value_dir in value_dirs:
        value = os.path.basename(value_dir.rstrip("/"))
        series = load_elimination_timeseries(value_dir)
        if series:
            results[value] = series
    return dict(sorted(results.items(), key=lambda kv: sort_key(kv[0])))


def write_elimination_timeseries_csv(out_path, results):
    fieldnames = ["value", "time", "mean_frac_killed"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for value, series in results.items():
            for point in series:
                writer.writerow({"value": value, "time": point["time"], "mean_frac_killed": point["mean_frac_killed"]})
    print(f"Cumulative-proportion-killed-over-time series written to {out_path}")


def write_csv(out_path, param, results):
    fieldnames = ["value", "cell_type", "n_orgs", "mean_n_killed", "frac_extinct", "mean_extinction_time", "mean_frac_eliminated"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for value, by_type in results.items():
            for cell_type, stats in by_type.items():
                writer.writerow({"value": value, "cell_type": cell_type, **stats})
    print(f"Summary written to {out_path}")


def write_death_causes_csv(out_path, results):
    fieldnames = [
        "value", "cell_type", "n_orgs", "mean_lonely_killed", "mean_signal_killed", "frac_signal",
        "mean_frac_lonely_of_born", "mean_frac_signal_of_born",
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for value, by_type in results.items():
            for cell_type, stats in by_type.items():
                total = stats["mean_lonely"] + stats["mean_signal"]
                frac_signal = stats["mean_signal"] / total if total else ""
                writer.writerow({
                    "value": value,
                    "cell_type": cell_type,
                    "n_orgs": stats["n_orgs"],
                    "mean_lonely_killed": stats["mean_lonely"],
                    "mean_signal_killed": stats["mean_signal"],
                    "frac_signal": frac_signal,
                    "mean_frac_lonely_of_born": stats["mean_frac_lonely_of_born"],
                    "mean_frac_signal_of_born": stats["mean_frac_signal_of_born"],
                })
    print(f"Death-cause summary written to {out_path}")


def print_death_causes_table(param, results):
    if not results:
        return
    print(f"\n=== Cause of death vs {param} (lonely/blastocoel vs signal-based apoptosis) ===")
    for cell_type in DEATH_CAUSE_TYPES:
        if not any(cell_type in by_type for by_type in results.values()):
            continue
        print(f"\n  {cell_type}:")
        print(
            f"  {'value':>12}  {'n_orgs':>6}  {'mean_lonely':>11}  {'mean_signal':>11}  {'frac_signal':>11}  "
            f"{'lonely/born':>11}  {'signal/born':>11}"
        )
        for value, by_type in results.items():
            stats = by_type.get(cell_type)
            if stats is None:
                continue
            total = stats["mean_lonely"] + stats["mean_signal"]
            frac_signal = f"{stats['mean_signal'] / total:.2f}" if total else "-"
            lonely_of_born = f"{stats['mean_frac_lonely_of_born']:.2f}" if stats["mean_frac_lonely_of_born"] != "" else "-"
            signal_of_born = f"{stats['mean_frac_signal_of_born']:.2f}" if stats["mean_frac_signal_of_born"] != "" else "-"
            print(
                f"  {value:>12}  {stats['n_orgs']:>6}  {stats['mean_lonely']:>11.2f}  "
                f"{stats['mean_signal']:>11.2f}  {frac_signal:>11}  {lonely_of_born:>11}  {signal_of_born:>11}"
            )


def print_table(param, results):
    print(f"\n=== Average outcome vs {param} ===")
    for cell_type in CELL_TYPES:
        if not any(cell_type in by_type for by_type in results.values()):
            continue
        print(f"\n  {cell_type}:")
        header = f"  {'value':>12}  {'n_orgs':>6}  {'mean_killed':>11}  {'frac_extinct':>12}  {'mean_ext_t':>10}"
        if cell_type == "total":
            header += f"  {'frac_eliminated':>15}"
        print(header)
        for value, by_type in results.items():
            stats = by_type.get(cell_type)
            if stats is None:
                continue
            ext_t = f"{stats['mean_extinction_time']:.1f}" if stats["mean_extinction_time"] != "" else "-"
            line = (
                f"  {value:>12}  {stats['n_orgs']:>6}  {stats['mean_n_killed']:>11.2f}  "
                f"{stats['frac_extinct']:>12.2f}  {ext_t:>10}"
            )
            if cell_type == "total":
                frac_elim = stats.get("mean_frac_eliminated", "")
                line += f"  {frac_elim:>15.2f}" if frac_elim != "" else f"  {'-':>15}"
            print(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("sweep_dir", help="sweep_results/PARAM_NAME directory produced by run_param_sweep.sh")
    parser.add_argument("--out", default=None, help="CSV file to write (default: SWEEP_DIR/summary.csv)")
    args = parser.parse_args()

    sweep_dir = args.sweep_dir.rstrip("/")
    param = os.path.basename(sweep_dir)
    out_path = args.out or os.path.join(sweep_dir, "summary.csv")

    results = collect(sweep_dir)
    print_table(param, results)
    write_csv(out_path, param, results)

    death_cause_results = collect_death_causes(sweep_dir)
    print_death_causes_table(param, death_cause_results)
    if death_cause_results:
        death_cause_out = os.path.join(os.path.dirname(out_path), "death_causes_summary.csv")
        write_death_causes_csv(death_cause_out, death_cause_results)

    elimination_results = collect_elimination_timeseries(sweep_dir)
    if elimination_results:
        elimination_out = os.path.join(os.path.dirname(out_path), "elimination_timeseries.csv")
        write_elimination_timeseries_csv(elimination_out, elimination_results)


if __name__ == "__main__":
    main()

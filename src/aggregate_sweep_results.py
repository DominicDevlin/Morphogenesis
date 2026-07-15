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
for each swept value, the mean number of cells killed by lonely/blastocoel
extrusion (ToxictoLonelyCells) vs. neighbour-competition signalling
(NeighbourBasedApoptosis), read from each organism's
org-data/death_causes-org-*.dat. Sweeps run before this data existed are
skipped for this part.

Usage:
    python3 aggregate_sweep_results.py SWEEP_DIR [--out FILE]

SWEEP_DIR is the sweep_results/PARAM_NAME directory created by
run_param_sweep.sh (containing one subdirectory per tested value, each with
its own death_summary.csv).
"""
import argparse
import csv
import glob
import os
import sys

CELL_TYPES = ["zona_pellucida", "sox2_high", "sox17_high", "loser", "undifferentiated", "differentiated"]


def sort_key(value):
    try:
        return (0, float(value))
    except ValueError:
        return (1, value)


def load_value_dir(path):
    """Reads death_summary.csv in `path` and returns per-cell-type stats."""
    csv_path = os.path.join(path, "death_summary.csv")
    rows_by_type = {ct: [] for ct in CELL_TYPES}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
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
        stats[cell_type] = {
            "n_orgs": n_orgs,
            "mean_n_killed": mean_n_killed,
            "frac_extinct": frac_extinct,
            "mean_extinction_time": mean_extinction_time,
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


def load_death_causes(value_dir):
    """Reads each org's death_causes-org-*.dat (time series of cumulative
    lonely/signal kill counts) and returns the final (total) counts per org."""
    org_files = sorted(glob.glob(os.path.join(value_dir, "org-data", "death_causes-org-*.dat")))
    per_org = []
    for path in org_files:
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows:
            continue
        last = rows[-1]
        per_org.append({"lonely": int(last["lonely_killed"]), "signal": int(last["signal_killed"])})
    return per_org


def collect_death_causes(sweep_dir):
    """Returns {value: {n_orgs, mean_lonely, mean_signal}}; skips values with
    no death_causes-org-*.dat (e.g. sweeps run before this data existed)."""
    value_dirs = [d for d in sorted(glob.glob(os.path.join(sweep_dir, "*"))) if os.path.isdir(d)]
    results = {}
    for value_dir in value_dirs:
        value = os.path.basename(value_dir.rstrip("/"))
        per_org = load_death_causes(value_dir)
        if not per_org:
            continue
        n_orgs = len(per_org)
        results[value] = {
            "n_orgs": n_orgs,
            "mean_lonely": sum(o["lonely"] for o in per_org) / n_orgs,
            "mean_signal": sum(o["signal"] for o in per_org) / n_orgs,
        }
    return dict(sorted(results.items(), key=lambda kv: sort_key(kv[0])))


def write_csv(out_path, param, results):
    fieldnames = ["value", "cell_type", "n_orgs", "mean_n_killed", "frac_extinct", "mean_extinction_time"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for value, by_type in results.items():
            for cell_type, stats in by_type.items():
                writer.writerow({"value": value, "cell_type": cell_type, **stats})
    print(f"Summary written to {out_path}")


def write_death_causes_csv(out_path, results):
    fieldnames = ["value", "n_orgs", "mean_lonely_killed", "mean_signal_killed", "frac_signal"]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for value, stats in results.items():
            total = stats["mean_lonely"] + stats["mean_signal"]
            frac_signal = stats["mean_signal"] / total if total else ""
            writer.writerow({
                "value": value,
                "n_orgs": stats["n_orgs"],
                "mean_lonely_killed": stats["mean_lonely"],
                "mean_signal_killed": stats["mean_signal"],
                "frac_signal": frac_signal,
            })
    print(f"Death-cause summary written to {out_path}")


def print_death_causes_table(param, results):
    if not results:
        return
    print(f"\n=== Cause of death vs {param} (lonely/blastocoel vs signal-based apoptosis) ===")
    print(f"  {'value':>12}  {'n_orgs':>6}  {'mean_lonely':>11}  {'mean_signal':>11}  {'frac_signal':>11}")
    for value, stats in results.items():
        total = stats["mean_lonely"] + stats["mean_signal"]
        frac_signal = f"{stats['mean_signal'] / total:.2f}" if total else "-"
        print(
            f"  {value:>12}  {stats['n_orgs']:>6}  {stats['mean_lonely']:>11.2f}  "
            f"{stats['mean_signal']:>11.2f}  {frac_signal:>11}"
        )


def print_table(param, results):
    print(f"\n=== Average outcome vs {param} ===")
    for cell_type in CELL_TYPES:
        if not any(cell_type in by_type for by_type in results.values()):
            continue
        print(f"\n  {cell_type}:")
        print(f"  {'value':>12}  {'n_orgs':>6}  {'mean_killed':>11}  {'frac_extinct':>12}  {'mean_ext_t':>10}")
        for value, by_type in results.items():
            stats = by_type.get(cell_type)
            if stats is None:
                continue
            ext_t = f"{stats['mean_extinction_time']:.1f}" if stats["mean_extinction_time"] != "" else "-"
            print(
                f"  {value:>12}  {stats['n_orgs']:>6}  {stats['mean_n_killed']:>11.2f}  "
                f"{stats['frac_extinct']:>12.2f}  {ext_t:>10}"
            )


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


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Analyses the celltypes-org-*.dat files produced by embryo_multi (columns:
time, zona_pellucida, sox2_high, sox17_high, loser, undifferentiated,
total) to find, for every cell type, when cells of that type are killed
and when they have all died out (if that happens).

A "death" is inferred from a decrease in a type's cell count between two
consecutive samples. "differentiated" is a derived type combining
sox2_high and sox17_high.

Results are printed to the console and also logged to a CSV file (one row
per organism/cell type).

Usage:
    python3 analyse_differentiated_deaths.py [path] [--csv FILE]

`path` can be a .dat file, or a directory containing celltypes-org-*.dat
files (defaults to org-data, next to this script).
"""
import argparse
import csv
import glob
import os
import sys

RAW_TYPES = ["zona_pellucida", "sox2_high", "sox17_high", "loser", "undifferentiated"]
CELL_TYPES = RAW_TYPES + ["differentiated"]


def load_counts(path):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            counts = {col: int(row[col]) for col in RAW_TYPES}
            counts["differentiated"] = counts["sox2_high"] + counts["sox17_high"]
            rows.append({"time": int(row["time"]), "counts": counts})
    return rows


def find_death_events(rows, cell_type):
    """List of (time, n_killed) where a type's cell count dropped from one
    sample to the next."""
    events = []
    prev = None
    for row in rows:
        cur = row["counts"][cell_type]
        if prev is not None and cur < prev:
            events.append((row["time"], prev - cur))
        prev = cur
    return events


def find_extinction_time(rows, cell_type):
    """First time a type's cell count falls back to 0 after having been
    positive. None if that never happens."""
    seen_positive = False
    for row in rows:
        cur = row["counts"][cell_type]
        if cur > 0:
            seen_positive = True
        elif seen_positive and cur == 0:
            return row["time"]
    return None


def summarise_type(path_rows, cell_type):
    events = find_death_events(path_rows, cell_type)
    extinction_time = find_extinction_time(path_rows, cell_type)
    return events, extinction_time


def format_summary(cell_type, rows, events, extinction_time):
    label = f"{cell_type:>16}"
    if not events:
        return f"  {label}: no deaths"

    n_killed = sum(n for _, n in events)
    line = f"  {label}: {n_killed} killed (first t={events[0][0]}, last t={events[-1][0]})"

    if extinction_time is not None:
        still_extinct = rows[-1]["counts"][cell_type] == 0
        note = "" if still_extinct else ", recovered afterwards"
        line += f", extinct at t={extinction_time}{note}"

    return line


def report(path, rows):
    """Prints the console summary and returns the matching CSV rows."""
    org = os.path.basename(path)
    print(f"== {org} ==")

    csv_rows = []
    for cell_type in CELL_TYPES:
        events, extinction_time = summarise_type(rows, cell_type)
        print(format_summary(cell_type, rows, events, extinction_time))

        recovered = extinction_time is not None and rows[-1]["counts"][cell_type] != 0
        csv_rows.append({
            "org": org,
            "cell_type": cell_type,
            "n_killed": sum(n for _, n in events),
            "first_death_time": events[0][0] if events else "",
            "last_death_time": events[-1][0] if events else "",
            "extinction_time": extinction_time if extinction_time is not None else "",
            "recovered_after_extinction": recovered,
        })
    print()
    return csv_rows


def write_csv(csv_path, csv_rows):
    fieldnames = [
        "org", "cell_type", "n_killed",
        "first_death_time", "last_death_time",
        "extinction_time", "recovered_after_extinction",
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"Results written to {csv_path}")


def resolve_files(path):
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, "celltypes-org-*.dat")))
        if not files:
            sys.exit(f"No celltypes-org-*.dat file found in {path}")
        return files
    if os.path.isfile(path):
        return [path]
    sys.exit(f"Path not found: {path}")


def main():
    default_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "org-data")

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "path",
        nargs="?",
        default=default_path,
        help=f"celltypes-org-*.dat file, or directory containing them (default: {default_path})",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="ignore samples before this time, e.g. to skip the pre-differentiation burn-in (default: 0)",
    )
    parser.add_argument(
        "--csv",
        default="death_summary.csv",
        help="CSV file to write the results to (default: death_summary.csv)",
    )
    args = parser.parse_args()

    csv_rows = []
    for path in resolve_files(args.path):
        rows = [row for row in load_counts(path) if row["time"] >= args.start]
        csv_rows.extend(report(path, rows))

    write_csv(args.csv, csv_rows)


if __name__ == "__main__":
    main()

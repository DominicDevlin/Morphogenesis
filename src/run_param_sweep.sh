#!/usr/bin/env bash
#
# Runs embryo_multi several times in a row, each time with a different value
# of one parameter from parameter-files/parameter_embryo_multi.cpp, so the
# resulting runs can be compared.
#
# Usage:
#   ./run_param_sweep.sh PARAM_NAME VALUE1 [VALUE2 ...]
#
# Examples:
#   ./run_param_sweep.sh starting_fraction_losers 0.1 0.25 0.4
#   ./run_param_sweep.sh set_loser_colours true false
#   ./run_param_sweep.sh T 0.8 1 1.2
#
# For each VALUE the script:
#   1. patches PARAM_NAME's initial value in parameter_embryo_multi.cpp
#      (only the first assignment line - any later line that derives another
#      parameter from PARAM_NAME, e.g. motility_zero from motility_strength,
#      picks up the new value automatically since it's computed afterwards)
#   2. rebuilds embryo_multi (incremental `make`)
#   3. runs it (par.n_orgs organisms per value, in parallel)
#   4. moves the resulting org-data/ and photos/ into
#      sweep_results/PARAM_NAME/VALUE/
#   5. runs analyse_differentiated_deaths.py on that run's org-data,
#      writing sweep_results/PARAM_NAME/VALUE/death_summary.csv
#
# Once every value has been run, aggregate_sweep_results.py averages each
# value's per-organism results (death counts, extinction times) into
# sweep_results/PARAM_NAME/summary.csv, so the effect of PARAM_NAME "on
# average" is visible directly instead of having to compare per-value CSVs
# by hand. generate_sweep_report.py then turns that into a single browsable
# sweep_results/PARAM_NAME/report.html (charts + full data tables).
#
# parameter_embryo_multi.cpp is restored to its original content when the
# script exits (even on error/Ctrl-C), and embryo_multi is rebuilt one last
# time so the binary on disk matches the checked-in source again.
#
# NOTE: values are substituted verbatim into the C++ source, so pass them in
# valid C++ literal syntax: booleans as true/false, numbers as e.g. 0.35 or
# -0.7, strings quoted e.g. '"foo"'.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PARAM_FILE="parameter-files/parameter_embryo_multi.cpp"
PRO_FILE="CellularPotts2.pro"
BINARY="./embryo_multi"
ANALYSE_SCRIPT="./analyse_differentiated_deaths.py"
AGGREGATE_SCRIPT="./aggregate_sweep_results.py"
REPORT_SCRIPT="./generate_sweep_report.py"

qmake_regen() {
  local qmake_log
  qmake_log="$(mktemp)"
  if command -v qmake-qt5 >/dev/null 2>&1; then
    qmake-qt5 > "$qmake_log" 2>&1
  else
    qmake > "$qmake_log" 2>&1
  fi || { cat "$qmake_log"; rm -f "$qmake_log"; echo "Error: qmake failed." >&2; exit 1; }
  rm -f "$qmake_log"
}

if [ $# -lt 2 ]; then
  echo "Usage: $0 PARAM_NAME VALUE1 [VALUE2 ...]" >&2
  exit 1
fi

PARAM="$1"
shift
VALUES=("$@")

case "$PARAM" in
  data_file|pic_dir)
    echo "Error: '$PARAM' is managed automatically by this script (it names the output dirs), sweep a different parameter." >&2
    exit 1
    ;;
esac

LINE_NO="$(grep -n -m1 -E "^[[:space:]]*${PARAM}[[:space:]]*=" "$PARAM_FILE" | cut -d: -f1)"
if [ -z "$LINE_NO" ]; then
  echo "Error: no assignment to '$PARAM' found in $PARAM_FILE" >&2
  exit 1
fi

if [ -e org-data ] || [ -e photos ]; then
  echo "Error: $SCRIPT_DIR/org-data and/or photos already exist." >&2
  echo "They look like results from a previous run and would be overwritten. Move or remove them first, then re-run." >&2
  exit 1
fi

BACKUP="$(mktemp)"
cp "$PARAM_FILE" "$BACKUP"
PRO_BACKUP="$(mktemp)"
cp "$PRO_FILE" "$PRO_BACKUP"
restore() {
  cp "$BACKUP" "$PARAM_FILE"
  rm -f "$BACKUP"
  cp "$PRO_BACKUP" "$PRO_FILE"
  rm -f "$PRO_BACKUP"
  echo "Restored $PARAM_FILE and $PRO_FILE, rebuilding original binary..."
  qmake_regen
  RESTORE_BUILD_LOG="$(mktemp)"
  if ! make -j"$(nproc)" > "$RESTORE_BUILD_LOG" 2>&1; then
    cat "$RESTORE_BUILD_LOG"
    echo "Error: rebuilding the original binary failed." >&2
  fi
  rm -f "$RESTORE_BUILD_LOG"
}
trap restore EXIT

# CellularPotts2.pro's TARGET decides which binary `make` produces (and which
# parameter file it compiles in) - force it to embryo_multi regardless of
# what's currently checked in, since this script always needs the
# multi-organism binary.
sed -i "s/^TARGET[[:space:]]*=.*/TARGET = embryo_multi/" "$PRO_FILE"
echo "Regenerating Makefile for embryo_multi..."
qmake_regen

RESULTS_DIR="sweep_results/${PARAM}"
mkdir -p "$RESULTS_DIR"

RUN_PREFIX=""
if [ -z "${DISPLAY:-}" ] && command -v xvfb-run >/dev/null 2>&1; then
  RUN_PREFIX="xvfb-run -a"
fi

for VALUE in "${VALUES[@]}"; do
  echo "=== $PARAM = $VALUE ==="

  cp "$BACKUP" "$PARAM_FILE"
  sed -i "${LINE_NO}s#\(^[[:space:]]*${PARAM}[[:space:]]*=[[:space:]]*\)[^;]*\(;.*\)#\1${VALUE}\2#" "$PARAM_FILE"

  echo "Building..."
  BUILD_LOG="$(mktemp)"
  if ! make -j"$(nproc)" > "$BUILD_LOG" 2>&1; then
    cat "$BUILD_LOG"
    rm -f "$BUILD_LOG"
    echo "Error: build failed for $PARAM = $VALUE" >&2
    exit 1
  fi
  rm -f "$BUILD_LOG"

  RUN_DIR="$RESULTS_DIR/$VALUE"
  rm -rf "$RUN_DIR"
  mkdir -p "$RUN_DIR"

  echo "Running..."
  $RUN_PREFIX "$BINARY" 2>&1 | tee "$RUN_DIR/run.log"

  mv org-data "$RUN_DIR/org-data"
  [ -d photos ] && mv photos "$RUN_DIR/photos"

  python3 "$ANALYSE_SCRIPT" --csv "$RUN_DIR/death_summary.csv" --start "800" "$RUN_DIR/org-data"

  echo "Done: results in $RUN_DIR"
done

python3 "$AGGREGATE_SCRIPT" "$RESULTS_DIR"
python3 "$REPORT_SCRIPT" "$RESULTS_DIR"

echo "Sweep complete. Results under $RESULTS_DIR/"

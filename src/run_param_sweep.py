#!/usr/bin/env python3
"""
Runs embryo_multi several times in a row, each time with a different value
(or combination of values, for 2+ parameters swept at once) from
parameter-files/parameter_embryo_multi.cpp, so the resulting runs can be
compared - including as a phase diagram when exactly two parameters are
swept together.

Usage:
    ./run_param_sweep.py PARAM1 v1,v2,... [PARAM2 v1,v2,... ...]

Examples:
    ./run_param_sweep.py starting_fraction_losers 0.1,0.25,0.4
    ./run_param_sweep.py set_loser_colours true,false
    ./run_param_sweep.py apop_threshold 0,5,10,20 loser_sox2_adhesion -0.3,0,0.3

With a single PARAM, results land exactly as before:
    sweep_results/PARAM_NAME/VALUE/
With two or more PARAMs, every combination (full cross product) is run, and
results land in a joined directory named after every parameter and value:
    sweep_results/PARAM1+PARAM2/PARAM1=v1,PARAM2=v2/

For each combination the script:
  1. patches every PARAM's initial value in parameter_embryo_multi.cpp
     (only the first assignment line for each - any later line that derives
     another parameter from it, e.g. motility_zero from motility_strength,
     picks up the new value automatically since it's computed afterwards)
  2. rebuilds embryo_multi (incremental `make`) and keeps a private copy of
     the resulting binary, since the next combination's build would
     otherwise overwrite it before it gets to run
  3. runs it (par.n_orgs organisms per combination, in parallel via OpenMP
     within the process)
  4. runs analyse_differentiated_deaths.py on that run's org-data, writing
     death_summary.csv alongside it

par.n_orgs is hardcoded in embryo_multi.cpp, so when it's low enough that
one process leaves cores idle, several combinations are instead built and
run concurrently (one process per combination, each with its own private
binary and working directory) so the available cores stay busy either way -
combinations are grouped into chunks of that size, built one at a time
(the build shares one object directory/binary path so can't itself run
concurrently) and then run together.

Once every combination has been run, aggregate_sweep_results.py averages
each combination's per-organism results (death counts, extinction times,
proportion of the starting population eliminated) into summary.csv, so the
effect of the swept parameter(s) "on average" is visible directly instead of
having to compare per-combination CSVs by hand. generate_sweep_report.py
then turns that into a single browsable report.html (charts + full data
tables - a phase-diagram heatmap per cell type when exactly two parameters
were swept).

parameter_embryo_multi.cpp is restored to its original content when the
script exits (even on error/Ctrl-C), and embryo_multi is rebuilt one last
time so the binary on disk matches the checked-in source again.

NOTE: values are substituted verbatim into the C++ source, so pass them in
valid C++ literal syntax: booleans as true/false, numbers as e.g. 0.35 or
-0.7, strings quoted e.g. '"foo"'. Values for one parameter are comma-
separated (no spaces) since spaces separate one parameter's group from the
next.
"""
import itertools
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PARAM_FILE = SCRIPT_DIR / "parameter-files/parameter_embryo_multi.cpp"
PRO_FILE = SCRIPT_DIR / "CellularPotts2.pro"
BINARY = SCRIPT_DIR / "embryo_multi"
ANALYSE_SCRIPT = SCRIPT_DIR / "analyse_differentiated_deaths.py"
AGGREGATE_SCRIPT = SCRIPT_DIR / "aggregate_sweep_results.py"
REPORT_SCRIPT = SCRIPT_DIR / "generate_sweep_report.py"

RESERVED_PARAMS = {"data_file", "pic_dir"}


def die(message):
    print(f"Error: {message}", file=sys.stderr)
    sys.exit(1)


def qmake_regen():
    qmake = shutil.which("qmake-qt5") or shutil.which("qmake")
    if qmake is None:
        die("no qmake or qmake-qt5 found on PATH.")
    result = subprocess.run([qmake], cwd=SCRIPT_DIR, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout + result.stderr)
        die("qmake failed.")


def make(context):
    result = subprocess.run(
        ["make", f"-j{os.cpu_count() or 1}"], cwd=SCRIPT_DIR, capture_output=True, text=True
    )
    if result.returncode != 0:
        print(result.stdout + result.stderr)
        die(f"build failed{context}")


def restore(param_backup, pro_backup):
    """Restores PARAM_FILE/PRO_FILE from their backups and rebuilds the
    original binary - best-effort, since this runs during cleanup and must
    not itself raise or mask whatever error is already propagating."""
    shutil.copy(param_backup, PARAM_FILE)
    os.unlink(param_backup)
    shutil.copy(pro_backup, PRO_FILE)
    os.unlink(pro_backup)
    print(f"Restored {PARAM_FILE} and {PRO_FILE}, rebuilding original binary...")
    try:
        qmake_regen()
        make(" while rebuilding the original binary")
    except SystemExit:
        print("Error: rebuilding the original binary failed.", file=sys.stderr)


def find_assignment_line(param):
    """1-based line number of param's first top-level assignment in
    PARAM_FILE, e.g. matching '  apop_threshold = 10;'."""
    pattern = re.compile(rf"^\s*{re.escape(param)}\s*=")
    with open(PARAM_FILE) as f:
        for line_no, line in enumerate(f, start=1):
            if pattern.match(line):
                return line_no
    return None


def patch_param_line(lines, line_no, param, value):
    """Replaces the assigned value on lines[line_no - 1], keeping everything
    up to and including '=' and the trailing ';...' (e.g. a comment)
    untouched - mirrors run_param_sweep.sh's sed substitution."""
    pattern = re.compile(rf"^(\s*{re.escape(param)}\s*=\s*)[^;]*(;.*)$")
    line = lines[line_no - 1]
    new_line, n = pattern.subn(lambda m: m.group(1) + value + m.group(2), line)
    if n == 0:
        die(f"could not patch line {line_no} for '{param}' - unexpected format: {line!r}")
    lines[line_no - 1] = new_line


def run_python(script, args):
    result = subprocess.run([sys.executable, str(script), *args])
    if result.returncode != 0:
        die(f"{script.name} exited with status {result.returncode}")


def read_n_orgs():
    """Reads the hardcoded 'par.n_orgs = N;' assignment in embryo_multi.cpp's
    main() - it isn't exposed via the parameter file. Used to tell how many
    organisms one process already simulates concurrently via OpenMP (see
    process_population()'s omp_set_num_threads call), and therefore how many
    cores are left idle for building/running further combinations at once.
    None if the assignment can't be found (e.g. the source changed shape) -
    callers then fall back to one combination at a time, same as before this
    existed."""
    match = re.search(r"^\s*par\.n_orgs\s*=\s*(\d+)\s*;", (SCRIPT_DIR / "embryo_multi.cpp").read_text(), re.MULTILINE)
    return int(match.group(1)) if match else None


def combo_leaf_and_label(param_names, combo):
    if len(param_names) == 1:
        return combo[0], f"{param_names[0]} = {combo[0]}"
    leaf_name = ",".join(f"{name}={value}" for name, value in zip(param_names, combo))
    label = ",".join(f"{name} = {value}" for name, value in zip(param_names, combo))
    return leaf_name, label


def build_combo(param_names, line_nos, param_backup, combo, label, run_dir):
    """Patches PARAM_FILE with this combination's values, rebuilds, and
    leaves a private copy of the resulting binary at
    run_dir/embryo_multi_bin - the next combination in the same batch also
    needs to rebuild before this one runs, which would otherwise overwrite
    the shared BINARY path first."""
    lines = Path(param_backup).read_text().splitlines(keepends=True)
    for line_no, param, value in zip(line_nos, param_names, combo):
        patch_param_line(lines, line_no, param, value)
    PARAM_FILE.write_text("".join(lines))

    print("Building...")
    make(f" for {label}")

    combo_binary = run_dir / "embryo_multi_bin"
    shutil.copy(BINARY, combo_binary)
    return combo_binary


def run_combo_group(run_prefix, group):
    """group: [(label, run_dir, binary), ...]. Launches each combination's
    private binary concurrently, each with cwd=its own run_dir (so its
    relative org-data/photos outputs land directly there), waits for all,
    and dies if any failed. A single-entry group behaves like one plain
    run."""
    if len(group) > 1:
        print(f"Running {len(group)} combinations concurrently: {', '.join(label for label, _, _ in group)}")
    else:
        print("Running...")

    procs = []
    for label, run_dir, binary in group:
        cmd = run_prefix + [str(binary)]
        log_file = open(run_dir / "run.log", "w")
        procs.append((label, binary, subprocess.Popen(cmd, cwd=run_dir, stdout=log_file, stderr=subprocess.STDOUT), log_file))

    failed = []
    for label, binary, proc, log_file in procs:
        proc.wait()
        log_file.close()
        binary.unlink(missing_ok=True)
        if len(group) > 1:
            print(f"  {label}: {'done' if proc.returncode == 0 else 'FAILED'}")
        if proc.returncode != 0:
            failed.append(label)
    if failed:
        die(f"embryo_multi failed for: {', '.join(failed)} - see run.log inside each combination's results directory")


def parse_args(argv):
    if len(argv) < 2 or len(argv) % 2 != 0:
        print(f"Usage: {sys.argv[0]} PARAM1 v1,v2,... [PARAM2 v1,v2,... ...]", file=sys.stderr)
        sys.exit(1)

    param_names = argv[0::2]
    value_lists = [csv.split(",") for csv in argv[1::2]]

    for param in param_names:
        if param in RESERVED_PARAMS:
            die(f"'{param}' is managed automatically by this script (it names the output dirs), sweep a different parameter.")

    return param_names, value_lists


def main():
    os.chdir(SCRIPT_DIR)
    param_names, value_lists = parse_args(sys.argv[1:])

    line_nos = []
    for param in param_names:
        line_no = find_assignment_line(param)
        if line_no is None:
            die(f"no assignment to '{param}' found in {PARAM_FILE}")
        line_nos.append(line_no)

    param_backup = tempfile.mkstemp()[1]
    shutil.copy(PARAM_FILE, param_backup)
    pro_backup = tempfile.mkstemp()[1]
    shutil.copy(PRO_FILE, pro_backup)

    try:
        # CellularPotts2.pro's TARGET decides which binary `make` produces
        # (and which parameter file it compiles in) - force it to
        # embryo_multi regardless of what's currently checked in, since this
        # script always needs the multi-organism binary.
        # newline='' + preserving each matched line's own trailing \r keeps
        # the file's existing line endings (CRLF in this repo) untouched -
        # read_text()/write_text()'s default universal-newline translation
        # would otherwise silently rewrite every line ending to LF.
        pro_text = PRO_FILE.read_text(newline="")
        pro_text = re.sub(r"^TARGET\s*=.*?(\r?)$", lambda m: "TARGET = embryo_multi" + m.group(1),
                           pro_text, flags=re.MULTILINE)
        PRO_FILE.write_text(pro_text, newline="")

        print("Regenerating Makefile for embryo_multi...")
        qmake_regen()

        if len(param_names) == 1:
            results_dir = SCRIPT_DIR / "sweep_results" / param_names[0]
        else:
            results_dir = SCRIPT_DIR / "sweep_results" / "+".join(param_names)
        results_dir.mkdir(parents=True, exist_ok=True)

        run_prefix = []
        if not os.environ.get("DISPLAY") and shutil.which("xvfb-run"):
            run_prefix = ["xvfb-run", "-a"]

        n_orgs = read_n_orgs()
        n_cores = os.cpu_count() or 1
        n_concurrent = max(1, n_cores // n_orgs) if n_orgs else 1
        if n_concurrent > 1:
            print(
                f"embryo_multi runs {n_orgs} organisms per process (OpenMP); {n_cores} cores available, "
                f"so building/running up to {n_concurrent} parameter combinations at once."
            )

        combos = list(itertools.product(*value_lists))
        for chunk_start in range(0, len(combos), n_concurrent):
            group = []
            for combo in combos[chunk_start:chunk_start + n_concurrent]:
                leaf_name, label = combo_leaf_and_label(param_names, combo)
                print(f"=== {label} ===")

                run_dir = results_dir / leaf_name
                if run_dir.exists():
                    shutil.rmtree(run_dir)
                run_dir.mkdir(parents=True)

                binary = build_combo(param_names, line_nos, param_backup, combo, label, run_dir)
                group.append((label, run_dir, binary))

            run_combo_group(run_prefix, group)

            for label, run_dir, binary in group:
                run_python(ANALYSE_SCRIPT, ["--csv", str(run_dir / "death_summary.csv"),
                                            "--start", "800", str(run_dir / "org-data")])
                print(f"Done: results in {run_dir}")

        run_python(AGGREGATE_SCRIPT, [str(results_dir)])
        run_python(REPORT_SCRIPT, [str(results_dir)])

        print(f"Sweep complete. Results under {results_dir}/")
    finally:
        restore(param_backup, pro_backup)


if __name__ == "__main__":
    main()

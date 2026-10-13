import os
import re
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_phase_diagram(
    base_dir='motchangebig',
    y_param=3,                 # 2 for (Param1 vs Param2), 3 for (Param1 vs Param3)
    fixed_value=None,          # Optional: target value for the remaining fixed parameter
    tolerance=1e-4,            # Float comparison tolerance
    output_file=None
):
    if y_param not in [2, 3]:
        raise ValueError("y_param must be either 2 or 3.")

    # Determine which parameter is held fixed
    fixed_param = 3 if y_param == 2 else 2
    y_col = f'Param{y_param}'
    fixed_col = f'Param{fixed_param}'

    phase_data = []

    # Regex: (param1)-(param2)-(param3)
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    if not os.path.exists(base_dir):
        print(f"Error: Directory '{base_dir}' does not exist.")
        return

    print(f"Scanning '{base_dir}' and reading data...")
    for folder_name in sorted(os.listdir(base_dir)):
        folder_path = os.path.join(base_dir, folder_name)

        if not os.path.isdir(folder_path):
            continue

        match = folder_pattern.match(folder_name)
        if not match:
            continue

        # Extract all three parameters
        x = float(match.group(1))
        p2 = round(float(match.group(2)), 4)
        p3 = round(float(match.group(3)), 4)

        # --- TRANSFORMATION FOR PARAMETER 1 ---
        k = (x + 0.3) / 1.1
        y = 2 * (1.2 - x - ((1.2 + 0.8 - (k * 1.6) + 0.4) / 2.0))
        param1_transformed = round(y, 2)

        # Check for data file (with typo fallback)
        file_path = os.path.join(folder_path, "migration_proportion.dat")
        if not os.path.exists(file_path):
            file_path = os.path.join(folder_path, "migartion_proportion.dat")

        if os.path.exists(file_path):
            try:
                with open(file_path, 'r') as f:
                    content = f.read().strip()
                    if content:
                        phase_data.append({
                            'Param1': param1_transformed,
                            'Param2': p2,
                            'Param3': p3,
                            'Proportion': float(content)
                        })
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

    if not phase_data:
        print("No valid data found. Check directory name and files.")
        return

    df = pd.DataFrame(phase_data)

    # --- FILTER / HANDLE FIXED PARAMETER ---
    if fixed_value is not None:
        df = df[df[fixed_col].apply(lambda v: math.isclose(v, fixed_value, abs_tol=tolerance))].copy()
        if df.empty:
            print(f"Error: No data matching Parameter {fixed_param} = {fixed_value}.")
            return
        disp_fixed = f"Param {fixed_param} = {fixed_value}"
    else:
        unique_vals = sorted(df[fixed_col].unique())
        if len(unique_vals) == 1:
            disp_fixed = f"Param {fixed_param} = {unique_vals[0]}"
        else:
            print(f"Notice: Multiple values found for fixed Parameter {fixed_param}: {unique_vals}")
            disp_fixed = f"Param {fixed_param} (multiple)"

    # --- CREATE HEATMAP PIVOT TABLE ---
    # Safe aggregation using mean in case of duplicate runs
    pivot_table = df.pivot_table(index=y_col, columns='Param1', values='Proportion', aggfunc='mean')

    # Sort Y descending (standard graph coordinates: high values top)
    pivot_table = pivot_table.sort_index(ascending=False)
    # Sort X ascending
    pivot_table = pivot_table.reindex(sorted(pivot_table.columns), axis=1)

    # --- PLOTTING ---
    print(f"Generating Phase Diagram (Y = Parameter {y_param}, X = Parameter 1)...")
    plt.figure(figsize=(9, 7))

    sns.heatmap(
        pivot_table,
        cmap='viridis',
        vmin=0.0,
        vmax=1.0,
        cbar_kws={'label': 'Migration Proportion'},
        square=True
    )

    plt.title(f'Phase Diagram ({base_dir}) | {disp_fixed}', pad=18, fontsize=13)
    plt.xlabel('Transformed Parameter 1', fontsize=12)
    plt.ylabel(f'Parameter {y_param}', fontsize=12)

    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300)
        print(f"Plot saved successfully as '{output_file}'")

    plt.show()


if __name__ == "__main__":
    # -------------------------------------------------------------
    # TOGGLE BETWEEN DATASETS HERE:
    # -------------------------------------------------------------

    # === DATASET 1: Param 1 vs Param 2 (Param 3 is constant) ===
    # plot_phase_diagram(
    #     base_dir='mot0',
    #     y_param=2,              # Y-axis is Param 2
    #     fixed_value=None,       # Or specify e.g. 0.2 if needed
    #     output_file='phase_diagram_dataset1.png'
    # )

    # === DATASET 2: Param 1 vs Param 3 (Param 2 is constant) ===
    plot_phase_diagram(
        base_dir='motchange',
        y_param=3,              # Y-axis is Param 3
        fixed_value=0.09,       # Fixed constant Param 2 value
        output_file='phase_diagram_dataset2.png'
    )
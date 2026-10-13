import os
import re
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_boundary_file(filepath, default_columns):
    """
    Reads a boundary length .dat file, handling both commented (#) 
    headers and headerless files.
    """
    with open(filepath, 'r') as f:
        first_line = f.readline()
    
    # Check if header is present in the first line
    if 'time' in first_line:
        cols = first_line.lstrip('#').strip().split()
        header_setting = 0 if not first_line.startswith('#') else None
        return pd.read_csv(filepath, sep=r'\s+', comment='#', names=cols, header=header_setting)
    else:
        return pd.read_csv(filepath, sep=r'\s+', comment='#', names=default_columns, header=None)


def plot_boundary_vs_param(
    base_dir='mot0big',
    target_col='sox2_loser',    # Column to plot ('sox2_med', 'sox2_loser', 'loser_med', 'total', etc.)
    last_n_steps=100,           # Number of final time steps to average over
    x_param=1,                  # Parameter on the X-axis (1, 2, or 3)
    line_param=2,               # Parameter representing each line (1, 2, or 3)
    fixed_value=0.0,            # Constant value for the remaining 3rd parameter
    line_values=None,           # Optional list: specific line values to plot (e.g., [0.1, 0.2]). None = plot all
    tolerance=1e-4,             # Float matching tolerance
    file_pattern="boundary_lengths-org-*.dat",
    output_file=None            # Set filename (e.g., 'boundary_param.png') to save
):
    # Validate parameter selections
    if x_param == line_param or not {x_param, line_param}.issubset({1, 2, 3}):
        raise ValueError("x_param and line_param must be distinct numbers chosen from {1, 2, 3}.")

    # Identify remaining constant parameter
    fixed_param = list({1, 2, 3} - {x_param, line_param})[0]
    print(fixed_param)
    x_col = f'Param{x_param}'
    line_col = f'Param{line_param}'
    fixed_col = f'Param{fixed_param}'

    # Default fallback columns if headers are missing
    default_columns = ['time', 'sox2_med', 'sox2_loser', 'loser_med', 'total']
    
    # Human-friendly labels for known columns
    col_labels = {
        'sox2_sox17': 'Sox2 / Sox17',
        'sox2_med':   'Sox2 / Medium',
        'sox2_loser': 'Sox2 / Loser',
        'sox17_med':  'Sox17 / Medium',
        'loser_med':  'Loser / Medium',
        'loser_sox17': 'Loser / Sox17',
        'total':      'Total Boundary'
    }
    label_name = col_labels.get(target_col, target_col)

    if not os.path.exists(base_dir):
        print(f"Error: Directory '{base_dir}' does not exist.")
        return

    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    # --- STEP 1: SCAN AND PRE-FILTER FOLDERS (OPTIMIZATION) ---
    matching_folders = []
    available_fixed_values = set()

    for folder_name in sorted(os.listdir(base_dir)):
        folder_path = os.path.join(base_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        match = folder_pattern.match(folder_name)
        if not match:
            continue

        x = float(match.group(1))
        p2 = float(match.group(2))
        p3 = float(match.group(3))

        k = (x + 0.3) / 1.1
        y = 2 * (1.2 - x - ((1.2 + 0.8 - (k * 1.6) + 0.4) / 2.0))
        p1 = round(y, 2)

        p_vals = {
            1: p1,
            2: p2,
            3: p3
        }
        available_fixed_values.add(p_vals[fixed_param])

        # Filter by fixed_param upfront to avoid reading unwanted data
        if not math.isclose(p_vals[fixed_param], fixed_value, abs_tol=tolerance):
            continue

        # Optional filter by specific line parameter values
        if line_values is not None:
            if not any(math.isclose(p_vals[line_param], target, abs_tol=tolerance) for target in line_values):
                continue

        matching_folders.append((folder_path, p_vals))

    if not matching_folders:
        print(f"\nNo folders matched Parameter {fixed_param} = {fixed_value}.")
        print(f"Available unique values for Parameter {fixed_param}:")
        print(sorted(available_fixed_values))
        return

    print(f"Found {len(matching_folders)} matching folders. Processing data...")

    # --- STEP 2: LOAD FILES AND COMPUTE FINAL AVERAGE ---
    all_data = []

    for folder_path, p_vals in matching_folders:
        files = glob.glob(os.path.join(folder_path, file_pattern))
        if not files:
            continue

        df_list = []
        for file in files:
            try:
                df = read_boundary_file(file, default_columns)
                if 'time' in df.columns and target_col in df.columns:
                    # Keep only required columns to reduce memory overhead
                    df_list.append(df[['time', target_col]])
            except Exception as e:
                print(f"Error reading {file}: {e}")

        if not df_list:
            continue

        # Combine all replicate runs for this parameter combination
        combined = pd.concat(df_list, ignore_index=True)

        # Average replicates across each time step, then average the last N steps
        mean_over_time = combined.groupby('time')[target_col].mean().sort_index()
        final_avg = mean_over_time.tail(last_n_steps).mean()

        all_data.append({
            x_col: p_vals[x_param],
            line_col: p_vals[line_param],
            fixed_col: p_vals[fixed_param],
            'FinalValue': final_avg
        })

    if not all_data:
        print(f"No valid data found for column '{target_col}'.")
        return

    df_plot = pd.DataFrame(all_data)

    # --- STEP 3: PLOTTING ---
    plt.figure(figsize=(9, 6))

    unique_lines = sorted(df_plot[line_col].unique())
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(unique_lines)))

    for val, color in zip(unique_lines, colors):
        subset = df_plot[df_plot[line_col] == val].sort_values(by=x_col)
        
        plt.plot(
            subset[x_col],
            subset['FinalValue'],
            marker='o',
            linestyle='-',
            linewidth=2,
            markersize=5,
            color=color,
            label=f'Param {line_param} = {val}'
        )

    # Formatting
    plt.title(
        f'{label_name} vs. Parameter {x_param}\n'
        f'(Grouped by Param {line_param} | Param {fixed_param} = {fixed_value} | Avg of last {last_n_steps} steps)',
        fontsize=13,
        pad=12
    )
    plt.xlabel(f'Parameter {x_param}', fontsize=12)
    plt.ylabel(f'{label_name} Boundary Length (px)', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)

    plt.legend(title=f'Parameter {line_param}', bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved successfully as '{output_file}'")

    plt.show()


if __name__ == "__main__":
    # Example usage:
    # - Directory: 'mot0'
    # - Metric: 'sox2_loser'
    # - X-axis: Parameter 1
    # - Separate curves: Parameter 2
    # - Constant: Parameter 3 held at 0.0
    plot_boundary_vs_param(
        # base_dir='mot0',
        # target_col='sox2_loser',
        # last_n_steps=100,
        # x_param=1,
        # line_param=2,
        # fixed_value=0.0
    )
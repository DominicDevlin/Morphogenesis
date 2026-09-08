import os
import re
import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_migration_multiline(
    base_dir='mot0',
    x_param=2,               # Parameter on the X-axis (1, 2, or 3)
    line_param=1,            # Parameter representing each line (1, 2, or 3)
    fixed_value=0.0,         # Constant value for the 3rd remaining parameter
    line_values=None,        # Optional list: specific line values to plot (e.g., [0.0, 0.2, 0.4]). None = plot all
    tolerance=1e-4,          # Float matching tolerance
    output_file='migration_multiline.png'
):
    # Validate parameter selections
    if x_param == line_param or not {x_param, line_param}.issubset({1, 2, 3}):
        raise ValueError("x_param and line_param must be distinct numbers chosen from {1, 2, 3}.")

    # Determine which parameter is held constant
    fixed_param = list({1, 2, 3} - {x_param, line_param})[0]
    
    x_col = f'Param{x_param}'
    line_col = f'Param{line_param}'
    fixed_col = f'Param{fixed_param}'

    # Regex to extract: (param1)-(param2)-(param3)
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    if not os.path.exists(base_dir):
        print(f"Error: Directory '{base_dir}' does not exist.")
        return

    all_data = []
    print(f"Scanning '{base_dir}' and reading data...")
    for folder_name in sorted(os.listdir(base_dir)):
        folder_path = os.path.join(base_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        match = folder_pattern.match(folder_name)
        if not match:
            continue

        p1 = float(match.group(1))
        p2 = float(match.group(2))
        p3 = float(match.group(3))

        # Check for data file (including typo fallback)
        file_path = os.path.join(folder_path, "migration_proportion.dat")
        if not os.path.exists(file_path):
            file_path = os.path.join(folder_path, "migartion_proportion.dat")

        if os.path.exists(file_path):
            try:
                with open(file_path, 'r') as f:
                    content = f.read().strip()
                    if content:
                        all_data.append({
                            'Param1': p1,
                            'Param2': p2,
                            'Param3': p3,
                            'Proportion': float(content)
                        })
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

    if not all_data:
        print("No valid data found in directory.")
        return

    df = pd.DataFrame(all_data)

    # --- FILTERING BY THE FIXED PARAMETER ---
    filtered_df = df[df[fixed_col].apply(lambda val: math.isclose(val, fixed_value, abs_tol=tolerance))].copy()

    if filtered_df.empty:
        print(f"\nNo data found where Parameter {fixed_param} = {fixed_value}.")
        print(f"Available unique values for Parameter {fixed_param}:")
        print(sorted(df[fixed_col].unique()))
        return

    # Filter specific lines if user supplied a list
    if line_values is not None:
        filtered_df = filtered_df[
            filtered_df[line_col].apply(lambda v: any(math.isclose(v, target, abs_tol=tolerance) for target in line_values))
        ]

    # --- PLOTTING ---
    plt.figure(figsize=(9, 6))

    # Get sorted unique values for the line parameter
    unique_lines = sorted(filtered_df[line_col].unique())
    
    # Use a colormap to assign distinct colors to each line
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(unique_lines)))

    for val, color in zip(unique_lines, colors):
        subset = filtered_df[filtered_df[line_col] == val].sort_values(by=x_col)
        
        plt.plot(
            subset[x_col],
            subset['Proportion'],
            marker='o',
            linestyle='-',
            linewidth=2,
            markersize=5,
            color=color,
            label=f'Param {line_param} = {val}'
        )

    # Labels and Titles
    plt.title(
        f'Migration Proportion vs. Parameter {x_param}\n'
        f'(Grouped by Parameter {line_param} | Parameter {fixed_param} = {fixed_value})',
        fontsize=13,
        pad=12
    )
    plt.xlabel(f'Parameter {x_param}', fontsize=12)
    plt.ylabel('Migration Proportion', fontsize=12)
    
    plt.ylim(-0.05, 1.05)
    plt.grid(True, linestyle='--', alpha=0.5)
    
    # Legend formatting
    plt.legend(title=f'Parameter {line_param}', bbox_to_anchor=(1.04, 1), loc="upper left")
    plt.tight_layout()

    # if output_file:
    #     plt.savefig(output_file, dpi=300, bbox_inches='tight')
    #     print(f"Plot saved successfully as '{output_file}'")

    plt.show()


if __name__ == "__main__":
    # Example:
    # - X-axis: Parameter 1
    # - Lines:  Parameter 2 (a curve for each value of Param 2)
    # - Fixed:  Parameter 3 held constant at 0.2
    plot_migration_multiline()
import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt

def plot_total_cells_param_sweep(base_dir="filtered-data", fixed_param="a", fixed_val=6):
    """
    Plots Total Cells over Time for an explored parameter while holding the other fixed.
    
    Parameters:
        base_dir (str): Directory containing the parameter folders (e.g., '6-0.3').
        fixed_param (str): Which parameter to hold constant ('a' or 'b').
        fixed_val (float): The specific value to fix 'fixed_param' to.
    """
    # Standardize input
    fixed_param = fixed_param.lower()
    if fixed_param not in ['a', 'b']:
        print("Error: fixed_param must be either 'a' or 'b'.")
        return

    # Determine which parameter is being explored (varied)
    var_param = 'b' if fixed_param == 'a' else 'a'

    # Regex to match folder names formatted as "a_val-b_val" (e.g., "6-0.3", "9--0.1", "-5--2")
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    if not os.path.exists(base_dir):
        print(f"Directory '{base_dir}' does not exist.")
        return

    # Dictionary to store mean time series: { var_val: pandas.Series }
    data_by_var = {}

    print(f"Scanning '{base_dir}' for folders with {fixed_param} = {fixed_val}...")

    # Iterate through folders in base directory
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)
        
        # Skip if not a directory
        if not os.path.isdir(folder_path):
            continue

        match = folder_pattern.match(folder_name)
        if not match:
            continue

        # Folder convention: group(1) = a_val, group(2) = b_val
        a_val = float(match.group(1))
        b_val = float(match.group(2))

        # Assign values based on chosen fixed parameter
        current_fixed = a_val if fixed_param == 'a' else b_val
        current_var = b_val if fixed_param == 'a' else a_val

        # Filter by target fixed_val (using small tolerance for float equality)
        if abs(current_fixed - fixed_val) > 1e-5:
            continue

        # Find all matching replicate data files in this parameter directory
        file_pattern = os.path.join(folder_path, "*org-*.dat")
        files = glob.glob(file_pattern)

        if not files:
            print(f"No matching files found in: {folder_path}")
            continue

        # Read and combine all replicate files for this folder
        df_list = []
        for file in files:
            try:
                # sep=r'\s+' handles tabs or space-separated files
                df = pd.read_csv(file, sep=r'\s+')
                df_list.append(df)
            except Exception as e:
                print(f"Error reading {file}: {e}")

        if not df_list:
            continue

        all_data = pd.concat(df_list, ignore_index=True)

        # Calculate mean 'total' cells grouped by 'time'
        mean_data = all_data.groupby('time')['total'].mean()
        data_by_var[current_var] = mean_data

    if not data_by_var:
        print(f"No valid data found matching {fixed_param} = {fixed_val}.")
        return

    # --- PLOTTING ---
    plt.figure(figsize=(10, 6))

    # Sort keys so lines are added in numerical order of the explored parameter
    for var_val in sorted(data_by_var.keys()):
        mean_series = data_by_var[var_val]
        plt.plot(
            mean_series.index, 
            mean_series.values, 
            linewidth=2, 
            label=f"{var_param} = {var_val}"
        )

    # Plot formatting (dynamically updates titles and labels based on choice)
    plt.title(f"Average Total Cells over Time (Fixed {fixed_param} = {fixed_val})", fontsize=14)
    plt.xlabel("Time", fontsize=12)
    plt.ylabel("Total Cells", fontsize=12)
    plt.legend(loc="best", title=f"Parameter {var_param}")
    plt.xlim(left=0)
    plt.grid(True, linestyle="--", alpha=0.5)

    # Save and display plot
    plt.tight_layout()
    output_filename = f"total_cells_fixed_{fixed_param}_{fixed_val}.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved successfully as '{output_filename}'")
    plt.show()

if __name__ == "__main__":
    # Example 1: Fix b_val at 0.3 and plot multiple lines for different a_val options
    plot_total_cells_param_sweep(base_dir="filtered-data", fixed_param="b", fixed_val=0.3)

    # Example 2: Fix a_val at 6.0 and plot multiple lines for different b_val options
    # plot_total_cells_param_sweep(base_dir="filtered-data", fixed_param="a", fixed_val=6.0)
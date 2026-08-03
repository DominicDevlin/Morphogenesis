import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt

def plot_boundary_lengths_sweep(base_dir="newest-data/", 
                                fixed_param="b", 
                                fixed_val=0.4, 
                                plot_type="both"):
    """
    Plots Boundary Lengths over Time for an explored parameter while holding the other fixed.
    
    Parameters:
        base_dir (str): Directory containing the parameter folders (e.g., '6-0.3').
        fixed_param (str): Which parameter to hold constant ('a' or 'b').
        fixed_val (float): The specific value to fix 'fixed_param' to.
        plot_type (str): What to plot. Options: 'loserwinner', 'sox2sox17', or 'both'.
    """
    
    # Standardize inputs
    fixed_param = fixed_param.lower()
    plot_type = plot_type.lower()
    
    if fixed_param not in ['a', 'b']:
        print("Error: fixed_param must be either 'a' or 'b'.")
        return
    if plot_type not in ['loserwinner', 'sox2sox17', 'both']:
        print("Error: plot_type must be 'loserwinner', 'sox2sox17', or 'both'.")
        return

    # Determine which parameter is being explored (varied)
    var_param = 'b' if fixed_param == 'a' else 'a'

    # Regex to match folder names formatted as "a_val-b_val"
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    if not os.path.exists(base_dir):
        print(f"Directory '{base_dir}' does not exist.")
        return

    # Dictionary to store mean time series: 
    # { var_val: {'loserwinner': pd.Series, 'sox2sox17': pd.Series} }
    data_by_var = {}

    print(f"Scanning '{base_dir}' for folders with {fixed_param} = {fixed_val}...")

    # Iterate through folders in base directory
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)
        
        if not os.path.isdir(folder_path):
            continue

        match = folder_pattern.match(folder_name)
        if not match:
            continue

        # Extract parameters from folder name
        a_val = float(match.group(1))
        b_val = float(match.group(2))

        # Assign values based on chosen fixed parameter
        current_fixed = a_val if fixed_param == 'a' else b_val
        current_var = b_val if fixed_param == 'a' else a_val

        # Filter by target fixed_val (using small tolerance for float equality)
        if abs(current_fixed - fixed_val) > 1e-5:
            continue

        # Find all matching boundary length replicate files
        file_pattern = os.path.join(folder_path, "boundary_lengths-org-*.dat")
        files = glob.glob(file_pattern)

        if not files:
            continue

        # Read and combine all replicate files for this folder
        df_list = []
        for file in files:
            try:
                df = pd.read_csv(file, sep=r'\s+')
                df_list.append(df)
            except Exception as e:
                print(f"Error reading {file}: {e}")

        if not df_list:
            continue

        all_data = pd.concat(df_list, ignore_index=True)

        # Calculate mean boundary lengths grouped by 'time'
        mean_lw = all_data.groupby('time')['loserwinner'].mean()
        mean_sox = all_data.groupby('time')['sox2sox17'].mean()
        
        # Store in our dictionary
        data_by_var[current_var] = {
            'loserwinner': mean_lw,
            'sox2sox17': mean_sox
        }

    if not data_by_var:
        print(f"No valid data found matching {fixed_param} = {fixed_val}.")
        return

    # --- PLOTTING ---
    
    # Sort keys so legend reflects numerical order of the explored parameter
    sorted_vars = sorted(data_by_var.keys())

    if plot_type == "both":
        # Create a side-by-side plot for both metrics
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        for var_val in sorted_vars:
            d = data_by_var[var_val]
            ax1.plot(d['loserwinner'].index, d['loserwinner'].values, linewidth=2, label=f"{var_param} = {var_val}")
            ax2.plot(d['sox2sox17'].index, d['sox2sox17'].values, linewidth=2, label=f"{var_param} = {var_val}")

        ax1.set_title("Average Loser/Winner Boundary Length", fontsize=14)
        ax1.set_xlabel("Time", fontsize=12)
        ax1.set_ylabel("Boundary Length", fontsize=12)
        ax1.set_xlim(left=800)
        ax1.legend(title=f"Parameter {var_param}")
        ax1.set_ylim(0.3)

        
        ax2.set_title("Average Sox2/Sox17 Boundary Length", fontsize=14)
        ax2.set_xlabel("Time", fontsize=12)
        ax2.set_ylabel("Boundary Length", fontsize=12)
        ax2.set_xlim(left=800)
        ax2.legend(title=f"Parameter {var_param}")

        plt.suptitle(f"Boundary Lengths over Time (Fixed {fixed_param} = {fixed_val})", fontsize=16)

    else:
        # Create a single plot for just one metric
        plt.figure(figsize=(10, 6))
        for var_val in sorted_vars:
            d = data_by_var[var_val]
            series_to_plot = d[plot_type]
            plt.plot(series_to_plot.index, series_to_plot.values, linewidth=2, label=f"{var_param} = {var_val}")
        
        title_metric = "Loser/Winner" if plot_type == 'loserwinner' else "Sox2/Sox17"
        plt.title(f"Average {title_metric} Boundary Length (Fixed {fixed_param} = {fixed_val})", fontsize=14)
        plt.xlabel("Time", fontsize=12)
        plt.ylabel("Boundary Length", fontsize=12)
        plt.xlim(left=800)
        plt.legend(title=f"Parameter {var_param}", loc="best")

    # Save and display plot
    plt.tight_layout()
    output_filename = f"boundary_lengths_{plot_type}_fixed_{fixed_param}_{fixed_val}.png"
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved successfully as '{output_filename}'")
    plt.show()


if __name__ == "__main__":
    # Example 1: Plot BOTH metrics side-by-side (fixes b at 0.3, varies a)
    plot_boundary_lengths_sweep()

    # Example 2: Plot ONLY Loser/Winner boundary lengths
    # plot_boundary_lengths_sweep(base_dir="newest-data", fixed_param="b", fixed_val=0.3, plot_type="loserwinner")

    # Example 3: Plot ONLY Sox2/Sox17 boundary lengths (fixes a at 6.0, varies b)
    # plot_boundary_lengths_sweep(base_dir="newest-data", fixed_param="a", fixed_val=6.0, plot_type="sox2sox17")
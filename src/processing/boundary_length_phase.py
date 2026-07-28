import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress warnings when calculating mean of empty slices
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

def plot_phase_diagram_sox(base_dir='newest-data/'):
    # Data structure to hold aggregated results
    phase_data = []

    # Regex to match folder names like "6-0.3", "9--0.1", "-5--2"
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    print("Scanning directories and processing boundary length files...")
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)
        
        # Skip if it's not a directory
        if not os.path.isdir(folder_path):
            continue
            
        # Extract X and Y parameter values from folder name
        match = folder_pattern.match(folder_name)
        if not match:
            continue
            
        x_val = float(match.group(1))
        y_val = float(match.group(2))
        
        replicate_means = []
        
        # Loop through up to 60 replicate files
        for i in range(1, 61):
            file_name = f"boundary_lengths-org-{i}.dat"
            file_path = os.path.join(folder_path, file_name)
            
            if os.path.exists(file_path):
                try:
                    # Read boundary length file using whitespace as separator
                    df = pd.read_csv(file_path, sep=r'\s+')
                    
                    if 'loserwinner' in df.columns and not df.empty:
                        # Extract the last 100 rows (time points) for 'sox2sox17'
                        # (If fewer than 100 rows exist, .tail(100) safely takes all available)
                        last_100_sox = df['loserwinner'].tail(100)
                        
                        # Compute the average over the last 100 rows for this replicate
                        rep_mean = last_100_sox.mean()
                        replicate_means.append(rep_mean)
                    else:
                        replicate_means.append(np.nan)
                        
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")
        
        # Calculate average across all valid replicates for this (X, Y) parameter pair
        if replicate_means:
            avg_sox = np.nanmean(replicate_means)
            phase_data.append({
                'X': x_val, 
                'Y': y_val, 
                'Sox2Sox17': avg_sox
            })

    if not phase_data:
        print("No valid data found. Please check your directory path and folder names.")
        return

    # Convert processed results to DataFrame
    df_plot = pd.DataFrame(phase_data)
    
    # Pivot DataFrame into 2D grid format: Y values as rows, X values as columns
    pivot_table = df_plot.pivot(index='Y', columns='X', values='Sox2Sox17')
    
    # Sort Y-axis descending so larger values appear at top of plot
    pivot_table = pivot_table.sort_index(ascending=False)

    # --- PLOTTING ---
    print("Generating Phase Diagram...")
    plt.figure(figsize=(10, 8))
    
    # Heatmap color scaling will automatically adapt to min/max of Sox2Sox17 values
    ax = sns.heatmap(
        pivot_table, 
        cmap='viridis', 
        vmin=0.4,
        cbar_kws={'label': 'Mean Sox2/Sox17 Boundary Length (Last 100 Time Points)'},
        annot=False, 
        square=True
    )

    plt.title('Phase Diagram: Average Sox2/Sox17 Boundary Length', pad=20, fontsize=14)
    plt.xlabel('Parameter X', fontsize=12)
    plt.ylabel('Parameter Y', fontsize=12)
    
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)

    plt.tight_layout()
    
    output_filename = 'phase_diagram_sox2sox17.png'
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved successfully as '{output_filename}'")
    plt.show()

if __name__ == "__main__":
    plot_phase_diagram_sox()
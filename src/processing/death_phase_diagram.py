import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress warnings from calculating mean of empty arrays (when all deaths are 0)
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

def plot_phase_diagram(base_dir='t1/'):
    # Data structure to hold the aggregated results
    # List of dictionaries: [{'X': x_val, 'Y': y_val, 'Proportion': avg_prop}, ...]
    phase_data = []

    # Regex to match folder names like "6-0.3", "9--0.1", "-5--2"
    # Matches: (number) - (number). Numbers can optionally start with + or -
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    # Iterate through all items in the base directory
    print("Scanning directories and processing files...")
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)
        print(folder_path)
        # Skip if it's not a directory
        if not os.path.isdir(folder_path):
            continue
            
        # Check if the folder matches our expected naming convention
        match = folder_pattern.match(folder_name)
        if not match:
            continue
            
        # Extract X and Y values from the folder name
        x_val = float(match.group(1))
        y_val = float(match.group(2))
        
        replicate_proportions = []
        
        # Loop through the 60 replicate files
        for i in range(1, 61):
            file_name = f"death_causes-org-{i}.dat"
            file_path = os.path.join(folder_path, file_name)
            
            if os.path.exists(file_path):
                try:
                    # Read the file. sep='\s+' handles tabs or variable spaces
                    df = pd.read_csv(file_path, sep='\s+')
                    
                    # Get the final line (last time point)
                    last_row = df.iloc[-1]
                    
                    # Extract the required values
                    lonely = last_row['undifferentiated_lonely']
                    signal = last_row['undifferentiated_signal']
                    
                    total_undiff = lonely + signal
                    
                    # Calculate proportion. Guard against division by zero.
                    if total_undiff > 0:
                        proportion = lonely / total_undiff
                        replicate_proportions.append(proportion)
                    else:
                        # If no undifferentiated deaths occurred, we append NaN
                        # so it doesn't skew the average of actual events.
                        replicate_proportions.append(np.nan)
                        
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")
        
        # Calculate the average proportion across all replicates for this parameter pair
        if replicate_proportions:
            avg_prop = np.nanmean(replicate_proportions)
            phase_data.append({
                'X': x_val, 
                'Y': y_val, 
                'Proportion': avg_prop
            })

    if not phase_data:
        print("No valid data found. Please check your directory path and folder names.")
        return

    # Convert the processed data into a Pandas DataFrame
    df_plot = pd.DataFrame(phase_data)
    
    # Create a 2D pivot table for the heatmap: Index is Y, Columns are X, Values are Proportion
    # This automatically sorts the X and Y axes numerically
    pivot_table = df_plot.pivot(index='Y', columns='X', values='Proportion')
    
    # Sort the Y-axis so larger values are at the top of the plot (standard graphing behavior)
    pivot_table = pivot_table.sort_index(ascending=False)

    # --- PLOTTING ---
    print("Generating Phase Diagram...")
    plt.figure(figsize=(10, 8))
    
    # Create the heatmap
    # vmin=0 and vmax=1 keep the color scale strictly between 0 and 1 (as it's a proportion)
    # cmap='viridis' is a colorblind-friendly, visually appealing sequential colormap
    ax = sns.heatmap(pivot_table, cmap='viridis', vmin=0, vmax=1, 
                     cbar_kws={'label': 'Proportion: Lonely / (Lonely + Signal)'},
                     annot=False, square=True)

    plt.title('Phase Diagram: Undifferentiated Death Causes at Final Time Point', pad=20, fontsize=14)
    plt.xlabel('Parameter X', fontsize=12)
    plt.ylabel('Parameter Y', fontsize=12)
    
    # Rotate the x-axis labels if they overlap
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)

    # Save and show the plot
    plt.tight_layout()
    plt.savefig('phase_diagram.png', dpi=300)
    print("Plot saved as 'phase_diagram.png'")
    plt.show()

if __name__ == "__main__":
    # If the script is placed in the folder containing '6-0.3', '9--0.1', etc.
    # Just run it as is. Otherwise, replace '.' with your target directory path.
    plot_phase_diagram()
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress warnings from calculating mean of empty arrays
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

def plot_normalized_phase_diagram(base_dir='filtered-data/', target_time=10000.0):
    # Data structure to hold the raw aggregated results
    raw_data = []

    # Regex to match folder names like "3-0.3", "9--0.1", "-5--2"
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    print(f"Scanning directories and processing files for time = {target_time}...")
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)
        
        if not os.path.isdir(folder_path):
            continue
            
        match = folder_pattern.match(folder_name)
        if not match:
            continue
            
        a_val = float(match.group(1))
        b_val = float(match.group(2))
        
        replicate_totals = []
        
        # Loop through the 60 replicate files
        for i in range(1, 61):
            file_name = f"celltypes-org-{i}.dat"
            file_path = os.path.join(folder_path, file_name)
            
            if os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path, sep='\s+')
                    
                    # Find the row where the 'time' column matches the target_time.
                    # np.isclose is used here to avoid issues with floating-point precision.
                    matched_rows = df[np.isclose(df['time'], target_time, atol=1e-2)]
                    
                    if not matched_rows.empty:
                        # Extract the 'total' value from the first matching row
                        total_val = matched_rows.iloc[0]['total']
                        replicate_totals.append(total_val)
                    else:
                        # Log if the specific time was not reached or found in this file
                        pass
                        
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")
        
        # Calculate the average total across all valid replicates for this parameter pair
        if replicate_totals:
            avg_total = np.nanmean(replicate_totals)
            raw_data.append({
                'A': a_val, 
                'B': b_val, 
                'Avg_Total': avg_total
            })
            print(a_val, b_val, avg_total)


    if not raw_data:
        print(f"No valid data found for time = {target_time}. Please check the time value and directory path.")
        return

    # Convert raw data into a DataFrame
    df_raw = pd.DataFrame(raw_data)
    
    # --- NORMALIZATION STEP ---
    print("Normalizing data against A = 3...")
    
    # Isolate rows where A is 3
    ref_df = df_raw[np.isclose(df_raw['A'], 24.0)]
    print(ref_df)
    
    # Map B -> Avg_Total (where A=3)
    ref_dict = dict(zip(ref_df['B'], ref_df['Avg_Total']))
    
    def normalize_total(row):
        b = row['B']
        if b in ref_dict and ref_dict[b] != 0:
            return row['Avg_Total'] / ref_dict[b]
        else:
            return np.nan
            
    df_raw['Normalized_Total'] = df_raw.apply(normalize_total, axis=1)

    # --- PIVOTING ---
    pivot_table = df_raw.pivot(index='B', columns='A', values='Normalized_Total')
    pivot_table = pivot_table.sort_index(ascending=False)

    # --- PLOTTING ---
    print("Generating Phase Diagram...")
    plt.figure(figsize=(10, 8))
    
    # Using 'coolwarm' centered at 1.0
    ax = sns.heatmap(pivot_table, cmap='coolwarm', center=1, 
                     cbar_kws={'label': f'Normalized Total (Relative to a=3 at t={target_time})'}, vmin=0.75,
                     annot=False, square=True)

    plt.title(f'Phase Diagram: Normalized Cell Totals at Time = {target_time}', pad=20, fontsize=14)
    plt.xlabel('Parameter a', fontsize=12)
    plt.ylabel('Parameter b', fontsize=12)
    
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)

    plt.tight_layout()
    plt.savefig(f'normalized_phase_diagram_t{int(target_time)}.png', dpi=300)
    print(f"Plot saved as 'normalized_phase_diagram_t{int(target_time)}.png'")
    plt.show()

if __name__ == "__main__":
    # Specify the target time point here (e.g., 10000)
    plot_normalized_phase_diagram()
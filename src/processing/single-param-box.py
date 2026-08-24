import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress warnings
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

def plot_a_boxplot(base_dir='newest-data/', target_time=12000.0, target_b=0.4):
    # List to hold individual replicate values
    raw_data = []

    # Regex to match folder names like "3-0.3", "9--0.1", "-5--2"
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    print(f"Scanning directories for B = {target_b} at time = {target_time}...")
    
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)
        
        if not os.path.isdir(folder_path):
            continue
            
        match = folder_pattern.match(folder_name)
        if not match:
            continue
            
        a_val = float(match.group(1))
        b_val = float(match.group(2))
        
        # Filter: Only process directories matching the target value of parameter B
        if not np.isclose(b_val, target_b, atol=1e-2):
            continue
            
        # Loop through all 60 replicate files
        for i in range(1, 61):
            file_name = f"celltypes-org-{i}.dat"
            file_path = os.path.join(folder_path, file_name)
            
            if os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path, sep=r'\s+')
                    
                    # Match target time
                    matched_rows = df[np.isclose(df['time'], target_time, atol=1e-2)]
                    
                    if not matched_rows.empty:
                        total_val = matched_rows.iloc[0]['total']
                        # Collect individual replicate point
                        raw_data.append({
                            'A': a_val, 
                            'B': b_val, 
                            'Replicate': i,
                            'Total': total_val
                        })
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")

    if not raw_data:
        print(f"No valid data found for B = {target_b} at time = {target_time}. Check directory path or parameter target.")
        return

    # Convert raw data into a DataFrame
    df_raw = pd.DataFrame(raw_data)

    # Sort numerical values of A for proper x-axis ordering
    df_raw = df_raw.sort_values(by='A')

    # --- PLOTTING ---
    print("Generating Box Plot...")
    plt.figure(figsize=(10, 6))
    
    # Create box plot: A on x-axis, raw Total values on y-axis
    sns.boxplot(data=df_raw, x='A', y='Total', color='lightskyblue', width=0.5)

    # Styling
    plt.title(f'Cell Totals Distribution Across Parameter A (B = {target_b}, t = {target_time})', pad=15, fontsize=14)
    plt.xlabel('Parameter A', fontsize=12)
    plt.ylabel(f'Total Count at t = {target_time}', fontsize=12)
    
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    
    # Save image
    output_filename = f'boxplot_A_for_B{target_b}_t{int(target_time)}.png'
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved as '{output_filename}'")
    plt.show()

if __name__ == "__main__":
    # Specify target parameter B and target time point here
    plot_a_boxplot()
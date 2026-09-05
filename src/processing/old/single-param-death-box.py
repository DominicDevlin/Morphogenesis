import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

# Suppress warnings
warnings.filterwarnings(action='ignore', message='Mean of empty slice')

def plot_death_causes_boxplot(base_dir='newest-data/', target_b=0.3):
    # Data structure to hold individual replicate results
    raw_data = []

    # Regex to match folder names like "6-0.3", "9--0.1", "-5--2"
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    print(f"Scanning directories for Parameter B (Y) = {target_b}...")
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)
        
        if not os.path.isdir(folder_path):
            continue
            
        match = folder_pattern.match(folder_name)
        if not match:
            continue
            
        # Extract X (Parameter A) and Y (Parameter B) values
        x_val = float(match.group(1))
        y_val = float(match.group(2))
        
        # Filter: Only process directories matching the target value of Parameter B (Y)
        if not np.isclose(y_val, target_b, atol=1e-2):
            continue
            
        # Loop through the 60 replicate files
        for i in range(1, 61):
            file_name = f"death_causes-org-{i}.dat"
            file_path = os.path.join(folder_path, file_name)
            
            if os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path, sep=r'\s+')
                    
                    # Get the final line (last time point)
                    last_row = df.iloc[-1]
                    
                    lonely = last_row['undifferentiated_lonely']
                    signal = last_row['undifferentiated_signal']
                    
                    total_undiff = lonely + signal
                    
                    # Calculate proportion for this replicate
                    if total_undiff > 0:
                        proportion = lonely / total_undiff
                        raw_data.append({
                            'X': x_val, # Parameter A
                            'Y': y_val, # Parameter B
                            'Replicate': i,
                            'Proportion': proportion
                        })
                        
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")

    if not raw_data:
        print(f"No valid data found for Parameter B = {target_b}. Check directory path and folder names.")
        return

    # Convert the raw data into a DataFrame
    df_raw = pd.DataFrame(raw_data)

    # Sort numerical values of X (Parameter A) for proper x-axis ordering
    df_raw = df_raw.sort_values(by='X')

    # --- PLOTTING ---
    print("Generating Box Plot...")
    plt.figure(figsize=(10, 6))
    
    # Create box plot: Parameter A (X) on x-axis, Proportion on y-axis
    sns.boxplot(data=df_raw, x='X', y='Proportion', color='lightskyblue', width=0.5)

    # Styling
    plt.title(f'Death Causes Proportion across Parameter A (B = {target_b})', pad=15, fontsize=14)
    plt.xlabel('Parameter A (X)', fontsize=12)
    plt.ylabel('Proportion: Lonely / (Lonely + Signal)', fontsize=12)
    
    # Set y-axis bounds strictly between 0 and 1 since it is a proportion
    plt.ylim(-0.05, 1.05)
    
    plt.xticks(rotation=45)
    #plt.grid(axis='y', linestyle='--', alpha=0.7)

    plt.tight_layout()
    
    output_filename = f'boxplot_death_causes_B{target_b}.png'
    plt.savefig(output_filename, dpi=300)
    print(f"Plot saved as '{output_filename}'")
    plt.show()

if __name__ == "__main__":
    # Specify the target parameter B (Y) value here
    plot_death_causes_boxplot(target_b=0.3)
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

def plot_average_total_cells(data_dir="mot0.0/0.3-0.12-0.0"):
    # 1. Find all relevant .dat files in the directory
    # Adjust the pattern if your files have a specific prefix (e.g., "celltypes-org-*.dat")
    file_pattern = os.path.join(data_dir, "boundary_lengths-org-*.dat")
    print(file_pattern)
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No files matching '{file_pattern}' found in directory: {data_dir}")
        return

    print(f"Found {len(files)} files. Processing...")

    # 2. Read all files and combine them into a single DataFrame
    df_list = []
    for file in files:
        # \s+ is used as separator to handle any combination of spaces or tabs
        df = pd.read_csv(file, sep='\s+') 
        df_list.append(df)
        
    # Concatenate all dataframes together
    all_data = pd.concat(df_list, ignore_index=True)

    # 3. Group the data by 'time' and calculate statistics
    # Calculate the mean for both columns
    mean_data_wl = all_data.groupby('time')['loserwinner'].mean()
    mean_data_sox = all_data.groupby('time')['sox2sox17'].mean()
    
    # Calculate Standard Deviation (for plotting) for both columns
    std_data_wl = all_data.groupby('time')['loserwinner'].std()
    std_data_sox = all_data.groupby('time')['sox2sox17'].std()

    # Fill NaN values with 0 in case there's only 1 file (which makes std NaN)
    std_data_wl = std_data_wl.fillna(0)
    std_data_sox = std_data_sox.fillna(0)

    # 4. Create the plot
    plt.figure(figsize=(10, 6))
    
    # --- PLOT 1: Loser/Winner ---
    # Plot the mean line
    plt.plot(mean_data_wl.index, mean_data_wl.values, 
             label='Average Loser/Winner', color='blue', linewidth=2)
    
    # Plot the shaded region for Standard Deviation
    plt.fill_between(mean_data_wl.index, 
                     mean_data_wl.values - std_data_wl.values, 
                     mean_data_wl.values + std_data_wl.values, 
                     color='blue', alpha=0.2, 
                     label='± 1 Std Dev (Loser/Winner)')

    # --- PLOT 2: Sox2/Sox17 ---
    # Plot the mean line
    plt.plot(mean_data_sox.index, mean_data_sox.values, 
             label='Average Sox2/Sox17', color='orange', linewidth=2)
    
    # Plot the shaded region for Standard Deviation
    plt.fill_between(mean_data_sox.index, 
                     mean_data_sox.values - std_data_sox.values, 
                     mean_data_sox.values + std_data_sox.values, 
                     color='orange', alpha=0.2, 
                     label='± 1 Std Dev (Sox2/Sox17)')

    # Formatting the plot
    plt.title('Average Cell Populations over Time Across Replicates', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Cell Count / Ratio', fontsize=12) # Adjusted label to be more generic
    
    # Place legend outside or upper left so it doesn't cover data
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    
    # Ensure time axis starts at 0 
    plt.xlim(left=0)

    # Show the plot
    plt.tight_layout() # Ensures everything fits, especially if legend is moved
    plt.xlim(800)
    plt.show()

if __name__ == "__main__":
    plot_average_total_cells()
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

def plot_average_total_cells(data_dir="newest-data/2.5-0.3"):
    # 1. Find all relevant .dat files in the directory
    # Adjust the pattern if your files have a specific prefix (e.g., "celltypes-org-*.dat")
    file_pattern = os.path.join(data_dir, "*org-*.dat")
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

    # 3. Group the data by 'time' and calculate statistics for 'total'
    # Calculate the mean
    mean_data = all_data.groupby('time')['total'].mean()
    
    # Calculate Variance (as requested) and Standard Deviation (for plotting)
    var_data = all_data.groupby('time')['total'].var()
    std_data = all_data.groupby('time')['total'].std()

    # Fill NaN values with 0 in case there's only 1 file (which makes std/var NaN)
    std_data = std_data.fillna(0)

    # 4. Create the plot
    plt.figure(figsize=(10, 6))
    
    # Plot the mean line
    plt.plot(mean_data.index, mean_data.values, 
             label='Average Total Cells', color='blue', linewidth=2)
    
    # Plot the variance as a shaded region using Standard Deviation
    # (Since variance is "cells squared", standard deviation gives the correct Y-axis units)
    plt.fill_between(mean_data.index, 
                     mean_data.values - std_data.values, 
                     mean_data.values + std_data.values, 
                     color='blue', alpha=0.2, 
                     label='± 1 Std Dev (Variance spread)')

    # Formatting the plot
    plt.title('Average Total Cells over Time Across Replicates', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Total Cells', fontsize=12)
    plt.legend(loc='upper left')
    
    # Ensure axes start at 0 if appropriate for your data
    #plt.ylim(bottom=0) 
    plt.xlim(left=0)

    # Show the plot
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # If your scripts are in the same folder as the script, leave as "."
    # Otherwise, replace "." with the path to your data folder, e.g., "./data/"
    plot_average_total_cells()
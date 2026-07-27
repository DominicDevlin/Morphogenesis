import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

def plot_individual_total_cells(data_dir="t1/12--0.1/"):
    # 1. Find all relevant .dat files in the directory
    file_pattern = os.path.join(data_dir, "celltypes-org-*.dat")
    print(file_pattern)
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No files matching '{file_pattern}' found in directory: {data_dir}")
        return

    print(f"Found {len(files)} files. Processing and plotting...")

    # 2. Create the plot
    plt.figure(figsize=(10, 6))
    
    # 3. Read each file and plot its data individually
    for i, file in enumerate(files):
        file_name = os.path.basename(file)
        try:
            # \s+ handles space or tab delimiters
            df = pd.read_csv(file, sep='\s+') 
            
            if 'time' in df.columns and 'total' in df.columns:
                # To avoid cluttering the legend, we only assign a label to the first line.
                # If you prefer to show every file name, use: label=file_name
                label = 'Individual Replicates' if i == 0 else ""
                
                plt.plot(df['time'], df['total'], alpha=0.7, linewidth=1.5, label=label)
            else:
                print(f"Warning: 'time' or 'total' columns not found in {file_name}")
                
        except Exception as e:
            print(f"Error reading {file_name}: {e}")

    # Formatting the plot
    plt.title('Total Cells over Time for Individual Replicates', fontsize=14)
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('Total Cells', fontsize=12)
    plt.legend(loc='upper left')
    
    # Ensure axes start at 0 if appropriate
    plt.xlim(left=0)

    # Show the plot
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_individual_total_cells()
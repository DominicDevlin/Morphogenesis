import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def plot_average_boundary_lengths(data_dir="motchange/0.3-0.06-0.72"):
    # 1. Find all relevant .dat files in the directory
    file_pattern = os.path.join(data_dir, "boundary_lengths-org-*.dat")
    files = glob.glob(file_pattern)
    
    if not files:
        print(f"No files matching '{file_pattern}' found in directory: {data_dir}")
        return

    print(f"Found {len(files)} files. Processing...")

    # Expected column names based on the C++ header
    columns = [
        'time', 
        # 'sox2_sox17', 
        'sox2_med', 
        'sox2_loser', 
        'sox17_med', 
        'loser_med', 
        'loser_sox17', 
        'total'
    ]

    # 2. Read all files and combine them
    df_list = []
    for file in files:
        # Read the file. If your C++ code already writes the header line, remove names=columns, header=0
        # If your C++ file DOES NOT write headers, keep names=columns, header=None
        df = pd.read_csv(file, sep=r'\s+', comment='#')
        
        # Fallback: assign columns if header wasn't present in the file
        if len(df.columns) == len(columns) and not isinstance(df.columns[0], str):
            df.columns = columns
            
        df_list.append(df)
        
    all_data = pd.concat(df_list, ignore_index=True)

    # 3. Group by 'time' and compute mean and std simultaneously
    grouped = all_data.groupby('time')
    means = grouped.mean()
    stds = grouped.std().fillna(0)

    # 4. Create the plot
    plt.figure(figsize=(11, 7))
    
    # Boundary pairs to plot with descriptive labels and distinct colors
    boundaries_to_plot = [
        # ('sox2_sox17',   'Sox2 / Sox17',   '#e41a1c'),  # Red
        ('sox2_med',     'Sox2 / Medium',  '#377eb8'),  # Blue
        ('sox2_loser',   'Sox2 / Loser',   '#4daf4a'),  # Green
        ('sox17_med',    'Sox17 / Medium', '#984ea3'),  # Purple
        ('loser_med',    'Loser / Medium', '#ff7f00'),  # Orange
        ('loser_sox17',  'Loser / Sox17',  '#a65628')   # Brown
    ]

    for col, label, color in boundaries_to_plot:
        if col not in means.columns:
            continue
        
        # Plot Mean Line
        plt.plot(means.index, means[col], label=label, color=color, linewidth=2)
        
        # Plot Standard Deviation Envelope
        plt.fill_between(
            means.index, 
            means[col] - stds[col], 
            means[col] + stds[col], 
            color=color, 
            alpha=0.15
        )

    # Formatting
    plt.title('Cell-Cell & Cell-Medium Boundary Lengths Over Time', fontsize=14, fontweight='bold')
    plt.xlabel('Monte Carlo Steps (Time)', fontsize=12)
    plt.ylabel('Boundary Length (pixels)', fontsize=12)
    # plt.grid(True, linestyle='--', alpha=0.5)

    # Time boundaries (starts at 0, goes up to 800 or autoscale)
    plt.xlim(left=800)#, right=800)
    plt.ylim(bottom=0)

    # Legend outside to prevent overlapping lines
    plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1), frameon=True)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    plot_average_boundary_lengths()
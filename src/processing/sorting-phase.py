import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_phase_diagram(base_dir='mot0.2', output_file='phase_diagram.png'):
    phase_data = []

    # Regex to match folder names like "0.4-0.0-0.2" or "-0.5--0.1-0.2"
    # Matches: (param1)-(param2)-(param3)
    folder_pattern = re.compile(r'^([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)-([-+]?\d*\.?\d+)$')

    if not os.path.exists(base_dir):
        print(f"Error: Directory '{base_dir}' does not exist.")
        return

    print(f"Scanning '{base_dir}' and reading data...")
    for folder_name in sorted(os.listdir(base_dir)):
        folder_path = os.path.join(base_dir, folder_name)

        # Skip if not a directory
        if not os.path.isdir(folder_path):
            continue

        match = folder_pattern.match(folder_name)
        if not match:
            continue

        # Extract the first two varying parameters
        param1 = float(match.group(1))
        param2 = float(match.group(2))

        # Check for the data file (includes fallback for typo in "migration")
        file_path = os.path.join(folder_path, "migartion_proportion.dat")
        if not os.path.exists(file_path):
            file_path = os.path.join(folder_path, "migration_proportion.dat")

        if os.path.exists(file_path):
            try:
                with open(file_path, 'r') as f:
                    content = f.read().strip()
                    if content:
                        value = float(content)
                        phase_data.append({
                            'Param1': param1,
                            'Param2': param2,
                            'Proportion': value
                        })
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
        else:
            print(f"Warning: File missing in {folder_name}")

    if not phase_data:
        print("No valid data found. Check your folder names and paths.")
        return

    # Convert to DataFrame
    df = pd.DataFrame(phase_data)

    # Pivot table: Param2 as Y-axis, Param1 as X-axis
    pivot_table = df.pivot(index='Param2', columns='Param1', values='Proportion')

    # Sort Y descending so larger values are at the top (standard graph coordinates)
    pivot_table = pivot_table.sort_index(ascending=False)
    pivot_table = pivot_table.reindex(sorted(pivot_table.columns), axis=1)

    # --- PLOTTING ---
    print("Generating Phase Diagram...")
    plt.figure(figsize=(9, 7))

    # Heatmap setup:
    # vmin=0.0 and vmax=1.0 assume the proportion ranges from 0 to 1.
    # Remove vmin/vmax if your numbers exceed this range.
    ax = sns.heatmap(
        pivot_table,
        cmap='viridis',
        vmin=0.0,
        vmax=1.0,
        cbar_kws={'label': 'Migration Proportion'},
        square=True
    )

    # Customize titles and labels as appropriate
    plt.title('Phase Diagram (mot0.0)', pad=20, fontsize=14)
    plt.xlabel('Parameter 1', fontsize=12)
    plt.ylabel('Parameter 2', fontsize=12)

    plt.xticks(rotation=45)
    plt.yticks(rotation=0)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Plot saved successfully as '{output_file}'")
    plt.show()

if __name__ == "__main__":
    # Ensure this script is in the directory containing 'mot0.0'
    # Or pass the absolute path: plot_phase_diagram('/path/to/mot0.0')
    plot_phase_diagram()
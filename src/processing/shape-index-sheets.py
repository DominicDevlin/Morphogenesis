import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# 1. Configuration Settings
# ==========================================
# Paths to your two data files
FILE_1 = 'shape-index-data/shape_index-pre.dat'
FILE_2 = 'shape-index-data/shape_index-epiblast.dat'  # Change to your second file's name

# Labels used in the legend to distinguish the two files
LABEL_FILE_1 = 'PrE'
LABEL_FILE_2 = 'Epiblast'

WINDOW_SIZE = 1       # Rolling window size for smoothing
PLOT_RAW_DATA = True   # Set to True to see noisy raw data faintly in the background

# Styling configuration per cell type
SHAPE_CONFIG = {
    'Sox2_ShapeIndex': {'name': 'Sox2 Cells', 'color': 'pink'},
    'Sox17_ShapeIndex': {'name': 'Sox17 Cells', 'color': 'forestgreen'},
    'Loser_ShapeIndex': {'name': 'Loser Cells', 'color': 'crimson'}
}

# ==========================================
# 2. Helper Function to Process Data
# ==========================================
def load_and_process(file_path):
    print(f"Loading and processing {file_path}...")
    df = pd.read_csv(file_path, sep='\t', na_values='NA')

    columns = list(SHAPE_CONFIG.keys())
    # Group by 'Time' and compute the mean across cell types
    df_avg = df.groupby('Time')[columns].mean()

    # Apply rolling moving average
    df_smoothed = df_avg.rolling(window=WINDOW_SIZE, center=True, min_periods=1).mean()
    return df_avg, df_smoothed

# Load both datasets
df_avg_1, df_smoothed_1 = load_and_process(FILE_1)
df_avg_2, df_smoothed_2 = load_and_process(FILE_2)

# ==========================================
# 3. Plotting
# ==========================================
print("Generating plot...")
fig, ax = plt.subplots(figsize=(10, 6))

# --- Plot File 1 (Solid lines) ---
for col, style in SHAPE_CONFIG.items():
    if PLOT_RAW_DATA:
        ax.plot(df_avg_1.index, df_avg_1[col], color=style['color'], alpha=0.15)
    ax.plot(
        df_smoothed_1.index, df_smoothed_1[col],
        label=f"{style['name']} ({LABEL_FILE_1})",
        color=style['color'],
        linestyle='-',
        linewidth=2.5
    )

# --- Plot File 2 (Dashed lines) ---
for col, style in SHAPE_CONFIG.items():
    if PLOT_RAW_DATA:
        ax.plot(df_avg_2.index, df_avg_2[col], color=style['color'], alpha=0.15, linestyle=':')
    ax.plot(
        df_smoothed_2.index, df_smoothed_2[col],
        label=f"{style['name']} ({LABEL_FILE_2})",
        color=style['color'],
        # linestyle='--',
        linewidth=2.5
    )

# --- Reference Line & Formatting ---
# ax.axhline(y=3.81, color='black', linestyle=':', alpha=0.6, label='Solid-Fluid Transition (~3.81)')
ax.set_title('Average Cell Shape Index Over Time', fontsize=15, fontweight='bold')
ax.set_xlabel('Time Step', fontsize=13)
ax.set_ylabel('Shape Index (Average)', fontsize=13)
ax.legend(fontsize=10, loc='best')
ax.grid(True, linestyle=':', alpha=0.7)
ax.set_xlim(0,10000)
ax.set_ylim(3.5, 5)

plt.tight_layout()

plt.show()
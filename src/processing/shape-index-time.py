import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# 1. Configuration Settings
# ==========================================
FILE_NAME = 'shape-index-data/shape_index.dat'  # Change to your actual file name
WINDOW_SIZE = 50                      # Rolling window size for smoothing
PLOT_RAW_DATA = True                  # Set to True to see noisy raw data faintly in the background

# ==========================================
# 2. Load and Process the Data
# ==========================================
print("Loading data...")
# Read tab-separated file and treat 'NA' as NaN
df = pd.read_csv(FILE_NAME, sep='\t', na_values='NA')

print("Separating target_perim by cell type...")
# Assign target_perim to the respective cell type where its ShapeIndex is not NaN
df['Sox2_TargetPerim'] = df['target_perim'].where(df['Sox2_ShapeIndex'].notna())
df['Sox17_TargetPerim'] = df['target_perim'].where(df['Sox17_ShapeIndex'].notna())
df['Loser_TargetPerim'] = df['target_perim'].where(df['Loser_ShapeIndex'].notna())

# Columns to analyze
columns_to_average = [
    'Sox2_ShapeIndex', 'Sox17_ShapeIndex', 'Loser_ShapeIndex',
    'Sox2_TargetPerim', 'Sox17_TargetPerim', 'Loser_TargetPerim'
]

print("Calculating averages per time step...")
# Group by 'Time' and compute the mean (NaN values are ignored automatically)
df_avg = df.groupby('Time')[columns_to_average].mean()

print(f"Applying smoothing (Window Size = {WINDOW_SIZE})...")
# Rolling moving average for smoothing
df_smoothed = df_avg.rolling(window=WINDOW_SIZE, center=True, min_periods=1).mean()

# ==========================================
# 3. Plotting
# ==========================================
print("Generating plots...")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

# Styling configurations
shape_config = {
    'Sox2_ShapeIndex': {'label': 'Sox2 Cells', 'color': 'dodgerblue'},
    'Sox17_ShapeIndex': {'label': 'Sox17 Cells', 'color': 'forestgreen'},
    'Loser_ShapeIndex': {'label': 'Loser Cells', 'color': 'crimson'}
}

perim_config = {
    'Sox2_TargetPerim': {'label': 'Sox2 Target Perim', 'color': 'dodgerblue'},
    'Sox17_TargetPerim': {'label': 'Sox17 Target Perim', 'color': 'forestgreen'},
    'Loser_TargetPerim': {'label': 'Loser Target Perim', 'color': 'crimson'}
}

# --- Top Subplot: Shape Index ---
for col, style in shape_config.items():
    if PLOT_RAW_DATA:
        ax1.plot(df_avg.index, df_avg[col], color=style['color'], alpha=0.15)
    ax1.plot(df_smoothed.index, df_smoothed[col], 
             label=style['label'], color=style['color'], linewidth=2.5)

ax1.axhline(y=3.81, color='black', linestyle='--', alpha=0.5, label='Solid-Fluid Transition (~3.81)')
ax1.set_title('Average Cell Shape Index Over Time', fontsize=15, fontweight='bold')
ax1.set_ylabel('Shape Index (Average)', fontsize=13)
ax1.legend(fontsize=11, loc='best')
ax1.grid(True, linestyle=':', alpha=0.7)

# --- Bottom Subplot: Target Perimeter ---
for col, style in perim_config.items():
    if PLOT_RAW_DATA:
        ax2.plot(df_avg.index, df_avg[col], color=style['color'], alpha=0.15)
    ax2.plot(df_smoothed.index, df_smoothed[col], 
             label=style['label'], color=style['color'], linewidth=2.5)

ax2.set_title('Average Target Perimeter Over Time', fontsize=15, fontweight='bold')
ax2.set_xlabel('Time Step', fontsize=13)
ax2.set_ylabel('Target Perimeter (Average)', fontsize=13)
ax2.legend(fontsize=11, loc='best')
ax2.grid(True, linestyle=':', alpha=0.7)

plt.tight_layout()

# Save and show the plot
output_image = 'shape_index_and_perim_plot.png'
plt.savefig(output_image, dpi=300)
print(f"Plot saved successfully as '{output_image}'")

plt.show()
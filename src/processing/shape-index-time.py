import pandas as pd
import matplotlib.pyplot as plt

# ==========================================
# 1. Configuration Settings
# ==========================================
FILE_NAME = 'shape-index-data/shape_index.dat'  # Change to your actual file name
WINDOW_SIZE = 50                      # Increase this number for MORE smoothing, decrease for LESS
PLOT_RAW_DATA = True                  # Set to True to see the noisy data faintly in the background

# ==========================================
# 2. Load and Process the Data
# ==========================================
print("Loading data...")
# Read the file. 
# sep='\t' tells pandas the file is tab-separated.
# na_values='NA' automatically converts the 'NA' strings into mathematical NaNs (Not a Number).
df = pd.read_csv(FILE_NAME, sep='\t', na_values='NA')

print("Calculating averages per time step...")
# Group by 'Time' and calculate the mean for each column.
# Pandas automatically ignores 'NaN' values when calculating the mean, 
# which perfectly handles our staggered columns!
df_avg = df.groupby('Time')[['Sox2_ShapeIndex', 'Sox17_ShapeIndex', 'Loser_ShapeIndex']].mean()

print(f"Applying smoothing (Window Size = {WINDOW_SIZE})...")
# Apply a rolling moving average to massively smooth the data.
# center=True prevents the smoothed line from "lagging" behind the real time.
# min_periods=1 ensures we still get lines at the very beginning/end of the simulation.
df_smoothed = df_avg.rolling(window=WINDOW_SIZE, center=True, min_periods=1).mean()

# ==========================================
# 3. Plotting
# ==========================================
print("Generating plot...")
plt.figure(figsize=(10, 6))

# Colors for the different cell types
colors = {
    'Sox2_ShapeIndex': 'dodgerblue',
    'Sox17_ShapeIndex': 'forestgreen',
    'Loser_ShapeIndex': 'crimson'
}

labels = {
    'Sox2_ShapeIndex': 'Sox2 Cells',
    'Sox17_ShapeIndex': 'Sox17 Cells',
    'Loser_ShapeIndex': 'Loser Cells'
}

# Plot the lines
for column in df_smoothed.columns:
    # Optional: Plot the raw, noisy data in the background with low opacity
    if PLOT_RAW_DATA:
        plt.plot(df_avg.index, df_avg[column], color=colors[column], alpha=0.15)
    
    # Plot the smoothed data heavily
    plt.plot(df_smoothed.index, df_smoothed[column], 
             label=labels[column], color=colors[column], linewidth=2.5)

# Formatting the plot
plt.title('Average Cell Shape Index Over Time', fontsize=16, fontweight='bold')
plt.xlabel('Time Step', fontsize=14)
plt.ylabel('Shape Index (Average)', fontsize=14)

# Add a reference line for standard solid/fluid transition (Shape Index ~ 3.81 in 2D tissue)
# Feel free to remove this if it doesn't apply to your specific model
plt.axhline(y=3.81, color='black', linestyle='--', alpha=0.5, label='Solid-Fluid Transition (~3.81)')

plt.legend(fontsize=12, loc='best')
plt.grid(True, linestyle=':', alpha=0.7)
plt.tight_layout()

# Save and show the plot
output_image = 'shape_index_plot.png'
plt.savefig(output_image, dpi=300)
print(f"Plot saved successfully as '{output_image}'")

plt.show()
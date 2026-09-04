import os
import matplotlib.pyplot as plt
import pandas as pd

# ==========================================
# 1. Configuration Settings
# ==========================================
# Set directory path (easy to change parameters here)
PARAM_FOLDER = "0.5-0.1-0.0"
DATA_DIR = os.path.join("org-data", PARAM_FOLDER)

WINDOW_SIZE = 5  # Rolling window size for smoothing (adjust based on your total steps)
PLOT_INDIVIDUAL_SIMS = (
    True  # Set to True to see faded lines for each individual simulation
)

# File names and display styling
DATA_CONFIG = {
    "Sox2": {
        "file": "average_shape_indices_sox2.dat",
        "color": "dodgerblue",
        "label": "Sox2 Cells",
    },
    "Sox17": {
        "file": "average_shape_indices_sox17.dat",
        "color": "forestgreen",
        "label": "Sox17 Cells",
    },
    "Loser": {
        "file": "average_shape_indices_loser.dat",
        "color": "crimson",  # Red line
        "label": "Loser Cells",
    },
}

# ==========================================
# 2. Plotting
# ==========================================
plt.figure(figsize=(10, 6))

for cell_type, config in DATA_CONFIG.items():
    file_path = os.path.join(DATA_DIR, config["file"])

    if not os.path.exists(file_path):
        print(f"Warning: File not found: {file_path}. Skipping...")
        continue

    print(f"Loading {cell_type} data from {file_path}...")
    df = pd.read_csv(file_path, sep="\t", na_values="NA")

    time = df["time"]
    sim_cols = [col for col in df.columns if col.startswith("sim_")]

    # 1. Plot individual simulations (faded in the background)
    if PLOT_INDIVIDUAL_SIMS:
        for col in sim_cols:
            plt.plot(
                time,
                df[col],
                color=config["color"],
                alpha=0.15,
                linewidth=1,
                zorder=1,
            )

    # 2. Calculate average across simulations per time step
    mean_across_sims = df[sim_cols].mean(axis=1)

    # 3. Smooth the averaged line
    smoothed_mean = mean_across_sims.rolling(
        window=WINDOW_SIZE, center=True, min_periods=1
    ).mean()

    # 4. Plot the smoothed average line prominently
    plt.plot(
        time,
        smoothed_mean,
        color=config["color"],
        linewidth=2.5,
        label=config["label"],
        zorder=3,
    )

# Reference line for solid-fluid jamming transition
plt.axhline(
    y=3.81,
    color="black",
    linestyle="--",
    alpha=0.6,
    label="Solid-Fluid Transition (~3.81)",
    zorder=2,
)

# Plot formatting
plt.title(
    f"Cell Shape Index Over Time ({PARAM_FOLDER})",
    fontsize=14,
    fontweight="bold",
)
plt.xlabel("Monte Carlo Steps (MCS)", fontsize=12)
plt.ylabel("Average Shape Index ($P / \\sqrt{A}$)", fontsize=12)
plt.grid(True, linestyle=":", alpha=0.6)
plt.legend(fontsize=11, loc="best")
plt.tight_layout()

# Save output inside the same directory
output_image = os.path.join(DATA_DIR, "shape_index_by_cell_type.png")
# plt.savefig(output_image, dpi=300)
print(f"Plot saved successfully as '{output_image}'")

plt.show()
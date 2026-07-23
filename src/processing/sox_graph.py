import io
import matplotlib.pyplot as plt
import numpy as np



# Load data (using StringIO to simulate reading from a file)
data = np.loadtxt("sox_start_values.dat")
x = data[:, 0]
y = data[:, 1]

# Define the threshold
threshold = 0.2

# Define masks for the three categories
# Category 1: One above and one below the threshold
mask_mixed = ((x > threshold) & (y < threshold)) | (
    (x < threshold) & (y > threshold)
)

# Category 2: Both above the threshold
mask_both_high = (x >= threshold) & (y >= threshold)

# Category 3: Both below the threshold
mask_both_low = (x <= threshold) & (y <= threshold)

print(mask_mixed)

# Set up the plot
plt.figure(figsize=(8, 6))

# Plot each category with a different color
plt.scatter(
    x[mask_mixed],
    y[mask_mixed],
    color="orange",
    label="One above, one below",
    s=100,
    zorder=3,
)
plt.scatter(
    x[mask_both_high],
    y[mask_both_high],
    color="blue",
    label="Both above 0.2",
    s=100,
    zorder=3,
)
plt.scatter(
    x[mask_both_low],
    y[mask_both_low],
    color="gray",
    label="Both below 0.2",
    s=100,
    zorder=3,
)

# Add dashed lines to indicate the threshold boundaries
plt.axhline(
    y=threshold,
    color="red",
    linestyle="--",
    alpha=0.5,
    label=f"Threshold ({threshold})",
)
plt.axvline(x=threshold, color="red", linestyle="--", alpha=0.5)

# Labels and styling
plt.xlabel("X values")
plt.ylabel("Y values")
# Display the plot
plt.tight_layout()
plt.show()
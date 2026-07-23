import matplotlib.pyplot as plt
from collections import defaultdict

# Define the filename
filename = 'death_total.dat'

# Dictionary to hold our data. 
# Keys will be the 1st column (e.g., 5, 8, 12, 16)
# Values will be a list of the 2nd column values in order of appearance
data_groups = defaultdict(list)

# Read the .dat file
with open(filename, 'r') as file:
    for line in file:
        # Skip empty lines
        if not line.strip():
            continue

        # Split the line by whitespace (spaces or tabs)
        col1, col2 = line.split()
        
        # Convert to appropriate data types
        group_id = int(col1)
        value = float(col2)
        
        # Append the value to the specific group
        data_groups[group_id].append(value)

# Create the plot
plt.figure(figsize=(10, 6))

# Sort the keys so the legend is in numerical order (5, 8, 12, 16)
for group_id in sorted(data_groups.keys()):
    y_values = data_groups[group_id]
    
    # Generate the x-axis (time) based on the number of items in the list
    # This creates a sequence: 0, 1, 2, 3...
    x_time = range(len(y_values))
    
    # Plot the line for this group
    plt.plot(x_time, y_values, marker='o', label=f'Group {group_id}')

# Formatting the plot
plt.title('Data Grouped by Column 1 over Time')
plt.xlabel('Time (Sequential Order)')
plt.ylabel('Value (Column 2)')
# plt.xticks(range(max(len(v) for v in data_groups.values()))) # Ensure x-axis shows whole numbers
plt.legend()
# plt.grid(True, linestyle='--', alpha=0.7)
plt.xlim(0,100)

# Display the plot
plt.tight_layout()
plt.show()
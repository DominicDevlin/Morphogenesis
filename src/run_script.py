import os
import sys

def main():
    # Verify that an index argument was provided
    if len(sys.argv) < 2:
        print("Usage: python run_sweep.py <index>")
        sys.exit(1)

    try:
        index = int(sys.argv[1])
    except ValueError:
        print("Error: The index argument must be an integer.")
        sys.exit(1)

    # Define the parameter arrays
    adhesion_values = [-0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    apop_values = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 30]

    n_cols = len(adhesion_values)
    num_combinations = len(adhesion_values) * len(apop_values)

    # Ensure the index is within the valid range of combinations
    if index < 0 or index >= num_combinations:
        print(f"Error: Index {index} is out of range. Must be between 0 and {num_combinations - 1}.")
        sys.exit(1)

    # Determine J_L and J_S based on the index
    row = index // n_cols  # Row index for J_S
    col = index % n_cols   # Column index for J_L

    a_val = adhesion_values[col]
    apop_val = apop_values[row]

    print(f"Index: {index} => Selected adhesino val: {a_val}, apop val: {apop_val}")

    # Build and run the command to execute the C++ file
    # Adjust "./phase-optimize" to match your actual executable if needed
    command = f"./embryo_multi {apop_val} {a_val}"
    print(f"Running command: {command}")
    
    exit_status = os.system(command)
    
    if exit_status != 0:
        print(f"Warning: Executable exited with non-zero status code: {exit_status}")
        sys.exit(exit_status)

if __name__ == "__main__":
    main()
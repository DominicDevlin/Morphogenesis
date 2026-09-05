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
    adhesion_values = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8]
    perim_values = [0., 0.03, 0.06, 0.09, 0.12, 0.15, 0.18, 0.21, 0.24, 0.27]
    motility_values = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45]

    n_cols = len(adhesion_values)
    num_combinations = len(adhesion_values) * len(perim_values)

    # Ensure the index is within the valid range of combinations
    if index < 0 or index >= num_combinations:
        print(f"Error: Index {index} is out of range. Must be between 0 and {num_combinations - 1}.")
        sys.exit(1)

    row = index // n_cols  
    col = index % n_cols  

    a_val = adhesion_values[col]
    # perim_val = perim_values[row]
    perim_val = 0.15
    # motility=0.2
    motility = motility_values[row]

    print(f"Index: {index} => Selected adhesion val: {a_val}, added perim proportion val: {perim_val}, loser motility: {motility}")

    # Build and run the command to execute the C++ file
    # Adjust "./phase-optimize" to match your actual executable if needed
    command = f"./embryo_multi {a_val} {perim_val} {motility}"
    print(f"Running command: {command}")
    
    exit_status = os.system(command)
    
    if exit_status != 0:
        print(f"Warning: Executable exited with non-zero status code: {exit_status}")
        sys.exit(exit_status)

if __name__ == "__main__":
    main()
import pandas as pd
import glob
import os
import shutil
import numpy as np

def clean_high_variance_directories(root_dir="t1/", variance_threshold=30.0):
    """
    Iterates through subdirectories, calculates variance of the final 'total' 
    cells across celltypes-org replicates, and marks the directory for deletion if above threshold.
    """
    
    # Check if root directory existss
    if not os.path.exists(root_dir):
        print(f"Root directory '{root_dir}' not found.")
        return

    # Get all subdirectories within the root directory
    subdirs = [os.path.join(root_dir, d) for d in os.listdir(root_dir) 
               if os.path.isdir(os.path.join(root_dir, d))]
    
    print(f"Found {len(subdirs)} subdirectories in '{root_dir}'.\n")

    dirs_to_delete = []

    for subdir in subdirs:
        # Specifically match files starting with 'celltypes-org-'
        file_pattern = os.path.join(subdir, "celltypes-org-*.dat")
        files = glob.glob(file_pattern)
        
        if not files:
            print(f"Skipping {subdir}: No 'celltypes-org-*.dat' files found.")
            continue
            
        if len(files) < 2:
            print(f"Skipping {subdir}: Not enough replicates to calculate variance (Found {len(files)}).")
            continue

        # Extract the final 'total' value from each replicate
        replicate_final_totals = []
        
        for file in files:
            file_name = os.path.basename(file)
            try:
                df = pd.read_csv(file, sep='\s+') 
                if 'total' in df.columns:
                    # Get the very last 'total' value from the time series
                    final_total = df['total'].iloc[-1]
                    replicate_final_totals.append(final_total)
                else:
                    print(f"Warning: 'total' column not found in {file_name}")
            except Exception as e:
                print(f"Error reading {file_name}: {e}")

        # Calculate variance across the replicates
        if len(replicate_final_totals) > 1:
            # ddof=1 calculates sample variance
            variance = np.var(replicate_final_totals, ddof=1)
            
            if variance > variance_threshold:
                print(f"[MARK FOR DELETION] {subdir} | Variance: {variance:.2f}")
                dirs_to_delete.append(subdir)
            else:
                print(f"[SAFE] {subdir} | Variance: {variance:.2f}")

    # --- DELETION SUMMARY & EXECUTION ---
    print("\n" + "="*50)
    print(f"SUMMARY: {len(dirs_to_delete)} out of {len(subdirs)} directories marked for deletion.")
    print("="*50)
    
    for d in dirs_to_delete:
        print(f"Pending deletion: {d}")
        
        # =========================================================
        # UNCOMMENT THE LINES BELOW TO ACTUALLY DELETE DIRECTORIES
        # =========================================================
        # try:
        #     shutil.rmtree(d)
        #     print(f"Successfully deleted {d}")
        # except Exception as e:
        #     print(f"Failed to delete {d}. Reason: {e}")

if __name__ == "__main__":
    # Adjust your variance threshold here
    clean_high_variance_directories()
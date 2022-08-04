#!/bin/bash -e
#SBATCH --job-name=CPM_ev_40  
#SBATCH --time=23:00:00      # Walltime (HH:MM:SS)
#SBATCH --mem=4GB 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --account=nesi03638         
#SBATCH --output=ev_sim_out-%j.out 
#SBATCH --error=ev_sim_err-%j.out 

module load LegacySystemLibs/7
./evolution

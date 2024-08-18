#!/bin/bash -e
#SBATCH --job-name=CPM_evolution  
#SBATCH --time=00:20:00      # Walltime (HH:MM:SS)
#SBATCH --mem=150GB 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30  
#SBATCH --account=uoa02799         
#SBATCH --output=ev_sim_out-%j.out 
#SBATCH --error=ev_sim_err-%j.out 



find . -name "*.o" -type f -delete
module load Qt5/5.12.3-GCCcore-9.2.0 
qmake
make
module load LegacySystemLibs/7
xvfb-run -a ./phasemulti

## to output images on the cluster, prepend the output with "xvfb-run". e.g. "xvfb-run ./evolution"

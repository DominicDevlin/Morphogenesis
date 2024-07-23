#!/bin/bash -e
#SBATCH --job-name=CPM_evolution  
#SBATCH --time=70:00:00      # Walltime (HH:MM:SS)
#SBATCH --mem=2GB 
#SBATCH --array=0-8
#SBATCH --cpus-per-task=60  
#SBATCH --account=uoa02799         
#SBATCH --output=sim_out-%j.out 
#SBATCH --error=sim_err-%j.out 
#SBATCH --partition=milan



find . -name "*.o" -type f -delete
module load Qt5/5.12.3-GCCcore-9.2.0 
if [ ! -f build_done.flag ]; then
    touch build_done.flag
    qmake
    make
else
    sleep 60
fi
module load LegacySystemLibs/7
xvfb-run -a ./phase-evolution

## to output images on the cluster, prepend the output with "xvfb-run -a". e.g. "xvfb-run -a ./evolution"

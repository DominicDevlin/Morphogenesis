#!/bin/bash -e
#SBATCH --job-name=CPM_evolution
#SBATCH --time=12:00:00      # Walltime (HH:MM:SS)
#SBATCH --mem=2GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0-11
#SBATCH --account=uoa02799
#SBATCH --output=ev_sim_out-%j.out
#SBATCH --error=ev_sim_err-%j.out
#SBATCH --partition=milan



find . -name "*.o" -type f -delete
module load Qt5/5.12.3-GCCcore-9.2.0
module load Python/3.9.5-gimkl-2020a
qmake
make
module load LegacySystemLibs/7
python3 optimize-script.py ${SLURM_ARRAY_TASK_ID}

## to output images on the cluster, prepend the output with "xvfb-run". e.g. "xvfb-run ./evolution"

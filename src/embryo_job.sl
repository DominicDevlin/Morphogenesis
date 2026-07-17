#!/bin/bash -e
#SBATCH --job-name=CPM_evolution  
#SBATCH --time=1:00:00      # Walltime (HH:MM:SS)
#SBATCH --mem=10GB 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --array=0-80  
#SBATCH --account=uoa02799         
#SBATCH --output=ev_sim_out-%j.out 
#SBATCH --error=ev_sim_err-%j.out 

module load GCCcore/11.3.0          #  or GCCcore/11.2.0
module load binutils/2.38-GCCcore-11.3.0   # linker that matches GCC 11


if [ ! -f build_done.flag ]; then
    touch build_done.flag
    which qmake-qt5 && qmake-qt5 || qmake            # or `qmake-qt5` if that exists on your system
    make  -j $SLURM_CPUS_ON_NODE
else
    sleep 40
fi

python3 run_script.py ${SLURM_ARRAY_TASK_ID}


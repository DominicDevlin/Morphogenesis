#!/bin/bash -e
#SBATCH --job-name=CPM_evolution  
#SBATCH --time=80:00:00      # Walltime (HH:MM:SS)
#SBATCH --mem=2GB 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30  
#SBATCH --account=uoa02799         
#SBATCH --output=ev_sim_out-%j.out 
#SBATCH --error=ev_sim_err-%j.out 


module load GCCcore/11.3.0          #  or GCCcore/11.2.0
module load binutils/2.38-GCCcore-11.3.0   # linker that matches GCC 11

#––– regenerate Makefile –––
which qmake-qt5 && qmake-qt5 || qmake            # or `qmake-qt5` if that exists on your system

#––– now the Makefile exists, so 'distclean' is defined –––
make distclean             # optional – only if you want to wipe leftovers

which qmake-qt5 && qmake-qt5 || qmake            # or `qmake-qt5` if that exists on your system

# full build
make  -j $SLURM_CPUS_ON_NODE
./evolution
## to output images on the cluster, prepend the output with "xvfb-run". e.g. "xvfb-run ./evolution"



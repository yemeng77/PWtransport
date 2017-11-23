#!/bin/bash -l
#SBATCH -N 2         #Use 2 nodes
#SBATCH -t 00:30:00  #Set 30 minute time limit
#SBATCH -p debug   #Submit to the regular 'partition'
#SBATCH -L SCRATCH   #Job requires $SCRATCH file system
#SBATCH -C haswell   #Use Haswell nodes
#SBATCH -J Cu-CENTRAL

srun -n 64 -c 2 ~/cori/bin/PEtot_trans1

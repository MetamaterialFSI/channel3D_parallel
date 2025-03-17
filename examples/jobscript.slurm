#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=2:00:00

#SBATCH --job-name=channel_trial
#SBATCH --output=status.out

# module load openmpi
module load intel-oneapi-mpi

input_file="input_parameters"
mkdir channel_output

cd $SLURM_SUBMIT_DIR
mpirun -np 2 ../build/channel3D "$input_file" 

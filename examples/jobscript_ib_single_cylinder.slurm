#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00

#SBATCH --job-name=channel_trial
#SBATCH --output=status_ib_single_cylinder.out

# module load openmpi
module load intel-oneapi-mpi

input_file="input_parameters_ib_single_cylinder"

mkdir -p channel_output_ib_single_cylinder

cd $SLURM_SUBMIT_DIR
mpirun -np 2 ../build/channel3D "$input_file" 

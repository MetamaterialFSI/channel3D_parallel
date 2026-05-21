#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00

#SBATCH --job-name=channel_trial
#SBATCH --output=status_ib_bottom_wall.out

module load openmpi

input_file="input_parameters_ib_bottom_wall"

cd $SLURM_SUBMIT_DIR
mpirun -np 2 ../build/channel3D "$input_file" 

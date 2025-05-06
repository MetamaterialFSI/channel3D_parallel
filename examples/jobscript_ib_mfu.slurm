#!/bin/bash

#SBATCH --nodes=1
#SBATCH -q debug
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00

#SBATCH --job-name=channel_trial
#SBATCH --output=status_ib_mfu.out

module load openmpi

input_file="input_parameters_ib_mfu"

mkdir -p channel_output_ib_mfu

cd $SLURM_SUBMIT_DIR
mpirun -np 2 ../build/channel3D "$input_file" 

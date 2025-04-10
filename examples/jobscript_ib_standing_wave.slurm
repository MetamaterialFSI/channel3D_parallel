#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00

#SBATCH --job-name=channel_trial
#SBATCH --output=status_standing_wave.out

module load openmpi

input_file="input_parameters_ib_standing_wave"

mkdir -p channel_output_ib_standing_wave

cd $SLURM_SUBMIT_DIR
mpirun -np 2 ../build/channel3D "$input_file" 

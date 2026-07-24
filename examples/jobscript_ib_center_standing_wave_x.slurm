#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

#SBATCH --job-name=center_standing_wave_x
#SBATCH --output=status_ib_center_standing_wave_x.out

# module load openmpi
module load intel-oneapi-mpi

input_file="input_parameters_ib_center_standing_wave_x"

mkdir -p channel_output_ib_center_standing_wave_x

cd $SLURM_SUBMIT_DIR
mpirun -np 2 ../build/channel3D "$input_file" 

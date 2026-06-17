#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00

#SBATCH --job-name=double_cylinders_x
#SBATCH --output=status_ib_double_rotating_cylinders_x.out

module load openmpi

input_file="input_parameters_ib_double_rotating_cylinders_x"

mkdir -p channel_output_ib_double_rotating_cylinders_x

cd $SLURM_SUBMIT_DIR
mpirun -np 2 ../build/channel3D "$input_file"

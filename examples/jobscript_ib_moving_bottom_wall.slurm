#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH -q normal
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=ctlin@caltech.edu
#SBATCH --job-name=channel_moving_bottom_wall
#SBATCH --output=status_ib_moving_bottom_wall_y+15_coarse_2ch.out

module load openmpi

input_file="input_parameters_ib_moving_bottom_wall"

cd $SLURM_SUBMIT_DIR
mpirun -np 32 ../build/channel3D "$input_file" 
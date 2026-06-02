#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH -q normal
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH --mail-user=ctlin@caltech.edu
#SBATCH --job-name=channel_IB_botom_wall
#SBATCH --output=status_ib_bottom_wall_test.out

module load openmpi

input_file="input_parameters_ib_bottom_wall"

cd $SLURM_SUBMIT_DIR
mpirun -np 64 ../build/channel3D "$input_file" 
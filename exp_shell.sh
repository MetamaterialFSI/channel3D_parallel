module load fftw
module load intel-oneapi-compilers
module load openmpi
rm -rf ./build
cmake -B build
cmake --build build
cd ./examples
sbatch jobscript_ib_planar_wall.slurm 


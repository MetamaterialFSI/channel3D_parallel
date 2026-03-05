module load fftw
module load intel-oneapi-compilers
module load openmpi
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build 
cd ./examples/
#rm -rf ./channel_output_ib_planar_wall/
sbatch jobscript_ib_planar_wall.slurm

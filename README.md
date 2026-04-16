# channel3D

[Documentation](https://metamaterialfsi.github.io/channel3D_parallel/)

## Building from source

### Prerequisites
- CMake
- LAPACK
- MPI (tested with version 3.1)
- FFTW with MPI enabled (tested with version 3.3.10)
- Intel Fortran compiler (extension to non-Intel compilers is in progress)

For example, on the Caltech HPC, you can run
```
module load fftw
module load intel-oneapi-compilers
module load openmpi
```
or, if you want to use Intel MPI,
```
module load fftw/3.3.10-oneapi-2023.2.1-7czoymn
module load intel-oneapi-compilers
module load intel-oneapi-mpi
```

### Building with CMake

```
cmake -B build
cmake --build build
```

or, similarly,

```
mkdir build
cd build
cmake ..
make
```

By default, CMake builds the `Release` version, which optimizes the Fortran code using the flags `-ipo -O3`. Alternatively, you can run CMake with `-DCMAKE_BUILD_TYPE=Debug` to enable backtraces and extra debugging information.

## Running tests

To run the built-in tests, configure CMake to build the test framework:

```

cmake -B build -DBUILD_TESTING=ON
cmake --build build
ctest --test-dir build
```

or, similarly,

```
mkdir build
cd build
cmake .. -DBUILD_TESTING=ON
make
ctest
```

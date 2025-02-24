# channel3D

## Building from source

### TODO: Prerequisites

#### TODO: Fortran compiler

#### TODO: Math libraries

### Manual compilation

```
export F_UFMTENDIAN=big
mpif90 -O3 -ipo mpi.f90 global.f90 mass_flow.f90 interpolation.f90 equations.f90 boundary_conditions.f90 pressure.f90 input_output.f90 statistics.f90 initialization.f90 finalization.f90 projection.f90 time_integration.f90 monitor.f90 main.f90 -o channel_oc -lfftw3_mpi -lfftw3 PATH_TO_LIBLAPACK -I PATH_TO_FFTW_INCLUDE -L PATH_TO_FFTW_LIB
```

or, using Intel oneAPI MKL,

```
export F_UFMTENDIAN=big
mpifort -O3 -ipo mpi.f90 global.f90 mass_flow.f90 interpolation.f90 equations.f90 boundary_conditions.f90 pressure.f90 input_output.f90 statistics.f90 initialization.f90 finalization.f90 projection.f90 time_integration.f90 monitor.f90 main.f90 -o channel_oc -lfftw3_mpi PATH_TO_LIBLAPACK -I PATH_TO_FFTW_INCLUDE -L PATH_TO_FFTW_LIB  -qmkl=parallel
```

### Building with CMake

Make sure that the Intel Fortran compiler (ifort or ifx), Math Kernel Library (MKL), and a compatible MPI implementation are installed and available in your system's PATH. If not, you might need to load it first by sourcing the Intel environment setup script `/path/to/intel/oneapi/setvars.sh`.

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

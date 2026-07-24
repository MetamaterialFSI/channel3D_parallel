# channel3D

[Documentation](https://metamaterialfsi.github.io/channel3D_parallel/)

## Building from source

### Prerequisites
- CMake
- LAPACK
- MPI (tested with version 3.1)
- FFTW with MPI enabled (tested with version 3.3.10)
- A Fortran compiler: Intel (`ifx`, `ifort`) or GNU (`gfortran`)

Note that FFTW must be built with MPI support; a plain FFTW install is not
enough. CMake reports this at configure time as
`Could NOT find FFTW (missing: DOUBLE_MPI_LIB)`.

#### Caltech HPC

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

#### Ubuntu

```
sudo apt-get install -y gfortran cmake \
    libopenmpi-dev openmpi-bin \
    libfftw3-dev libfftw3-mpi-dev \
    libblas-dev liblapack-dev
```

#### macOS

```
brew install gcc cmake open-mpi fftw
```

LAPACK comes from the Accelerate framework, so it does not need to be
installed. Homebrew's `fftw` already includes the MPI library, and `gcc`
provides `gfortran`.

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

By default, CMake builds the `Release` version, which compiles with `-O3` and enables interprocedural optimization (`-ipo` for Intel, `-flto` for GNU). Alternatively, you can run CMake with `-DCMAKE_BUILD_TYPE=Debug` to enable backtraces and extra debugging information.

To pick a specific compiler, set `FC` when configuring, e.g. `FC=gfortran cmake -B build`.

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

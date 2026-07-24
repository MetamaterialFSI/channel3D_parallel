# Installation

This guide walks you through building the project from source, including required dependencies and build steps.

---

## Prerequisites

Before building, make sure the following dependencies are installed on your system:

- **CMake**
- **LAPACK**
- **MPI** (tested with version 3.1)
- **FFTW with MPI support** (tested with version 3.3.10)
- **Fortran compiler** — Intel (`ifx`, `ifort`) or GNU (`gfortran`)

!!! warning "FFTW must have MPI support"
    A plain FFTW installation is not sufficient. If the MPI-enabled library is
    missing, CMake stops at configure time with
    `Could NOT find FFTW (missing: DOUBLE_MPI_LIB)`.

---

## Environment Setup

### Caltech HPC

On Caltech HPC, you can load the required modules with:

```bash
module load fftw
module load intel-oneapi-compilers
module load openmpi
```

Alternatively, if you prefer Intel MPI on Caltech HPC:

```bash
module load fftw/3.3.10-oneapi-2023.2.1-7czoymn
module load intel-oneapi-compilers
module load intel-oneapi-mpi
```

### Ubuntu

```bash
sudo apt-get install -y gfortran cmake \
    libopenmpi-dev openmpi-bin \
    libfftw3-dev libfftw3-mpi-dev \
    libblas-dev liblapack-dev
```

`libfftw3-mpi-dev` is the MPI-enabled FFTW, built against Open MPI.

### macOS

```bash
brew install gcc cmake open-mpi fftw
```

LAPACK is provided by the Accelerate framework and needs no package.
Homebrew's `fftw` already ships the MPI library, and `gcc` provides
`gfortran`.

---

## Choosing a Compiler

CMake picks up the compiler from the `FC` environment variable:

```bash
FC=gfortran cmake -B build
```

Both compiler families are exercised by CI: Intel `ifx` on Linux, and
`gfortran` on both Ubuntu and macOS.

---

## Building the Project

```bash
cmake -B build
cmake --build build
```

Alternatively, you can use
```bash
mkdir build
cd build
cmake ..
make
```

---

## Build Configuration

- **Default build type:** `Release`  
  - Compiles with `-O3` and enables interprocedural optimization
    (`-ipo` for Intel, `-flto` for GNU)

- **Debug build:**  
  To enable debugging symbols and backtraces:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
```

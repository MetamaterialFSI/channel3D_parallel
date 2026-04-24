# Installation

This guide walks you through building the project from source, including required dependencies and build steps.

---

## Prerequisites

Before building, make sure the following dependencies are installed on your system:

- **CMake**
- **LAPACK**
- **MPI** (tested with version 3.1)
- **FFTW with MPI support** (tested with version 3.3.10)
- **Fortran compiler** (currently, only Intel Fortran compilers are supported)

---

## Environment Setup (Example for Caltech HPC)

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
  - Applies optimization flags: `-ipo -O3`

- **Debug build:**  
  To enable debugging symbols and backtraces:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
```

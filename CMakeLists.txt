cmake_minimum_required(VERSION 3.10)

project(Channel3D VERSION 0.1
	DESCRIPTION "3D channel flow solver"
        LANGUAGES Fortran)

# set(CMAKE_Fortran_STANDARD 2018)
# set(CMAKE_Fortran_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Default to Release for single-config generators when nothing was requested.
# Multi-config generators (Ninja Multi-Config, Xcode) leave CMAKE_BUILD_TYPE empty
# and choose per build instead, so only set a default when neither is present.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Per-config optimization comes from CMake's own defaults (CMAKE_Fortran_FLAGS_RELEASE
# is -O3, CMAKE_Fortran_FLAGS_DEBUG is -g); we only add interprocedural optimization
# to the Release configuration on top of that.
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)

# Compiler-specific flags. These append to the flag variables rather than replacing
# them, so anything passed via -DCMAKE_Fortran_FLAGS=... or the FFLAGS environment
# variable is preserved.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # Some source lines exceed the 132-column free-form limit; Intel accepts them but
    # gfortran truncates at 132 (an error on some gfortran builds), so lift the limit.
    string(APPEND CMAKE_Fortran_FLAGS " -ffree-line-length-none")
    string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -fbacktrace -fcheck=bounds")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")   # matches Intel and IntelLLVM
    string(APPEND CMAKE_Fortran_FLAGS_DEBUG " -traceback -check bounds")
endif()

# ifx emits LLVM bitcode into the static library when IPO is enabled, which GNU ld
# cannot consume ("error adding symbols"); lld can. Other compilers do not need this,
# and lld is often absent (e.g. a default macOS toolchain), so only ask for it here.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    add_link_options("-fuse-ld=lld")
endif()
find_package(LAPACK REQUIRED)
# channel_lib links FFTW::Double and FFTW::DoubleMPI, so require both here:
# without the components an FFTW lacking the MPI library still configures and
# only fails later, when those targets turn out not to exist.
find_package(FFTW REQUIRED COMPONENTS DOUBLE_LIB DOUBLE_MPI_LIB)
find_package(MPI COMPONENTS Fortran REQUIRED)

message(STATUS "CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
message(STATUS "Fortran_COMPILER_NAME): ${Fortran_COMPILER_NAME}")
message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "FFTW_INCLUDE_DIRS: ${FFTW_INCLUDE_DIRS}")

#message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input

# Include subdirectories
add_subdirectory(src)

add_executable(channel3D app/main.f90)
target_link_libraries(channel3D PRIVATE channel_lib)

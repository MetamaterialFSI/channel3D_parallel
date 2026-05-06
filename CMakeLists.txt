cmake_minimum_required(VERSION 3.10)

if(NOT DEFINED CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER mpiifx)
endif()

project(Channel3D VERSION 0.1
        DESCRIPTION "3D channel flow solver"
        LANGUAGES Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)

if(uppercase_CMAKE_BUILD_TYPE STREQUAL "RELEASE")
    message(STATUS "Build type: Release")
    set(CMAKE_Fortran_FLAGS "-O3")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)

elseif(uppercase_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    message(STATUS "Build type: Debug")

    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR
       CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        set(CMAKE_Fortran_FLAGS "-g -O0 -traceback -check bounds")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_Fortran_FLAGS "-g -O0 -fbacktrace -fcheck=bounds -ffree-line-length-none")
    endif()

else()
    message(FATAL_ERROR "Invalid value for CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()

find_package(LAPACK REQUIRED)
find_package(MPI COMPONENTS Fortran REQUIRED)

set(FFTW_ROOT "$ENV{HOME}/fftw-mpi-intel")
set(FFTW_INCLUDE_DIRS "${FFTW_ROOT}/include")

if(EXISTS "${FFTW_ROOT}/lib64")
    set(FFTW_LIB_DIR "${FFTW_ROOT}/lib64")
else()
    set(FFTW_LIB_DIR "${FFTW_ROOT}/lib")
endif()

if(EXISTS "${FFTW_LIB_DIR}/libfftw3_mpi.so")
    set(FFTW_MPI_LIB "${FFTW_LIB_DIR}/libfftw3_mpi.so")
elseif(EXISTS "${FFTW_LIB_DIR}/libfftw3_mpi.a")
    set(FFTW_MPI_LIB "${FFTW_LIB_DIR}/libfftw3_mpi.a")
else()
    message(FATAL_ERROR "Could not find libfftw3_mpi in ${FFTW_LIB_DIR}")
endif()

if(EXISTS "${FFTW_LIB_DIR}/libfftw3.so")
    set(FFTW_LIB "${FFTW_LIB_DIR}/libfftw3.so")
elseif(EXISTS "${FFTW_LIB_DIR}/libfftw3.a")
    set(FFTW_LIB "${FFTW_LIB_DIR}/libfftw3.a")
else()
    message(FATAL_ERROR "Could not find libfftw3 in ${FFTW_LIB_DIR}")
endif()

set(FFTW_LIBRARIES
    "${FFTW_MPI_LIB}"
    "${FFTW_LIB}")

message(STATUS "CMAKE_Fortran_COMPILER_ID: ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "CMAKE_Fortran_COMPILER: ${CMAKE_Fortran_COMPILER}")
message(STATUS "MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "FFTW_ROOT: ${FFTW_ROOT}")
message(STATUS "FFTW_INCLUDE_DIRS: ${FFTW_INCLUDE_DIRS}")
message(STATUS "FFTW_LIB_DIR: ${FFTW_LIB_DIR}")
message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")

add_subdirectory(src)

add_executable(channel3D app/main.f90)
target_link_libraries(channel3D PRIVATE channel_lib)

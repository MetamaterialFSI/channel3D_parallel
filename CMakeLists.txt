cmake_minimum_required(VERSION 3.27)

project(Channel3D VERSION 0.1
	DESCRIPTION "3D channel flow solver"
        LANGUAGES Fortran)

# set(CMAKE_Fortran_STANDARD 2018)
# set(CMAKE_Fortran_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE) 
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)

if(uppercase_CMAKE_BUILD_TYPE STREQUAL "RELEASE")
    message(STATUS "Build type: Release")
    set(CMAKE_Fortran_FLAGS "-O3")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
elseif(uppercase_CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    message(STATUS "Build type: Debug")
    set(CMAKE_Fortran_FLAGS "-g -O0")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback -check bounds")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -traceback -check bounds")
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -fcheck=bounds")
    endif()
else()
    message(FATAL_ERROR "Invalid value for CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
endif()


# get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
# if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#   set (CMAKE_Fortran_FLAGS "<gfortran flags go here>")
# elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
# endif()

add_link_options("-fuse-ld=lld") # Necessary to avoid "error adding symbols" error on GitHub runners
find_package(MKL CONFIG REQUIRED PATHS $ENV{MKLROOT})
find_package(FFTW REQUIRED) # Not necessary if using Intel MKL (but then you would have to remove the include statements in the src files)
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

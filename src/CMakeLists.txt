# src/CMakeLists.txt

# Collect the module source files
set(SOURCES  
  global.f90
  differential_operators.f90
  mass_flow.f90
  interpolation.f90
  equations.f90
  boundary_conditions.f90
  pressure.f90
  input_output.f90
  statistics.f90
  initialization.f90
  finalization.f90
  poisson.f90
  projection.f90
  time_integration.f90
  monitor.f90
  immersed_boundary_geometry.f90
  immersed_boundary_operators.f90
  heaviside.f90
)

# Create a library from the modules
add_library(channel_lib STATIC "${SOURCES}")

target_link_libraries(channel_lib PUBLIC LAPACK::LAPACK)
target_link_libraries(channel_lib PUBLIC MPI::MPI_Fortran)
target_link_libraries(channel_lib PUBLIC FFTW::Double)
target_link_libraries(channel_lib PUBLIC FFTW::DoubleMPI)

# target_link_libraries(channel_lib PUBLIC -lifport)
# target_link_options(channel_lib PUBLIC -lifport)

# Make sure the include directories for the modules are added
target_include_directories(channel_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})  # Add the binary directory for .mod files

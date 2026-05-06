# src/CMakeLists.txt

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
  finalization.f90
  poisson.f90

  immersed_boundary_geometry.f90
  immersed_boundary_operators.f90
  heaviside.f90
  operators_FSI.f90

  initialization.f90
  projection.f90
  time_integration.f90
  monitor.f90
)

add_library(channel_lib STATIC "${SOURCES}")

target_include_directories(channel_lib PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  ${FFTW_INCLUDE_DIRS}
)

target_link_libraries(channel_lib PUBLIC
  LAPACK::LAPACK
  MPI::MPI_Fortran
  ${FFTW_LIBRARIES}
)

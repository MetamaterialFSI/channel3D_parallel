if(CMAKE_VERSION VERSION_LESS 3.12)
    find_package(PythonInterp 3.7 REQUIRED)
    set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}")
else()
    find_package(Python3 3.7 REQUIRED COMPONENTS Interpreter)
endif()

add_executable(test_bicgstab_metrics test_bicgstab_metrics.f90)
target_link_libraries(test_bicgstab_metrics PRIVATE channel_lib)
add_test(NAME unit_bicgstab_metrics COMMAND test_bicgstab_metrics)
set_tests_properties(unit_bicgstab_metrics PROPERTIES LABELS "unit;solver")

add_executable(test_ddf_kernels test_ddf_kernels.f90)
target_link_libraries(test_ddf_kernels PRIVATE channel_lib)
add_test(NAME unit_ddf_kernels COMMAND test_ddf_kernels)
set_tests_properties(unit_ddf_kernels PROPERTIES LABELS "unit;immersed_boundary")

set(channel3d_mpi_flag_args)
foreach(flag IN LISTS MPIEXEC_PREFLAGS)
    list(APPEND channel3d_mpi_flag_args "--mpi-preflag=${flag}")
endforeach()
foreach(flag IN LISTS MPIEXEC_POSTFLAGS)
    list(APPEND channel3d_mpi_flag_args "--mpi-postflag=${flag}")
endforeach()

function(add_bicgstab_diagnostic_test)
    set(one_value_args NAME INPUT CASE INTEGRATOR)
    cmake_parse_arguments(DIAGNOSTIC "BAO5" "${one_value_args}" "" ${ARGN})

    set(input_source "${CMAKE_CURRENT_SOURCE_DIR}/${DIAGNOSTIC_INPUT}")
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/work/${DIAGNOSTIC_NAME}")
    file(MAKE_DIRECTORY "${work_dir}")
    get_filename_component(input_name "${input_source}" NAME)
    configure_file("${input_source}" "${work_dir}/${input_name}" COPYONLY)
    if(DIAGNOSTIC_BAO5)
        file(READ "${input_source}" input_text)
        string(REPLACE "&params" "&params\n  use_bao5_ddf = .true." input_text "${input_text}")
        file(WRITE "${work_dir}/${input_name}" "${input_text}")
    endif()

    add_test(
        NAME "${DIAGNOSTIC_NAME}"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_bicgstab_diagnostics.py"
                --mpiexec "${MPIEXEC_EXECUTABLE}"
                "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
                --nprocs 2
                ${channel3d_mpi_flag_args}
                --driver "$<TARGET_FILE:channel3D>"
                --input "${work_dir}/${input_name}"
                --case "${DIAGNOSTIC_CASE}"
                --integrator "${DIAGNOSTIC_INTEGRATOR}"
    )
    set_tests_properties("${DIAGNOSTIC_NAME}" PROPERTIES
        WORKING_DIRECTORY "${work_dir}"
        TIMEOUT 300
        PROCESSORS 2
        LABELS "mpi;diagnostic;solver;${DIAGNOSTIC_CASE}"
    )
endfunction()

add_bicgstab_diagnostic_test(
    NAME diagnostic_bicgstab_output
    INPUT inputs/bicgstab_diagnostic.nml
    CASE center_wall
    INTEGRATOR rk3
)

add_bicgstab_diagnostic_test(
    NAME diagnostic_bicgstab_double_output
    INPUT inputs/bicgstab_double_diagnostic.nml
    CASE double_cylinder
    INTEGRATOR rk3
)

add_bicgstab_diagnostic_test(
    NAME diagnostic_bicgstab_euler_output
    INPUT inputs/bicgstab_euler_diagnostic.nml
    CASE center_wall
    INTEGRATOR euler
)

add_bicgstab_diagnostic_test(
    NAME diagnostic_bicgstab_output_bao5
    INPUT inputs/bicgstab_diagnostic.nml
    CASE center_wall
    INTEGRATOR rk3
    BAO5
)

add_bicgstab_diagnostic_test(
    NAME diagnostic_bicgstab_double_output_bao5
    INPUT inputs/bicgstab_double_diagnostic.nml
    CASE double_cylinder
    INTEGRATOR rk3
    BAO5
)

add_bicgstab_diagnostic_test(
    NAME diagnostic_bicgstab_euler_output_bao5
    INPUT inputs/bicgstab_euler_diagnostic.nml
    CASE center_wall
    INTEGRATOR euler
    BAO5
)

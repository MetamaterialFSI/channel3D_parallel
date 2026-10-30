include(CMakeParseArguments)

function(add_channel3d_mpi_test)
    set(one_value_args NAME INPUT NPROCS TIMEOUT)
    set(multi_value_args ARGS FILES LABELS)
    cmake_parse_arguments(CHANNEL3D_TEST "" "${one_value_args}" "${multi_value_args}" ${ARGN})

    if(CHANNEL3D_TEST_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR
            "Unknown arguments for add_channel3d_mpi_test: ${CHANNEL3D_TEST_UNPARSED_ARGUMENTS}")
    endif()
    if(NOT CHANNEL3D_TEST_NAME)
        message(FATAL_ERROR "add_channel3d_mpi_test requires NAME")
    endif()
    if(NOT CHANNEL3D_TEST_INPUT)
        message(FATAL_ERROR "add_channel3d_mpi_test requires INPUT")
    endif()
    if(NOT TARGET channel3D)
        message(FATAL_ERROR "add_channel3d_mpi_test requires the channel3D target")
    endif()
    if(NOT MPIEXEC_EXECUTABLE)
        message(FATAL_ERROR "add_channel3d_mpi_test requires an MPI launcher")
    endif()

    if(NOT CHANNEL3D_TEST_NPROCS)
        set(CHANNEL3D_TEST_NPROCS 2)
    endif()
    if(NOT CHANNEL3D_TEST_TIMEOUT)
        set(CHANNEL3D_TEST_TIMEOUT 300)
    endif()

    if(IS_ABSOLUTE "${CHANNEL3D_TEST_INPUT}")
        set(input_source "${CHANNEL3D_TEST_INPUT}")
    else()
        set(input_source "${CMAKE_CURRENT_SOURCE_DIR}/${CHANNEL3D_TEST_INPUT}")
    endif()
    if(NOT EXISTS "${input_source}")
        message(FATAL_ERROR "Test input does not exist: ${input_source}")
    endif()

    string(MAKE_C_IDENTIFIER "${CHANNEL3D_TEST_NAME}" test_id)
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/work/${test_id}")
    file(MAKE_DIRECTORY "${work_dir}")

    get_filename_component(input_name "${input_source}" NAME)
    configure_file("${input_source}" "${work_dir}/${input_name}" COPYONLY)

    foreach(extra_file IN LISTS CHANNEL3D_TEST_FILES)
        if(IS_ABSOLUTE "${extra_file}")
            set(extra_source "${extra_file}")
        else()
            set(extra_source "${CMAKE_CURRENT_SOURCE_DIR}/${extra_file}")
        endif()
        if(NOT EXISTS "${extra_source}")
            message(FATAL_ERROR "Test support file does not exist: ${extra_source}")
        endif()
        get_filename_component(extra_name "${extra_source}" NAME)
        configure_file("${extra_source}" "${work_dir}/${extra_name}" COPYONLY)
    endforeach()

    add_test(
        NAME "${CHANNEL3D_TEST_NAME}"
        COMMAND "${MPIEXEC_EXECUTABLE}"
                "${MPIEXEC_NUMPROC_FLAG}" "${CHANNEL3D_TEST_NPROCS}"
                ${MPIEXEC_PREFLAGS}
                "$<TARGET_FILE:channel3D>"
                ${MPIEXEC_POSTFLAGS}
                "${input_name}"
                ${CHANNEL3D_TEST_ARGS}
    )

    set_tests_properties("${CHANNEL3D_TEST_NAME}" PROPERTIES
        WORKING_DIRECTORY "${work_dir}"
        TIMEOUT "${CHANNEL3D_TEST_TIMEOUT}"
        PROCESSORS "${CHANNEL3D_TEST_NPROCS}"
        LABELS "mpi;${CHANNEL3D_TEST_LABELS}"
    )
endfunction()

if(CMAKE_VERSION VERSION_LESS 3.12)
    find_package(PythonInterp 3.6 REQUIRED)
    set(Python3_EXECUTABLE "${PYTHON_EXECUTABLE}")
else()
    find_package(Python3 3.6 REQUIRED COMPONENTS Interpreter)
endif()

add_executable(channel3d_verification verification_driver.f90)
target_link_libraries(channel3d_verification PRIVATE channel_lib)

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

set(bicgstab_diagnostic_input_source
    "${CMAKE_CURRENT_SOURCE_DIR}/inputs/bicgstab_diagnostic.nml")
set(bicgstab_diagnostic_work_dir
    "${CMAKE_CURRENT_BINARY_DIR}/work/diagnostic_bicgstab_output")
file(MAKE_DIRECTORY "${bicgstab_diagnostic_work_dir}")
configure_file(
    "${bicgstab_diagnostic_input_source}"
    "${bicgstab_diagnostic_work_dir}/bicgstab_diagnostic.nml"
    COPYONLY
)
add_test(
    NAME diagnostic_bicgstab_output
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_bicgstab_diagnostics.py"
            --mpiexec "${MPIEXEC_EXECUTABLE}"
            "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
            --nprocs 2
            ${channel3d_mpi_flag_args}
            --driver "$<TARGET_FILE:channel3D>"
            --input "${bicgstab_diagnostic_work_dir}/bicgstab_diagnostic.nml"
            --case center_wall
)
set_tests_properties(diagnostic_bicgstab_output PROPERTIES
    WORKING_DIRECTORY "${bicgstab_diagnostic_work_dir}"
    TIMEOUT 300
    PROCESSORS 2
    LABELS "mpi;diagnostic;solver"
)

set(bicgstab_double_diagnostic_input_source
    "${CMAKE_CURRENT_SOURCE_DIR}/inputs/bicgstab_double_diagnostic.nml")
set(bicgstab_double_diagnostic_work_dir
    "${CMAKE_CURRENT_BINARY_DIR}/work/diagnostic_bicgstab_double_output")
file(MAKE_DIRECTORY "${bicgstab_double_diagnostic_work_dir}")
configure_file(
    "${bicgstab_double_diagnostic_input_source}"
    "${bicgstab_double_diagnostic_work_dir}/bicgstab_double_diagnostic.nml"
    COPYONLY
)
add_test(
    NAME diagnostic_bicgstab_double_output
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_bicgstab_diagnostics.py"
            --mpiexec "${MPIEXEC_EXECUTABLE}"
            "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
            --nprocs 2
            ${channel3d_mpi_flag_args}
            --driver "$<TARGET_FILE:channel3D>"
            --input "${bicgstab_double_diagnostic_work_dir}/bicgstab_double_diagnostic.nml"
            --case double_cylinder
)
set_tests_properties(diagnostic_bicgstab_double_output PROPERTIES
    WORKING_DIRECTORY "${bicgstab_double_diagnostic_work_dir}"
    TIMEOUT 300
    PROCESSORS 2
    LABELS "mpi;diagnostic;solver;double_cylinder"
)

function(add_channel3d_convergence_test)
    set(one_value_args NAME CASE COARSE FINE TIMEOUT)
    set(multi_value_args CHECKS LABELS)
    cmake_parse_arguments(CONVERGENCE "" "${one_value_args}" "${multi_value_args}" ${ARGN})

    if(CONVERGENCE_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR
            "Unknown arguments for add_channel3d_convergence_test: ${CONVERGENCE_UNPARSED_ARGUMENTS}")
    endif()
    foreach(required NAME CASE COARSE FINE)
        if(NOT CONVERGENCE_${required})
            message(FATAL_ERROR "add_channel3d_convergence_test requires ${required}")
        endif()
    endforeach()
    if(NOT CONVERGENCE_TIMEOUT)
        set(CONVERGENCE_TIMEOUT 300)
    endif()

    set(coarse_input "${CMAKE_CURRENT_SOURCE_DIR}/${CONVERGENCE_COARSE}")
    set(fine_input "${CMAKE_CURRENT_SOURCE_DIR}/${CONVERGENCE_FINE}")
    foreach(input_file "${coarse_input}" "${fine_input}")
        if(NOT EXISTS "${input_file}")
            message(FATAL_ERROR "Convergence-test input does not exist: ${input_file}")
        endif()
    endforeach()

    string(MAKE_C_IDENTIFIER "${CONVERGENCE_NAME}" test_id)
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/work/${test_id}")
    file(MAKE_DIRECTORY "${work_dir}")

    add_test(
        NAME "${CONVERGENCE_NAME}"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_convergence.py"
                --mpiexec "${MPIEXEC_EXECUTABLE}"
                "--numproc-flag=${MPIEXEC_NUMPROC_FLAG}"
                --nprocs 2
                ${channel3d_mpi_flag_args}
                --driver "$<TARGET_FILE:channel3d_verification>"
                --case "${CONVERGENCE_CASE}"
                --coarse "${coarse_input}"
                --fine "${fine_input}"
                ${CONVERGENCE_CHECKS}
    )
    set_tests_properties("${CONVERGENCE_NAME}" PROPERTIES
        WORKING_DIRECTORY "${work_dir}"
        TIMEOUT "${CONVERGENCE_TIMEOUT}"
        PROCESSORS 2
        LABELS "mpi;verification;${CONVERGENCE_LABELS}"
    )
endfunction()

add_channel3d_convergence_test(
    NAME verification_center_wall_x
    CASE center_wall_x
    COARSE inputs/center_wall_x_coarse.nml
    FINE inputs/center_wall_x_fine.nml
    CHECKS
        --min-linf-bulk-order 1.8
        --min-linf-support-order 1.6
        --max-fine-l1 1.5e-3
        --max-fine-l2 5.0e-3
        --max-fine-linf 3.0e-2
        --max-fine-linf-bulk 1.2e-3
        --max-fine-linf-support 1.5e-3
    LABELS center_wall
)

add_channel3d_convergence_test(
    NAME verification_center_wall_z
    CASE center_wall_z
    COARSE inputs/center_wall_z_coarse.nml
    FINE inputs/center_wall_z_fine.nml
    CHECKS
        --min-linf-bulk-order 1.8
        --min-linf-support-order 1.6
        --max-fine-l1 1.5e-3
        --max-fine-l2 5.0e-3
        --max-fine-linf 3.0e-2
        --max-fine-linf-bulk 1.2e-3
        --max-fine-linf-support 1.5e-3
    LABELS center_wall
)

add_channel3d_convergence_test(
    NAME verification_double_cylinder
    CASE double_cylinder
    COARSE inputs/double_cylinder_coarse.nml
    FINE inputs/double_cylinder_fine.nml
    CHECKS
        --min-l1-order 1.5
        --min-l2-order 1.3
        --min-linf-order 0.8
        --max-fine-l1 1.5e-3
        --max-fine-l2 4.0e-3
        --max-fine-linf 3.0e-2
    LABELS double_cylinder
    TIMEOUT 600
)

add_channel3d_convergence_test(
    NAME verification_center_wall_x_bao5
    CASE center_wall_x
    COARSE inputs/center_wall_x_fine.nml
    FINE inputs/center_wall_x_bao5_fine.nml
    CHECKS
        --bao5-ddf
        --min-linf-bulk-order 1.8
        --min-linf-support-order 1.6
        --max-fine-l1 1.5e-3
        --max-fine-l2 5.0e-3
        --max-fine-linf 3.0e-2
        --max-fine-linf-bulk 1.2e-3
        --max-fine-linf-support 1.5e-3
    LABELS center_wall bao5
)

add_channel3d_convergence_test(
    NAME verification_center_wall_z_bao5
    CASE center_wall_z
    COARSE inputs/center_wall_z_fine.nml
    FINE inputs/center_wall_z_bao5_fine.nml
    CHECKS
        --bao5-ddf
        --min-linf-bulk-order 1.8
        --min-linf-support-order 1.6
        --max-fine-l1 1.5e-3
        --max-fine-l2 5.0e-3
        --max-fine-linf 3.0e-2
        --max-fine-linf-bulk 1.2e-3
        --max-fine-linf-support 1.5e-3
    LABELS center_wall bao5
)

add_channel3d_convergence_test(
    NAME verification_double_cylinder_bao5
    CASE double_cylinder
    COARSE inputs/double_cylinder_coarse.nml
    FINE inputs/double_cylinder_fine.nml
    CHECKS
        --bao5-ddf
        --min-l1-order 1.5
        --min-l2-order 1.3
        --min-linf-order 0.8
        --max-fine-l1 1.5e-3
        --max-fine-l2 4.0e-3
        --max-fine-linf 3.0e-2
    LABELS double_cylinder bao5
    TIMEOUT 600
)

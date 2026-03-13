!------------------------------------------------------!
! My channel :)                                        !
!                                                      !
! Solve incompression Navier-Stokes eqs                !
!       with spanwise rotation                         !   
!                                                      !    
! Spatial discretization:                              !
!        - 2nd order finite differences                !
!        - Staggered mesh                              !
!                                                      !
! Temporal discretization:                             !
!        - Explicit RK3                                !
!        - Fractional step method                      !
!                                                      !
! Boundary conditions:                                 !
!        - x and z periodic                            !
!        - y non-slip/slip                             !
!                                                      !
! Parallelization:                                     !
!        - MPI, z-slices                               !
!                                                      !
! Required:                                            !
!       - FFTW 3.X                                     !
!       - LAPACK 3.X                                   !
!                                                      !
! Parallel, Version 0.23                               !
!                                                      !
! Adrian Lozano Duran                                  !
! H. Jane Bae                                          !
! 2017                                                 !
!------------------------------------------------------!
Program channel_FD

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use input_output
  Use initialization
  Use time_integration
  Use monitor
  Use statistics
  Use finalization
  Use immersed_boundary_geometry
  Use immersed_boundary_operators
  Use heaviside
  Use mpi
  
  ! prevent implicit typing
  Implicit None

  ! initialize mpi
  call Mpi_init(ierr)
  call Mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call Mpi_comm_rank(MPI_COMM_WORLD,   myid, ierr)

  ! get the name of the parameters file from the command line arguments
  Call get_command_argument(1, fileparams)

  ! read the input from the parameters file
  Call read_input_parameters

  ! initialize flow variables and read input flow field
  Call initialize

  ! initialize IB variables
  Call initialize_ib_arrays

  ! initialize IB geometry
  Call setup_IB_geometry 

  ! small summary of input parameters
  Call summary

  ! initialize IB operators
  Call setup_IB_operators

  ! compute Heaviside fields
  If ( trim(body_type) /= 'none' ) Then
    Call compute_heaviside
  Else
    Hu_interior = 1.d0
    Hv_interior = 1.d0
    Hw_interior = 1.d0
    Hu_exterior = 0.d0
    Hv_exterior = 0.d0
    Hw_exterior = 0.d0
  End If

  ! recompute initial mass flow with heaviside masking
  Call compute_mean_mass_flow_U(U,Qflow_x_0)
  Call compute_mean_mass_flow_V(V,Qflow_y_0)
  Call compute_mean_mass_flow_W(W,Qflow_z_0)
  Qflow_y_0 = 0d0
  dPdy      = 0d0
  
  ! write snapshot if needed
  Call output_data
  Call compute_statistics 
  Call output_monitor

  ! temporal loop
  Do istep = 1, nsteps
    
    ! compute dt based on CFL
    Call compute_dt

    ! time step
    Call compute_time_step_RK3

    ! compute a few statistics
    Call compute_statistics 

    ! output some key values
    Call output_monitor

    ! write snapshot if needed
    Call output_data

  End Do

  ! finalize stuff
  Call finalize

End program channel_FD

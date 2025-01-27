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
  
  ! prevent implicit typing
  Implicit None

  ! initialize everything and read input file and input flow field
  Call initialize

  ! small summary of input parameters
  Call summary

  ! write snapshot if needed
  Call output_data
     
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

Program channel_FD

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
  Use operators_FSI
  Use heaviside
  Use mpi

  Implicit None
  Integer :: i, j
  Real(Int64) :: ramp_s, ramp_fac

  Call Mpi_init(ierr)
  Call Mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  Call Mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

  Call get_command_argument(1, fileparams)

  Call read_input_parameters

  ! initialize flow variables and read input flow field
  Call initialize

  ! initialize IB variables
  Call initialize_ib_arrays

  ! initialize IB geometry
  Call setup_IB_geometry

  ! If restarting, overwrite freshly initialized IB/FSI state
  ! with xb,yb,zb and FSI variables from the restart file.
  !uncomment this only from the restart onwards, doesn't make sense for 1st run
  If ( functionality .eq. 1 .and. init_type .eq. 0 ) Then
    Call read_fsi_body_restart_data
  End If

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

  ! build mass and stiffness matrix of structure in case FSI functionality is opted
  If ( functionality .eq. 1 ) Then

    Call preprocess_FSI
  ! Mmat_testcase_actual = Mmat_testcase
   !Cmat_testcase_actual = Cmat_testcase
  ! Kmat_testcase_actual = Kmat_testcase
    If ( init_type .eq. 0 ) Then
      ! restart: xb,yb,zb are already deformed
      ! do not set reference equal to deformed state
      xbref = xb
      ybref = yb - chi_k
      zbref = zb
    Else
      ! fresh start: current body position is the reference
      xbref = xb
      ybref = yb
      zbref = zb
    End If

  End If

  ! recompute initial mass flow with heaviside masking
  !Call compute_mean_mass_flow_U(U,Qflow_x_0)
  !Call compute_mean_mass_flow_V(V,Qflow_y_0)
  !Call compute_mean_mass_flow_W(W,Qflow_z_0)
!---------------------------------------------------------!
! Preserve prescribed constant-mass-flow targets.         !
!                                                         !
! Only compute the initial flow when that direction is    !
! NOT being controlled at a prescribed mass-flow rate.    !
!---------------------------------------------------------!

If (x_mass_cte /= 1) Then
  Call compute_mean_mass_flow_U(U,Qflow_x_0)
End If

If (y_mass_cte /= 1) Then
  Call compute_mean_mass_flow_V(V,Qflow_y_0)
End If

If (z_mass_cte /= 1) Then
  Call compute_mean_mass_flow_W(W,Qflow_z_0)
End If

Qflow_y_0 = 0.0d0
dPdy      = 0.0d0
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
! uncomment to generate turbulence if needed
    !if(istep .le. 2000) then
     !   nu= 0.5*3.33333e-04
    !end if
    !if(istep .ge. 2000) then
     ! nu= 3.33333e-04
    !end if
   !If ( functionality .eq. 1 .and. istep <= 1 ) Then
    !  ramp_s = Real(istep,Int64)/1.d0
    !  ramp_fac = 1.d0**(1.d0 - ramp_s)
    !  Mmat_testcase = ramp_fac*Mmat_testcase_actual
    !  Cmat_testcase = ramp_fac*Cmat_testcase_actual
    !  Kmat_testcase = ramp_fac*Kmat_testcase_actual

   !End If
   !If ( functionality .eq. 1 .and. istep > 1 ) Then
       !Mmat_testcase = Mmat_testcase_actual
       !Cmat_testcase = Cmat_testcase_actual
      ! Kmat_testcase = Kmat_testcase_actual
  !End If


    ! time step
    !Call compute_time_step_RK3
    Call compute_time_step_RK3
    ! compute a few statistics
    Call compute_statistics_topchannel_updated
    ! output some key values
    Call output_monitor

    ! write snapshot if needed
    Call output_data
    !If ( myid == 0 ) Then
     !If ( mod(istep,100) == 0 ) Then
      !Call writestuffsurfacestress(fb)
       !End If
     !End If
    ! FSI diagnostics at each time instance
    If ( functionality .eq. 1 ) Then
      
      Call check_slip_FSI

      If ( myid == 0 ) Then
        Call writechimax(chi_k)
        Call writestressmax(fb)
        If ( mod(istep,1) == 0 ) Then
         Call output_fsi_body_restart_data
          Call writestuffiterfsi(iter_FSI)
          Call writestuffsurfacestress(fb)
          Call writestuffsurfacestress_redist(fb_redist)
          Call writestuffchi(chi_k)
          Call writestuffzeta(zeta_k)
          !Call writestuffu(U)
        End If

      End If

    End If

  End Do

  Call finalize

End Program channel_FD
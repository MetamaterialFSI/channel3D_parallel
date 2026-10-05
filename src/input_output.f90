!--------------------------------------!
!          Module for I/O              !
!--------------------------------------!
Module input_output

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use ifport
  Use pressure

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------!
  !  Read input parameters from a txt file !
  !----------------------------------------!
  Subroutine read_input_parameters

        Character(200) :: dummy_line, msg
    Real(Int64)    :: Rossby_plus, utau_
    Integer(Int32) :: ioerr, iounit
    logical :: exterior_pressure_gradient
    namelist /params/ &
      Lxp, Lzp, Ly_channel, alpha_stretch, &
      nx_global, ny_global, nz_global, nd, &
      nxb, nzb, &
      CFL, &
      nu, &
      dPdx, dPdz, Qflow_x_0, Qflow_z_0, &
      x_mass_cte, y_mass_cte, z_mass_cte, &
      nsteps, nsave, nstats, nmonitor, &
      filein, fileout, &
      nstep_init, t_init, &
      init_type, grid_type, body_type, exterior_pressure_gradient, &
      body_param_3, body_param_1, body_param_2, body_ramp_up_time, &
      min_buffer_width, cg_tol, cg_max_iter, perturb_scale, &
      functionality, nblocks, &
      m_parm, k_parm, c_parm, &
      xtopmass, ytopmass, ztopmass, sigma, &
      m_parm_testcase, k_parm_testcase, c_parm_testcase, &
      Tx, Tz, A

    ! default values
    alpha_stretch = 2.6d0
    Lxp = 9.4248
    Lzp = 3.1416
    !Lxp = 1.718
    !Lzp = 0.859
    Ly_channel = 4d0
    min_buffer_width = 0.16d0
    nd=48
    cg_tol = 1d-3
    cg_max_iter = 3000
    t_init = 0d0
    body_type = 'none'
    exterior_pressure_gradient = .True.
    x_mass_cte = 0
    y_mass_cte = 0
    z_mass_cte = 0
    dPdx = 0d0
    dPdy = 0d0
    dPdz = 0d0
    body_ramp_up_time = 0d0
    perturb_scale = 0.5d0
    Qflow_x_0 = 0.666666d0
    Qflow_z_0 = 0d0

    functionality = 0
    nblocks = 1

    m_parm = 0d0
    k_parm = 0d0
    c_parm = 0d0

    xtopmass = 0d0
    ytopmass = 0d0
    ztopmass = 0d0
    sigma = 0d0

    m_parm_testcase = 0d0
    k_parm_testcase = 0d0
    c_parm_testcase = 0d0

    Tx = 0d0
    Tz = 0d0
    A  = 1d0

    ! processor 0 reads the data
    If ( myid==0 ) Then
      Write(*,*) 'reading input parameters...'

      Open(newunit=iounit, file=fileparams, status="old", action="read", iostat=ioerr, iomsg=msg)
      If (ioerr /= 0) then
          Print *, "Error opening input file. IOSTAT =", ioerr, " IOMSG = ", msg
          Stop 1
      End If

      Read(iounit, Nml=params, IOSTAT=ioerr, IOMSG=msg)
      If (ioerr /= 0) then
          Print *, "Error reading namelist. IOSTAT =", ioerr, " IOMSG = ", msg
          Stop 1
      End If

      Close(iounit)

      utau_    = dPdx ** 0.5d0
      dPdx_ref = dPdx
      Call to_lower(body_type)

      Print params

    End If

    ! broadcast data to all processors
    Call Mpi_bcast ( nx_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ny_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nz_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nd,       1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( Lxp,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( Lzp,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( Ly_channel,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( alpha_stretch,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( nxb,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nzb,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (      CFL,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (       nu,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (     dPdx,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (     dPdz,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (Qflow_x_0,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (Qflow_z_0,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( dPdx_ref,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   t_init,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (  nstep_init,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (      nsteps,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (       nsave,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (      nstats,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    nmonitor,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   init_type,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   grid_type,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   body_type,len(body_type),MPI_character,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (exterior_pressure_gradient,1,MPI_logical,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  x_mass_cte,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  y_mass_cte,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  z_mass_cte,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( min_buffer_width,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( cg_max_iter,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( cg_tol,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( perturb_scale,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (   body_param_1,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   body_param_2,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   body_param_3,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   body_ramp_up_time,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( functionality,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nblocks,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( m_parm,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( k_parm,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( c_parm,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( xtopmass,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ytopmass,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ztopmass,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( sigma,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( m_parm_testcase,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( k_parm_testcase,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( c_parm_testcase,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( Tx,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( Tz,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( A,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

  End Subroutine read_input_parameters

  !--------------------------------------!
  !    Generates a grid for channel      !
  !                                      !
  ! Output: x_global, y_global, z_global !
  !                                      !
  !--------------------------------------!
  Subroutine create_grid

    Integer(Int32) :: i
Integer(Int32) :: n_buffer, n_lower_outer
Integer(Int32) :: n_ref, i_ref_start
Integer(Int32) :: j_plate, j_last, ny_required

Real(Int64) :: dy
Real(Int64) :: dys
Real(Int64) :: buffer_target
Real(Int64) :: buffer_width
Real(Int64) :: lower_edge

Real(Int64), Allocatable, Dimension(:) :: y_single_unpadded

    n_uniform = 8

    Do i = 1, nx_global
      x_global(i) = Real(i-1,8)
    End Do
    x_global = Lxp * x_global / x_global(nx_global - 1)

    Do i = 1, nz_global
      z_global(i) = Real(i-1,8)
    End Do
    z_global = Lzp * z_global / z_global(nz_global - 1)

    Select Case (grid_type)
      Case (0) ! Uniform grid
        If ( myid==0 ) Write(*,*) 'Generating uniform y grid'
        Do i=1,ny_global
          y_global(i) = Real(i-1,8)
        End Do
        y_global = Ly_channel * y_global / Maxval(y_global)

      Case (1) ! Stretched grid wall to wall
        If ( myid==0 ) Write(*,*) 'Generating stretched y grid'
        Do i=1,ny_global
          y_global(i) = Real(i-1,8)
        End Do
        y_global = Ly_channel * y_global / Maxval(y_global) - 1d0

        If ( alpha_stretch > 0d0 ) Then
          Do i=1,ny_global
            y_global(i) = dtanh(alpha_stretch * y_global(i)) / dtanh(alpha_stretch)
          End Do
        End If

        y_global = y_global - Minval(y_global)
        y_global = y_global * Ly_channel / Maxval(y_global)

      Case (2) ! Stretched grid centered around 0 with uniform buffers on each end to account for IB (moving with amplitude body_param_1)
        If ( myid==0 ) Write(*,*) 'Generating stretched y grid with uniform buffers'

        ! Loop until buffer condition is met
        write(*,*) 'creating stretched grid using ny_global = ', ny_global
        Do
          Call create_stretched_grid(y_global, Ly_channel, ny_global, n_uniform, alpha_stretch)

          ! Compute dy between first two points (assume uniform region at beginning)
          dymin = y_global(2) - y_global(1)

          If ( myid==0 ) Write(*,'(A,I4,A,F12.6)') ' Number of buffer cells on each side = ', n_uniform, ' -> dymin = ', dymin 
          If ( myid==0 ) Write(*,'(A,F12.6,A,F12.6,A)') ' Buffer width = ', n_uniform * dymin, ' (minimum requirement = ', &
            min_buffer_width + (2 * suppy + 2) * dymin, ')'
          ! TODO: check if 2 * suppy + 2 is the correct amount
          If (n_uniform * dymin >= min_buffer_width + (2 * suppy + 2) * dymin) Exit

          n_uniform = n_uniform + 1

          If (2 * n_uniform >= Ny) Stop 'Number of buffer points exceeds the total number of grid points' 
        End Do
        ! Move entire channel in the positive y-direction to center around y = 1
        !y_global = y_global + 1d0
        y_global = y_global + 0.5d0*Ly_channel

      Case (3) ! Asymmetric Kim-Choi grid identical to serial guess generator

        If (myid == 0) Then
          Write(*,*) 'Generating asymmetric Kim-Choi-reference y grid'
        End If

        !---------------------------------------------------------!
        ! This grid is specifically constructed for the physical !
        ! domain 0 <= y <= 4 with the IB plate located at y = 2. !
        !---------------------------------------------------------!

        If (Abs(Ly_channel - 4.0d0) > 1.0d-13) Then
          If (myid == 0) Then
            Write(*,*) 'GRID_TYPE=3 requires Ly_channel = 4'
            Write(*,*) 'Current Ly_channel = ', Ly_channel
          End If
          Error Stop 'Incorrect Ly_channel for Kim-Choi grid'
        End If

        !---------------------------------------------------------!
        ! Kim and Choi 65-point reference grid.                  !
        !                                                       !
        ! Minimum spacing: dy_min^+ = 0.47 at Re_tau = 138.     !
        !---------------------------------------------------------!

        n_ref         = 65
        alpha_stretch = 2.23906084958144d0
        dy            = 0.47d0 / 138.0d0

        ! This is the minimum wall-normal spacing.
        dymin = dy

        !---------------------------------------------------------!
        ! Uniform buffer around the plate at y = 2.              !
        !                                                        !
        ! dy = 0.47/138                                          !
        !    = 0.003405797101449275                              !
        !                                                        !
        ! Required minimum buffer half-width = 0.16.             !
        !                                                        !
        ! Ceiling(0.16/dy) = 47 intervals per side.              !
        !                                                        !
        ! Therefore:                                             !
        !   nd       = 48                                        !
        !   n_buffer = 47                                        !
        !                                                        !
        ! Actual buffer half-width:                              !
        !   47*dy = 0.1600724637681159                           !
        !                                                        !
        ! Uniform fine-grid region therefore extends from        !
        !                                                        !
        !   y = 1.839927536231884                                !
        !                                                        !
        ! to                                                     !
        !                                                        !
        !   y = 2.160072463768116                                !
        !---------------------------------------------------------!

        buffer_target = 0.16d0
        n_buffer      = nd - 1

        !---------------------------------------------------------!
        ! Inactive lower region: 20 uniform intervals from y=0   !
        ! to the lower edge of the fine IB buffer.               !
        !---------------------------------------------------------!

        n_lower_outer = 20

        !---------------------------------------------------------!
        ! Require the smallest integer number of fine intervals  !
        ! that completely covers +/-0.16.                        !
        !                                                        !
        ! For the current grid this requires:                    !
        !                                                        !
        !   n_buffer = 47                                        !
        !   nd       = 48                                        !
        !---------------------------------------------------------!

        If (n_buffer /= Ceiling(buffer_target/dy)) Then

          If (myid == 0) Then
            Write(*,*) 'Incorrect nd for GRID_TYPE = 3'
            Write(*,*) 'Current nd                  = ', nd
            Write(*,*) 'Required nd                 = ', &
                       Ceiling(buffer_target/dy) + 1
            Write(*,*) 'Requested buffer half-width = ', buffer_target
            Write(*,*) 'Available uniform dy        = ', dy
            Write(*,*) 'Required uniform intervals  = ', &
                       Ceiling(buffer_target/dy)
            Write(*,*) 'Actual buffer half-width    = ', &
                 Real(Ceiling(buffer_target/dy),Int64)*dy
          End If

          Error Stop 'GRID_TYPE=3: incorrect ND'

        End If

        buffer_width = Real(n_buffer,Int64)*dy
        lower_edge   = 2.0d0 - buffer_width

        !---------------------------------------------------------!
        ! Construct Kim-Choi 65-point reference coordinates over !
        ! a channel extending from 0 to 2.                       !
        !                                                        !
        ! The retained coordinates will later be shifted upward  !
        ! by 2 into the physical upper channel.                  !
        !---------------------------------------------------------!

        Allocate(y_single_unpadded(n_ref))

        Do i = 1, n_ref

          dys = -1.0d0 + 2.0d0*Real(i-1,Int64) / &
                             Real(n_ref-1,Int64)

          y_single_unpadded(i) = 1.0d0 + &
               dtanh(alpha_stretch*dys) / dtanh(alpha_stretch)

        End Do

        !---------------------------------------------------------!
        ! Find first Kim-Choi reference coordinate strictly      !
        ! above the upper edge of the uniform IB buffer.         !
        !                                                        !
        ! For the current +/-0.16 buffer, this is reference      !
        ! index 17.                                              !
        !---------------------------------------------------------!

        i_ref_start = n_ref

        Do i = 2, n_ref

          If (y_single_unpadded(i) > buffer_width) Then
            i_ref_start = i
            Exit
          End If

        End Do

        !---------------------------------------------------------!
        ! Required number of y coordinates:                      !
        !                                                        !
        ! 1 bottom boundary                                      !
        ! + 20 lower-region intervals                            !
        ! + 47 lower-buffer intervals                            !
        ! + 47 upper-buffer intervals                            !
        ! + retained Kim-Choi points.                            !
        !                                                        !
        ! For the current configuration:                         !
        !                                                        !
        !   i_ref_start = 17                                     !
        !   NY_GLOBAL   = 164                                    !
        !---------------------------------------------------------!

        ny_required = 1 + n_lower_outer + 2*n_buffer + &
                      (n_ref - i_ref_start + 1)

        If (ny_global /= ny_required) Then

          If (myid == 0) Then
            Write(*,*) 'Incorrect NY_GLOBAL for GRID_TYPE = 3'
            Write(*,*) 'Current NY_GLOBAL           = ', ny_global
            Write(*,*) 'Required NY_GLOBAL          = ', ny_required
            Write(*,*) 'Current nd                  = ', nd
            Write(*,*) 'Reference start index       = ', i_ref_start
            Write(*,*) 'Lower outside intervals     = ', n_lower_outer
            Write(*,*) 'Uniform intervals per side  = ', n_buffer
          End If

          Error Stop 'GRID_TYPE=3: incorrect NY_GLOBAL'

        End If

        !---------------------------------------------------------!
        ! Lower inactive region.                                 !
        !                                                        !
        ! 20 uniformly spaced intervals from y=0 to the lower    !
        ! edge of the fine IB buffer.                            !
        !---------------------------------------------------------!

        y_global(1) = 0.0d0

        Do i = 1, n_lower_outer

          y_global(i+1) = lower_edge * Real(i,Int64) / &
                                      Real(n_lower_outer,Int64)

        End Do

        ! Force exact lower buffer boundary.

        y_global(n_lower_outer+1) = lower_edge

        !---------------------------------------------------------!
        ! Uniform lower IB buffer:                               !
        !                                                        !
        ! 2-buffer_width <= y <= 2                               !
        !---------------------------------------------------------!

        Do i = 1, n_buffer

          y_global(n_lower_outer+1+i) = lower_edge + &
                                        dy*Real(i,Int64)

        End Do

        j_plate = n_lower_outer + 1 + n_buffer

        ! Force plate coordinate to exactly y=2.

        y_global(j_plate) = 2.0d0

        !---------------------------------------------------------!
        ! Uniform upper IB buffer:                               !
        !                                                        !
        ! 2 <= y <= 2+buffer_width                               !
        !---------------------------------------------------------!

        Do i = 1, n_buffer

          y_global(j_plate+i) = 2.0d0 + &
                                dy*Real(i,Int64)

        End Do

        j_last = j_plate + n_buffer

        !---------------------------------------------------------!
        ! Add remaining Kim-Choi reference coordinates.          !
        !                                                        !
        ! Original reference positions span 0 <= y <= 2.        !
        ! Adding 2 shifts them into the physical upper channel.  !
        !---------------------------------------------------------!

        Do i = i_ref_start, n_ref

          j_last = j_last + 1

          y_global(j_last) = 2.0d0 + &
                             y_single_unpadded(i)

        End Do

        !---------------------------------------------------------!
        ! Safety checks.                                         !
        !---------------------------------------------------------!

        If (j_last /= ny_global) Then

          If (myid == 0) Then
            Write(*,*) 'Last assigned y index = ', j_last
            Write(*,*) 'Expected final index  = ', ny_global
          End If

          Error Stop 'GRID_TYPE=3: internal point-count error'

        End If

        If (Any(y_global(2:ny_global) <= &
                y_global(1:ny_global-1))) Then

          If (myid == 0) Then

            Do i = 1, ny_global-1

              If (y_global(i+1) <= y_global(i)) Then

                Write(*,*) 'Non-increasing grid at index ', i
                Write(*,*) 'y_global(i)   = ', y_global(i)
                Write(*,*) 'y_global(i+1) = ', y_global(i+1)

              End If

            End Do

          End If

          Error Stop 'GRID_TYPE=3: y grid is not increasing'

        End If

        If (Abs(y_global(1)) > 1.0d-13) Then

          If (myid == 0) Then
            Write(*,*) 'y_global(1) = ', y_global(1)
          End If

          Error Stop 'GRID_TYPE=3: lower wall is not at y=0'

        End If

        If (Abs(y_global(j_plate)-2.0d0) > 1.0d-13) Then

          If (myid == 0) Then
            Write(*,*) 'Plate coordinate = ', y_global(j_plate)
          End If

          Error Stop 'GRID_TYPE=3: plate is not at y=2'

        End If

        If (Abs(y_global(ny_global)-4.0d0) > 1.0d-12) Then

          If (myid == 0) Then
            Write(*,*) 'Upper-wall coordinate = ', &
                       y_global(ny_global)
          End If

          Error Stop 'GRID_TYPE=3: upper wall is not at y=4'

        End If

        !---------------------------------------------------------!
        ! Print grid summary from processor zero only.            !
        !---------------------------------------------------------!

        If (myid == 0) Then

          Write(*,*) 'Asymmetric Kim-Choi grid successfully generated'
          Write(*,*) '  NY_GLOBAL                    = ', ny_global
          Write(*,*) '  ND                           = ', nd
          Write(*,*) '  plate face index             = ', j_plate
          Write(*,*) '  top-channel face points      = ', &
                     ny_global-j_plate+1
          Write(*,*) '  lower outside intervals      = ', &
                     n_lower_outer
          Write(*,*) '  points strictly below buffer = ', &
                     n_lower_outer
          Write(*,*) '  lower outside uniform dy     = ', &
                     lower_edge/Real(n_lower_outer,Int64)
          Write(*,*) '  uniform intervals per side   = ', &
                     n_buffer
          Write(*,*) '  uniform buffer dy            = ', dy
          Write(*,*) '  uniform dy in plus units     = ', &
                     dy*138.0d0
          Write(*,*) '  target buffer half-width     = ', &
                     buffer_target
          Write(*,*) '  actual buffer half-width     = ', &
                     buffer_width
          Write(*,*) '  buffer lower edge            = ', &
                     lower_edge
          Write(*,*) '  buffer upper edge            = ', &
                     2.0d0+buffer_width
          Write(*,*) '  first retained KC index      = ', &
                     i_ref_start
          Write(*,*) '  first retained KC location   = ', &
                     y_single_unpadded(i_ref_start)
          Write(*,*) '  first retained physical y    = ', &
                     2.0d0+y_single_unpadded(i_ref_start)
          Write(*,*) '  y minimum                    = ', &
                     y_global(1)
          Write(*,*) '  y at plate                   = ', &
                     y_global(j_plate)
          Write(*,*) '  y maximum                    = ', &
                     y_global(ny_global)

        End If

        Deallocate(y_single_unpadded)















    End Select

  End Subroutine create_grid

  Subroutine create_stretched_grid(grid, L, n_total, n_uniform, alpha)
    Implicit None

    Integer(Int32), Intent(In) :: n_total, n_uniform
    Real(Int64), Intent(In) :: L, alpha
    Real(Int64), Dimension(ny), Intent(Out) :: grid

    Real(Int64), Dimension(:), Allocatable :: grid_unpadded
    Real(Int64) :: ds_min_unscaled, ds_min_scaled
    Integer(Int32) :: i

    ! We use n_uniform - 1 because the first cell of the stretched grid will be considered to be part of the uniform grid
    Allocate ( grid_unpadded ( n_total - 2 * (n_uniform - 1) ) )

    Do i = 1, (n_total - 2 * (n_uniform - 1))
      grid_unpadded(i) = Real(i - 1, 8)
    End Do

    ! Make the unpadded grid go from -1 to 1
    grid_unpadded = 2d0 * grid_unpadded / Maxval(grid_unpadded) - 1d0

    ! Stretch the unpadded grid
    If ( alpha_stretch > 0d0 ) Then
      grid_unpadded = dtanh(alpha * grid_unpadded) / dtanh(alpha)
    End If

    ! Compute the smallest spacing of the unpadded grid
    ds_min_unscaled = grid_unpadded(2) - grid_unpadded(1)
    
    ! Compute the smallest spacing of the final scaled grid such that the stretched portion plus two half-buffers cover L
    ds_min_scaled = L * ds_min_unscaled / (2 + (n_uniform - 1) * ds_min_unscaled)

    ! Create the scaled grid
    grid(n_uniform : n_total - n_uniform + 1) = grid_unpadded * ds_min_scaled / ds_min_unscaled
    Do i = 2, n_uniform
      grid(n_total - n_uniform + i) = grid(n_total - n_uniform + i - 1) + ds_min_scaled
      grid(n_uniform - i + 1) = grid(n_uniform - i + 2) - ds_min_scaled
    End Do

  End Subroutine create_stretched_grid

  !------------------------------------------------!
  !    Generates an initial condition for channel  !
  !                                                !
  ! Output: U,V,W                                  !
  !                                                !
  !------------------------------------------------!
  Subroutine init_flow

  Integer(Int32) :: ii, jj, kk
  Real(Int64)    :: ym_val
  Real(Int64)    :: r

  Select Case (init_type)

    !============================================================!
    ! Read an existing restart                                  !
    !============================================================!
    Case (0)

      If (myid == 0) Write(*,*) 'Reading input data'
      Call read_input_data


    !============================================================!
    ! Zero initial condition                                    !
    !============================================================!
    Case (1)

      If (myid == 0) Write(*,*) 'Generating zero initial condition'

      Call create_grid

      U = 0.0d0
      V = 0.0d0
      W = 0.0d0


    !============================================================!
    ! Random turbulent initial condition for the TOP CHANNEL     !
    !                                                            !
    ! Physical fluid channel:                                    !
    !                                                            !
    !                2 < y < 4                                   !
    !                                                            !
    ! Region y <= 2 is kept completely at rest.                  !
    !                                                            !
    ! This matches the serial rigid-guess initialization.        !
    !============================================================!
    Case (2)

      If (myid == 0) Then
        Write(*,*) 'Generating random top-channel initial condition'
      End If

      Call create_grid

      !----------------------------------------------------------!
      ! Start with absolutely zero velocity everywhere.          !
      ! This is important because the lower region y <= 2 must   !
      ! remain inactive.                                         !
      !----------------------------------------------------------!

      U = 0.0d0
      V = 0.0d0
      W = 0.0d0


      !==========================================================!
      ! U VELOCITY                                               !
      !                                                          !
      ! U is staggered in y, so its physical y-location is the   !
      ! midpoint between consecutive y_global face coordinates.  !
      !                                                          !
      ! Serial initialization:                                   !
      !                                                          !
      ! U = (y-2)*(4-y) + perturbation                           !
      !                                                          !
      ! only for y > 2.                                          !
      !==========================================================!

      Do jj = 1, ny_global-1

        ym_val = 0.5d0 * &
                 (y_global(jj) + y_global(jj+1))

        If (ym_val > 2.0d0 .and. ym_val < Ly_channel) Then

          ! Match the serial code: one random perturbation
          ! for this wall-normal plane.
          Call random_number(r)

          U(:,jj+1,:) = &
               (ym_val - 2.0d0) * (Ly_channel - ym_val) &
               + perturb_scale * (r - 0.5d0)

        Else

          U(:,jj+1,:) = 0.0d0

        End If

      End Do

      ! Ghost values at physical outer walls.
      U(:,1,:)             = -U(:,2,:)
      U(:,ny_global+1,:)   = -U(:,ny_global,:)


      !==========================================================!
      ! V VELOCITY                                               !
      !                                                          !
      ! V lives directly on y_global faces.                      !
      !                                                          !
      ! Add random perturbations only for physical y > 2.        !
      ! y = 2 itself remains zero.                               !
      !==========================================================!

      Do ii = 1, nxg_global

        Do jj = 1, ny_global

          If (y_global(jj) > 2.0d0 .and. &
              y_global(jj) < Ly_channel) Then

            Do kk = 1, nzg

              Call random_number(r)

              V(ii,jj,kk) = &
                   perturb_scale * (r - 0.5d0)

            End Do

          Else

            V(ii,jj,:) = 0.0d0

          End If

        End Do

      End Do

      ! Exact no-penetration at the two physical outer walls.
      V(:,1,:)         = 0.0d0
      V(:,ny_global,:) = 0.0d0


      !==========================================================!
      ! W VELOCITY                                               !
      !                                                          !
      ! W has the same wall-normal staggering as U.              !
      ! Therefore use the cell-midpoint coordinate ym_val.       !
      !==========================================================!

      Do jj = 1, ny_global-1

        ym_val = 0.5d0 * &
                 (y_global(jj) + y_global(jj+1))

        If (ym_val > 2.0d0 .and. ym_val < Ly_channel) Then

          Do ii = 1, nxg_global

            Do kk = 1, nz

              Call random_number(r)

              W(ii,jj+1,kk) = &
                   perturb_scale * (r - 0.5d0)

            End Do

          End Do

        Else

          W(:,jj+1,:) = 0.0d0

        End If

      End Do

      ! Ghost values at physical outer walls.
      W(:,1,:)             = -W(:,2,:)
      W(:,ny_global+1,:)   = -W(:,ny_global,:)


    !============================================================!
    ! Exact laminar initial condition for TOP CHANNEL only       !
    !                                                            !
    ! This is also corrected so that it no longer creates a      !
    ! parabolic flow through the inactive y < 2 region.          !
    !============================================================!
    Case (3)

      If (myid == 0) Then
        Write(*,*) 'Generating exact top-channel laminar condition'
      End If

      Call create_grid

      U = 0.0d0
      V = 0.0d0
      W = 0.0d0

      Do jj = 1, ny_global-1

        ym_val = 0.5d0 * &
                 (y_global(jj) + y_global(jj+1))

        If (ym_val > 2.0d0 .and. ym_val < Ly_channel) Then

          U(:,jj+1,:) = &
               dPdx / (2.0d0*nu) * &
               (ym_val - 2.0d0) * &
               (Ly_channel - ym_val)

        Else

          U(:,jj+1,:) = 0.0d0

        End If

      End Do

      U(:,1,:)           = -U(:,2,:)
      U(:,ny_global+1,:) = -U(:,ny_global,:)

      V = 0.0d0
      W = 0.0d0


    Case Default

      If (myid == 0) Then
        Write(*,*) 'Unknown init_type = ', init_type
      End If

      Error Stop 'Invalid init_type in init_flow'

  End Select


  !==============================================================!
  ! Initial-condition diagnostics                                !
  !==============================================================!

  If (myid == 0) Then

    Write(*,*) 'Max U  ', MaxVal(U)
    Write(*,*) 'Max V  ', MaxVal(V)
    Write(*,*) 'Max W  ', MaxVal(W)

    Write(*,*) 'Mean U ', &
         Sum(U) / Real(nx_global*nyg_global*nzg_global,8)

    Write(*,*) 'Mean V ', &
         Sum(V) / Real(nxg_global*ny_global*nzg_global,8)

    Write(*,*) 'Mean W ', &
         Sum(W) / Real(nxg_global*nyg_global*nz_global,8)

  End If

End Subroutine init_flow

  !--------------------------------------------!
  !    Read binary snapshot: mesh, U,V and W   !
  !                                            !
  ! Input:  filein                             !
  ! Output: U,V,W,x,y,z                        !
  !                                            !
  !--------------------------------------------!
  Subroutine read_input_data

    Integer(Int32) ::  nx_global_f,  ny_global_f,  nz_global_f, iproc, nze, nzge
    Integer(Int32) :: nxm_global_f, nym_global_f, nzm_global_f, nn(3), ndum
    Integer(Int64) :: pos_header, nsize_U, nsize_V, ii, jj, kk

    ! processor 0 Reads the all the data
    If ( myid==0 ) Then

      Write(*,*) 'reading ',Trim(Adjustl(filein)),'...'
      Open(1,file=filein,access='stream',form='unformatted',action='Read',convert='big_endian')

      ! mesh
      Read(1) nx_global_f
      If ( nx_global_f/=nx_global ) Stop 'nx_f/=nx'
      Read(1) x_global

      Read(1) ny_global_f
      If ( ny_global_f/=ny_global ) Stop 'ny_f/=ny'
      Read(1) y_global

      Read(1) nz_global_f
      If ( nz_global_f/=nz_global ) Stop 'nz_f/=nz'
      Read(1) z_global

      Read(1) nxm_global_f
      If ( nxm_global_f/=nxm_global ) Stop 'nxm_f/=nxm'
      Read(1) xm_global

      Read(1) nym_global_f
      If ( nym_global_f/=nym_global ) Stop 'nym_f/=nym'
      Read(1) ym_global

      Read(1) nzm_global_f
      If ( nzm_global_f/=nzm_global ) Stop 'nzm_f/=nzm'
      Read(1) zm_global

      ! get header position and size
      Inquire(1,pos=pos_header)
      pos_header = pos_header - 1
      nsize_U    = nx_global*nyg_global*nzg_global*8
      nsize_V    = nxg_global*ny_global*nzg_global*8

    End If

    ! U
    If ( myid==0 ) Then
      ! read dummy
      Read(1) nn 
      If ( nn(1)/=nx_global .or. nn(2)/=nyg_global .or. nn(3)/=nzg_global ) Then 
        Write(*,*) 'nn',nn
        Stop 'Error! wrong size in input file (U)'
      End If
      ! read data for processor 0
      nzge = kg2_global(myid) - kg1_global(myid) + 1
      Read(1) U(:,:,1:nzge)
      ! data for processor n>0    
      Do iproc = 1, nprocs-1
        nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
        If ( iproc<nprocs-1 ) Then
          ndum = fseek(1,-2*nx_global*nyg_global*8,seek_cur) ! ghost cell
          Read(1) Uo(:,:,1:nzge)
          Call Mpi_send(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        Else ! especial case: U has different size for last processor
          ndum = fseek(1,-2*nx_global*nyg_global*8,seek_cur) ! ghost cell
          Read(1) Uoo(:,:,1:nzge)
          Call Mpi_send(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        End If
      Enddo       
    Else
      Call Mpi_recv(U,nx*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! V
    If ( myid==0 ) Then
      ! go to correct position. I dont know, if I dont do this it gets lost sometimes
      ndum = fseek(1,pos_header+3*4+nsize_U,seek_set)
      ! read dummy
      Read(1) nn
      If ( nn(1)/=nxg_global .or. nn(2)/=ny_global .or. nn(3)/=nzg_global ) Then 
         Write(*,*) 'nn',nn
         Stop 'Error! wrong size in input file (V)'
      End If
      ! read data for processor 0
      nzge = kg2_global(myid) - kg1_global(myid) + 1
      Read(1) V(:,:,1:nzge)
      ! data for processor n>0    
      Do iproc = 1, nprocs-1
        nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
        If ( iproc<nprocs-1 ) Then
          ndum = fseek(1,-2*nxg_global*ny_global*8,seek_cur) ! ghost cell
          Read(1) Vo(:,:,1:nzge) 
          Call Mpi_send(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        Else ! especial case: V has different size for last processor
          ndum = fseek(1,-2*nxg_global*ny_global*8,seek_cur) ! ghost cell
          Read(1) Voo(:,:,1:nzge) 
          Call Mpi_send(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        End If
      Enddo       
    Else
      Call Mpi_recv(V,nxg*ny*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! W
    If ( myid==0 ) Then
      ! go to correct position. I dont know, if I dont do this it gets lost sometimes
      ndum = fseek(1,pos_header+3*4+nsize_U+3*4+nsize_V,seek_set)
      ! read dummy
      Read(1) nn
      If ( nn(1)/=nxg_global .or. nn(2)/=nyg_global .or. nn(3)/=nz_global ) Then 
        Write(*,*) 'nn',nn
        Stop 'Error! wrong size in input file (W)'
      End If
      ! read data for processor 0
      nzge = k2_global(myid) - k1_global(myid) + 1
      Read(1) W(:,:,1:nzge)
      ! data for processor n>0    
      Do iproc = 1, nprocs-1
        nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
        If ( iproc<nprocs-1 ) Then
          ndum = fseek(1,-2*nxg_global*nyg_global*8,seek_cur) ! ghost cell
          Read(1) Wo(:,:,1:nzge)
          Call Mpi_send(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        Else ! especial case: W has different size for last processor
          ndum = fseek(1,-2*nxg_global*nyg_global*8,seek_cur) ! ghost cell
          Read(1) Woo(:,:,1:nzge)
          Call Mpi_send(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        End If
      Enddo       
    Else
      Call Mpi_recv(W,nxg*nyg*nz,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! close file
    If (myid==0) Then
      Close(1)
    End If

    ! send data to all other processors
    ! mesh
    Call Mpi_bcast ( x_global,nx_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( y_global,ny_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( z_global,nz_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( xm_global,nxm_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ym_global,nym_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( zm_global,nzm_global,MPI_real8,0,MPI_COMM_WORLD,ierr ) 


    ! set solution for zero step
    Uo = U
    Vo = V
    Wo = W

  End Subroutine read_input_data

Subroutine read_fsi_body_restart_data

    Character(200) :: fname
    Integer(Int32) :: ndum
    Integer(Int32) :: nxb_f, nzb_f

    If ( functionality /= 1 ) Return
    If ( init_type /= 0 ) Return

    fname = Trim(Adjustl(filein))//'.fsi'

    If ( myid == 0 ) Then

      Write(*,*) 'reading FSI/body restart from ', Trim(Adjustl(fname)), '...'

      Open(99,file=fname,access='stream',form='unformatted', &
           action='read',convert='big_endian')

      ! body position
      Read(99) ndum
      If ( ndum /= Size(xb) ) Stop 'restart size mismatch: xb'
      Read(99) xb

      Read(99) ndum
      If ( ndum /= Size(yb) ) Stop 'restart size mismatch: yb'
      Read(99) yb

      Read(99) ndum
      If ( ndum /= Size(zb) ) Stop 'restart size mismatch: zb'
      Read(99) zb

      Read(99) nxb_f
      If ( nxb_f /= nxb ) Stop 'restart size mismatch: nxb'

      Read(99) nzb_f
      If ( nzb_f /= nzb ) Stop 'restart size mismatch: nzb'

      ! FSI variables
      Read(99) ndum
      If ( ndum /= Size(chi) ) Stop 'restart size mismatch: chi'
      Read(99) chi

      Read(99) ndum
      If ( ndum /= Size(zeta) ) Stop 'restart size mismatch: zeta'
      Read(99) zeta

      Read(99) ndum
      If ( ndum /= Size(zetadot) ) Stop 'restart size mismatch: zetadot'
      Read(99) zetadot

      Read(99) ndum
      If ( ndum /= Size(chi_k) ) Stop 'restart size mismatch: chi_k'
      Read(99) chi_k

      Read(99) ndum
      If ( ndum /= Size(zeta_k) ) Stop 'restart size mismatch: zeta_k'
      Read(99) zeta_k

      Read(99) ndum
      If ( ndum /= Size(zetadot_k) ) Stop 'restart size mismatch: zetadot_k'
      Read(99) zetadot_k

      Close(99)

    End If

    Call Mpi_bcast(xb,Size(xb),MPI_real8,0,MPI_COMM_WORLD,ierr)
    Call Mpi_bcast(yb,Size(yb),MPI_real8,0,MPI_COMM_WORLD,ierr)
    Call Mpi_bcast(zb,Size(zb),MPI_real8,0,MPI_COMM_WORLD,ierr)

    Call Mpi_bcast(chi,Size(chi),MPI_real8,0,MPI_COMM_WORLD,ierr)
    Call Mpi_bcast(zeta,Size(zeta),MPI_real8,0,MPI_COMM_WORLD,ierr)
    Call Mpi_bcast(zetadot,Size(zetadot),MPI_real8,0,MPI_COMM_WORLD,ierr)

    Call Mpi_bcast(chi_k,Size(chi_k),MPI_real8,0,MPI_COMM_WORLD,ierr)
    Call Mpi_bcast(zeta_k,Size(zeta_k),MPI_real8,0,MPI_COMM_WORLD,ierr)
    Call Mpi_bcast(zetadot_k,Size(zetadot_k),MPI_real8,0,MPI_COMM_WORLD,ierr)

End Subroutine read_fsi_body_restart_data
  !--------------------------------------------!
  !    write binary snapshot: mesh, U,V and W  !
  !                                            !
  ! Input: U,V,W,x,y,z,xm,ym,zm                !
  ! Output: fileout                            !
  !                                            !
  !--------------------------------------------!
  Subroutine output_data

    Character(200)   :: fname
    Character(8)     :: ext
    Integer  (Int32) :: iproc, nze, nzge
    
    If ( Mod(istep,nsave)==0 ) then

      ! P
      !If (pressure_computed==.False.) Then
      !   Call compute_pressure
      !End If

      ! processor 0 writes the data
      If ( myid==0 ) Then
        
        Write(ext,'(I8)') istep + nstep_init
        
        fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))
        Write(*,*) 'writing ',Trim(Adjustl(fname))
        Open(1,file=fname,access='stream',form='unformatted',action='write',convert='big_endian')
        
        ! mesh
        Write(1) Shape(x_global), x_global
        Write(1) Shape(y_global), y_global
        Write(1) Shape(z_global), z_global
        
        Write(1) Shape(xm_global), xm_global
        Write(1) Shape(ym_global), ym_global
        Write(1) Shape(zm_global), zm_global          
       
      End If

      ! U
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(U,nx*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write U size
        Write(1) nx_global,nyg_global,nzg_global
        ! processor 0 writes its data
        Write(1) U(:,:,1:nzg-1) 
        ! processor 0 receives and writes rest data
        Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Uo(:,:,2:nzge-1)
          Else
            Call Mpi_recv(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Uoo(:,:,2:nzge)
          End If
        End Do
      Endif

      ! V
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(V,nxg*ny*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write V size
        Write(1) nxg_global,ny_global,nzg_global
        ! processor 0 writes its data
        Write(1) V(:,:,1:nzg-1)
        ! processor 0 receives and write rest data
        Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Vo(:,:,2:nzge-1)
          Else
            Call Mpi_recv(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Voo(:,:,2:nzge)
          End If
        End Do
      Endif

      ! W
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(W,nxg*nyg*nz,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write W size
        Write(1) nxg_global,nyg_global,nz_global
        ! processor 0 writes its data
        Write(1) W(:,:,1:nz-1)
        ! processor 0 receives and writes rest data
        Do iproc = 1, nprocs-1
          nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Wo(:,:,2:nze-1)
          Else
            Call Mpi_recv(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Woo(:,:,2:nze)
          End If
        End Do
      Endif

      ! P
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(P,nxg*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write P size
        Write(1) nxg_global,nyg_global,nzg_global
        ! processor 0 writes its data
        Write(1) P(:,:,1:nzg-1)
        ! processor 0 receives and write rest data
        Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Po,nxg*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Po(:,:,2:nzge-1)
          Else
            Call Mpi_recv(Poo,nxg*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Poo(:,:,2:nzge)
          End If
        End Do
      Endif

      ! Hu_interior
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(Hu_interior,nx*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write Hu_interior size
        Write(1) nx_global,nyg_global,nzg_global
        ! processor 0 writes its data
        Write(1) Hu_interior(:,:,1:nzg-1) 
        ! processor 0 receives and writes rest data
        Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Hu_interior_o,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Hu_interior_o(:,:,2:nzge-1)
          Else
            Call Mpi_recv(Hu_interior_oo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Hu_interior_oo(:,:,2:nzge)
          End If
        End Do
      Endif

      ! Hv_interior
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(Hv_interior,nxg*ny*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write Hv_interior size
        Write(1) nxg_global,ny_global,nzg_global
        ! processor 0 writes its data
        Write(1) Hv_interior(:,:,1:nzg-1)
        ! processor 0 receives and write rest data
        Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Hv_interior_o,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Hv_interior_o(:,:,2:nzge-1)
          Else
            Call Mpi_recv(Hv_interior_oo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Hv_interior_oo(:,:,2:nzge)
          End If
        End Do
      Endif

      ! Hw_interior
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(Hw_interior,nxg*nyg*nz,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write W size
        Write(1) nxg_global,nyg_global,nz_global
        ! processor 0 writes its data
        Write(1) Hw_interior(:,:,1:nz-1)
        ! processor 0 receives and writes rest data
        Do iproc = 1, nprocs-1
          nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Hw_interior_o,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Hw_interior_o(:,:,2:nze-1)
          Else
            Call Mpi_recv(Hw_interior_oo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Hw_interior_oo(:,:,2:nze)
          End If
        End Do
      Endif

      If ( myid==0 ) Then
        ! Time
        Write(1) t

        ! Timestep
        Write(1) dt

        ! Pressure gradient
        Write(1) dpdx

        ! Viscosity
        Write(1) nu

        ! body
        Write(1) Shape(xb, Int32), xb
        Write(1) Shape(yb, Int32), yb
        Write(1) Shape(zb, Int32), zb
        Write(1) nxb
        Write(1) nzb

        ! FSI restart variables
        If ( functionality .eq. 1 ) Then
          Write(1) Shape(chi, Int32), chi
          Write(1) Shape(zeta, Int32), zeta
          Write(1) Shape(zetadot, Int32), zetadot

          Write(1) Shape(chi_k, Int32), chi_k
          Write(1) Shape(zeta_k, Int32), zeta_k
          Write(1) Shape(zetadot_k, Int32), zetadot_k
        End If

        ! surface stress
        Write(1) Shape(fb, Int32), fb

        ! body velocity
        Write(1) Shape(ub, Int32), ub

        ! body surface areas
        Write(1) Shape(sb, Int32), sb

        ! body normals and tangents
        Write(1) Shape(normals, Int32), normals
        Write(1) Shape(tangents_1, Int32), tangents_1
        Write(1) Shape(tangents_2, Int32), tangents_2

        ! array output for debugging
        Write(1) Shape(debug_surface_scalar, Int32), debug_surface_scalar

      End If
         
      ! close file
      If (myid==0) Then
        Close(1)
      End If
      
    End If
       
  End Subroutine output_data

Subroutine output_fsi_body_restart_data

    Character(200) :: fname
    Character(8)   :: ext

    If ( functionality /= 1 ) Return

    If ( myid == 0 ) Then

      Write(ext,'(I8)') istep + nstep_init
      fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))//'.fsi'

      Open(99,file=fname,access='stream',form='unformatted', &
           action='write',convert='big_endian')

      Write(99) Size(xb), xb
      Write(99) Size(yb), yb
      Write(99) Size(zb), zb
      Write(99) nxb
      Write(99) nzb

      Write(99) Size(chi), chi
      Write(99) Size(zeta), zeta
      Write(99) Size(zetadot), zetadot

      Write(99) Size(chi_k), chi_k
      Write(99) Size(zeta_k), zeta_k
      Write(99) Size(zetadot_k), zetadot_k

      Close(99)

    End If

End Subroutine output_fsi_body_restart_data
  !----------------------------------------------!
  !   Write some basic statistics in a txt file  !
  !----------------------------------------------!
  Subroutine output_statistics

    Character(200) :: fname
    Character(8)   :: ext
    Integer(Int32) :: jj

    If ( myid==0 ) Then

       Write(ext,'(I8)') istep + nstep_init
       
       fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))//'.stats.txt'
       Write(*,*) 'writing ',Trim(Adjustl(fname))
       Open(3,file=fname,form='formatted',action='write') 
       Write(3,'(A,4F15.8,4I)') '%',t, Retau_u, utau, nu, nx_global, ny_global, nz_global, istep
       Do jj=1,nyg
          Write(3,'(8F15.8)') yg(jj), Umean(jj), Vmean(jj), Wmean(jj), U2mean(jj), V2mean(jj), W2mean(jj), UVmean(jj)
       End Do
       Close(3)

    End If

  End Subroutine output_statistics

  !----------------------------------------------------------------!
  !   Subroutine to change the case of a string to all lower case  !
  !----------------------------------------------------------------!
  Subroutine to_lower(str)
    character(*), intent(in out) :: str
    integer :: i

    Do i = 1, len(str)
      Select Case(str(i:i))
        case("A":"Z")
          str(i:i) = achar(iachar(str(i:i))+32)
      End Select
    End Do  
  End Subroutine to_Lower

End Module input_output

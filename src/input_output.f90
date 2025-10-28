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

    namelist /params/ &
      Lxp, Lzp, alpha_stretch, &
      nx_global, ny_global, nz_global, &
      nxb, nzb, &
      CFL, &
      nu, &
      dPdx, dPdz, x_mass_cte, y_mass_cte, &
      nsteps, nsave, nstats, nmonitor, &
      filein, fileout, &
      nstep_init, t_init, &
      init_type, grid_type, body_type, &
      body_param_3, body_param_1, body_param_2, &
      min_buffer_width, cg_tol, cg_max_iter

    ! default values
    alpha_stretch = 2.6d0
    Lxp = 2d0 * pi
    Lzp = 1d0 * pi
    min_buffer_width = 0d0
    cg_tol = 1e-8
    cg_max_iter = 50
    t_init = 0d0
    body_type = 'none'

    ! processor 0 reads the data
    If ( myid==0 ) Then
      Write(*,*) 'reading input parameters...'

      Open(newunit=iounit, file=fileparams, status="old", action="read", iostat=ioerr, iomsg=msg)
      If (ierr /= 0) then
          Print *, "Error opening input file. IOSTAT =", ioerr, " IOMSG = ", msg
          Stop 1
      End If
      Read(iounit, Nml=params, IOSTAT=ioerr, IOMSG=msg)
      If (ierr /= 0) then
          Print *, "Error reading namelist. IOSTAT =", ioerr, " IOMSG = ", msg
          Stop 1
      End If

      utau_    = dPdx ** 0.5d0
      dPdx_ref = dPdx
      Call to_lower(body_type)

      Print params
       
    End If

    ! broadcast data to all processors
    Call Mpi_bcast ( nx_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ny_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nz_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( Lxp,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( Lzp,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( alpha_stretch,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( nxb,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nzb,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (      CFL,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (       nu,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (     dPdx,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (     dPdz,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( dPdx_ref,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (  nstep_init,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (      nsteps,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (       nsave,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (      nstats,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    nmonitor,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   init_type,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   grid_type,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   body_type,len(body_type),MPI_character,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  x_mass_cte,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  y_mass_cte,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( min_buffer_width,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( cg_max_iter,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( cg_tol,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (   body_param_1,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   body_param_2,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   body_param_3,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

  End Subroutine read_input_parameters

  !--------------------------------------!
  !    Generates a grid for channel      !
  !                                      !
  ! Output: x_global, y_global, z_global !
  !                                      !
  !--------------------------------------!
  Subroutine create_grid

    Integer(Int32) :: i

    n_uniform = 1

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
        y_global = 2d0 * y_global / Maxval(y_global)

      Case (1) ! Stretched grid wall to wall
        If ( myid==0 ) Write(*,*) 'Generating stretched y grid'
        Do i=1,ny_global
          y_global(i) = Real(i-1,8)
        End Do
        y_global = 2d0 * y_global / Maxval(y_global) - 1d0

        If ( alpha_stretch > 0d0 ) Then
          Do i=1,ny_global
            y_global(i) = dtanh(alpha_stretch * y_global(i)) / dtanh(alpha_stretch)
          End Do
        End If

        y_global = y_global - Minval(y_global)
        y_global = y_global * 2d0 / Maxval(y_global)

      Case (2) ! Stretched grid centered around 0 with uniform buffers on each end to account for IB (moving with amplitude body_param_1)
        If ( myid==0 ) Write(*,*) 'Generating stretched y grid with uniform buffers'

        ! Loop until buffer condition is met
        write(*,*) 'creating stretched grid using ny_global = ', ny_global
        Do
          Call create_stretched_grid(y_global, 2d0, ny_global, n_uniform, alpha_stretch)

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
        y_global = y_global + 1d0

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
    Real(Int64) :: ym_val

    Select Case (init_type)
      Case (0) ! read input data from file
        If ( myid==0 ) Write(*,*) 'Reading input data'
        Call read_input_data
      
      Case (1) ! create grid and initialize velocity to zero
        If ( myid==0 ) Write(*,*) 'Generating zero initial condition'
        Call create_grid
        U = 0d0; V = 0d0; W = 0d0

      Case (2) ! create grid and initialize velocity to the laminar parabolic profile with random perturbations
        If ( myid==0 ) Write(*,*) 'Generating random initial condition'
        Call create_grid

        ! U
        Do jj=1,ny_global-1
          ym_val = 0.5d0 * (y_global(jj) + y_global(jj + 1))
          U(:,jj+1,:) = dpdx / (2d0 * nu) * ym_val * (2d0 - ym_val)
        end Do
        Do ii=1,nx_global
          Do jj=1,nyg_global
             Do kk=1,nzg
               U(ii,jj,kk) = U(ii,jj,kk) + 0.5*(rand()-0.5)
             End Do
          End Do
        End Do
        U(:,1,:) = -U(:,2,:)
        U(:,ny_global+1,:) = -U(:,ny_global,:)

        ! V
        V = 0d0
        Do ii=1,nxg_global
          Do jj=1,ny_global
             Do kk=1,nzg
               V(ii,jj,kk) = V(ii,jj,kk) + 0.5*(rand()-0.5)
             End Do
          End Do
        End Do
        V(:,1,:) = 0d0
        V(:,ny_global,:) = 0d0

        ! W
        W = 0d0
        Do ii=1,nxg_global
          Do jj=1,ny_global
             Do kk=1,nz
               W(ii,jj,kk) = W(ii,jj,kk) + 0.5*(rand()-0.5)
             End Do
          End Do
        End Do
        W(:,1,:) = -W(:,2,:)
        W(:,ny_global+1,:) = -W(:,ny_global,:)
    End Select

    If ( myid==0 ) Then
      Write(*,*) 'Max U',MaxVal(U)
      Write(*,*) 'Max V',MaxVal(V)
      Write(*,*) 'Max W',MaxVal(W)

      Write(*,*) 'Mean U',sum(U)/Real(nx_global*nyg_global*nzg_global,8)
      Write(*,*) 'Mean V',sum(V)/Real(nxg_global*ny_global*nzg_global,8)
      Write(*,*) 'Mean W',sum(W)/Real(nxg_global*nyg_global*nz_global,8)
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

        ! surface stress
        Write(1) Shape(fb, Int32), fb

        ! interpolated body velocity
        Write(1) Shape(ub, Int32), ub

        ! body surface areas
        Write(1) Shape(sb, Int32), sb

        ! body normals and tangents
        Write(1) Shape(normals, Int32), normals
        Write(1) Shape(tangents_1, Int32), tangents_1
        Write(1) Shape(tangents_2, Int32), tangents_2

      End If
         
      ! close file
      If (myid==0) Then
        Close(1)
      End If
      
    End If
       
  End Subroutine output_data

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
       Write(3,'(A,4F15.8,4I)') '%',t, Retau, utau, nu, nx_global, ny_global, nz_global, istep
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

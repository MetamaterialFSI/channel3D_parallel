!------------------------------------------------!
! Module with initialization of global variables !
!------------------------------------------------!
Module initialization

  ! Modules
  Use, Intrinsic :: iso_c_binding
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use input_output
  Use mass_flow
  Use immersed_boundary_geometry

  ! prevent implicit typing
  Implicit None

  ! declarations
Contains

  !----------------------------------------!
  !      Initialize flow variables         !
  !----------------------------------------!
  Subroutine initialize
    
    Integer(Int32) :: i, j, k, kk, nzpe, pos, ipos, nze, nzme, count
    Real   (Int64) :: dy1, dy2, det, a, b, c, r, Qflow_ref
    Integer(Int32), Dimension(:,:), Allocatable :: A_kmodes, A_kmodes_local

    If (myid==0) Then 
       Write(*,*) '----------------------------------------------------------------------'       
       Write(*,*) '                                                                      '
       Write(*,*) '           My channel ^^, parallel version (clean) 1.0                '
       Write(*,*) '                                                                      '
       Write(*,*) '----------------------------------------------------------------------'
    End If

    ! time
    t = t_init
    
    !-------------------grid definitions-------------------------!
    Allocate (  k1_global(0:nprocs-1),  k2_global(0:nprocs-1) )
    Allocate ( kg1_global(0:nprocs-1), kg2_global(0:nprocs-1) )

    ! restrictions for FFTW mapping
    If ( Mod( nx_global   , 2     )/=0 ) Stop 'Error: nx must be even'
    If ( Mod( nz_global   , 2     )/=0 ) Stop 'Error: nz must be even'
    If ( Mod( nz_global-2 , nprocs)/=0 ) Stop 'nz-2 should be divisible by number of processors'

    ! number of interior z-planes per processor based on fftw decomposition
    nslices_z = Nint( Real((nz_global-2))/Real(nprocs) ) 
    
    ! restriction for MPI boundaries
    If ( nslices_z<2 ) Stop 'Error: nslices_z must be at least 2' 

    ! domain decomposition. Must be consistent with fftw
    Do i = 0, nprocs-1
       ! range index for faces in each processor
       k1_global(i)  = i*nslices_z  + 1
       k2_global(i)  = k1_global(i) + nslices_z + 1
       ! range index for centers in each processor
       kg1_global(i) = i*nslices_z   + 1
       kg2_global(i) = kg1_global(i) + nslices_z + 1
    End Do    

    ! remaining planes in last processor
    k2_global (nprocs-1) = nz_global 
    kg2_global(nprocs-1) = nz_global + 1

    ! face points
    nx = nx_global
    ny = ny_global
    nz = k2_global(myid) - k1_global(myid) + 1 

    ! middle points
    nxm_global = nx_global - 1
    nym_global = ny_global - 1
    nzm_global = nz_global - 1

    nxm = nx - 1
    nym = ny - 1
    nzm = kg2_global(myid) - kg1_global(myid) + 1 - 2  
    
    ! middle points + ghost cells
    nxg_global = nxm_global + 2
    nyg_global = nym_global + 2
    nzg_global = nzm_global + 2

    nxg = nxm + 2
    nyg = nym + 2
    nzg = kg2_global(myid) - kg1_global(myid) + 1 

    ! size for last proccesor nz and nzm -> nze and nzme
    nze  = nz
    nzme = nzm
    Call Mpi_bcast (  nze,1,MPI_integer,nprocs-1,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nzme,1,MPI_integer,nprocs-1,MPI_COMM_WORLD,ierr )

    ! Allocate main arrays
    If ( myid==0 ) Write(*,*) 'allocating main arrays...'
    Allocate ( x_global (  nx_global),  y_global (  ny_global),  z_global (  nz_global)  )
    Allocate ( xm_global( nxm_global),  ym_global( nym_global),  zm_global( nzm_global)  )
    Allocate ( xg_global(nxm_global+2), yg_global(nym_global+2), zg_global(nzm_global+2) )

    Allocate (  x (  nx),  y (  ny),  z (  nz) )
    Allocate (  xm( nxm),  ym( nym),  zm( nzm) )
    Allocate ( xg(nxm+2), yg(nym+2), zg(nzm+2) )

    Allocate ( yg_m (nyg-1) )
    Allocate ( yg_mm(nyg-2) )

    ! global interior + boundary + ghost points
    Allocate (U (    nx, nym+2, nzm+2) )
    Allocate (V ( nxm+2,    ny, nzm+2) )
    Allocate (W ( nxm+2, nym+2,    nz) )
    Allocate (P ( nxm+2, nym+2, nzm+2) )

    Allocate (Uo  (    nx,  nym+2, nzm+2) )
    Allocate (Vo  ( nxm+2,     ny, nzm+2) )
    Allocate (Wo  ( nxm+2,  nym+2,    nz) )
    Allocate (Po  ( nxm+2,  nym+2, nzm+2) )

    If (myid == 0) Then
       Allocate (Uoo (    nx,  nym+2, nzme+2) ) ! z-planes modified for I/O
       Allocate (Voo ( nxm+2,     ny, nzme+2) )
       Allocate (Woo ( nxm+2,  nym+2,    nze) )
       Allocate (Poo ( nxm+2,  nym+2, nzme+2) )
    End If

    ! arrays for accepting regularized body distributions
    Allocate (U_reg  (    nx,  nym+2, nzm+2) )
    Allocate (V_reg  ( nxm+2,     ny, nzm+2) )
    Allocate (W_reg  ( nxm+2,  nym+2,    nz) )

    ! arrays for storing intermediate velocity fields in the projection method
    Allocate (U_interim  (    nx,  nym+2, nzm+2) )
    Allocate (V_interim  ( nxm+2,     ny, nzm+2) )
    Allocate (W_interim  ( nxm+2,  nym+2,    nz) )
    Allocate (P_interim  ( nxm+2,  nym+2, nzm+2) )

    Allocate (Vw ( nxm+2, 2, nzm+2) )

    ! Auxiliary arrays
    Allocate ( term_1 ( nxg, nyg, nzg ) ) 
    Allocate ( term_2 ( nxg, nyg, nzg ) ) 

    ! RHS: interior points only
    Allocate ( rhs_uo ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) ) 
    Allocate ( rhs_vo ( 2:nxg-1, 2:ny-1,  2:nzg-1 ) )
    Allocate ( rhs_wo ( 2:nxg-1, 2:nyg-1, 2:nz-1  ) )
    Allocate ( rhs_p  ( 2:nxg-1, 2:nyg-1, 2:nzg   ) ) ! ONE EXTRA PLANE IN Z FOR GHOST CELL
    If (.False.) Then ! used to calculate pressure from eqs.
       Allocate ( rhs_uf  ( 1:nx  ,  1:nyg  , 1:nzg   ) ) 
       Allocate ( rhs_vf  ( 1:nxg  , 1:ny  ,  1:nzg   ) )
       Allocate ( rhs_wf  ( 1:nxg  , 1:nyg  , 1:nz    ) )
    End If

    ! read data 
    If ( myid==0 ) Write(*,*) 'preparing initial condition...'
    Call init_flow

    ! define global grids from x_global, y_global and z_global (face to centers)
    ! local faces
    x = x_global
    y = y_global
    z = z_global( k1_global(myid):k2_global(myid) )

    ! global interior centers
    Do i = 1, nxm_global
       xm_global(i) = 0.5d0*( x_global(i) + x_global(i+1) )
    End Do
    Do j=1,nym_global
       ym_global(j) = 0.5d0*( y_global(j) + y_global(j+1) )
    End Do
    Do k=1,nzm_global
       zm_global(k) = 0.5d0*( z_global(k) + z_global(k+1) )
    End Do

    ! local interior centers
    xm = xm_global
    ym = ym_global
    zm = zm_global( kg1_global(myid):kg2_global(myid)-2 )

    ! global 
    xg_global(2:nxm_global+1) = xm_global
    xg_global(1)              = xm_global(1)          - 2d0*(xm_global(1)-x_global(1))
    xg_global(nxm_global+2)   = xm_global(nxm_global) + 2d0*(x_global(nx_global)-xm_global(nxm_global))
    
    yg_global(2:nym_global+1) = ym_global
    yg_global(1)              = ym_global(1)          - 2d0*(ym_global(1)-y_global(1))
    yg_global(nym_global+2)   = ym_global(nym_global) + 2d0*(y_global(ny_global)-ym_global(nym_global))
    
    zg_global(2:nzm_global+1) = zm_global
    zg_global(1)              = zm_global(1)          - 2d0*(zm_global(1)-z_global(1))
    zg_global(nzm_global+2)   = zm_global(nzm_global) + 2d0*(z_global(nz_global)-zm_global(nzm_global))    

    xg = xg_global
    yg = yg_global
    zg = zg_global( kg1_global(myid):kg2_global(myid) )

    ! middle points for yg (.not. equal to y in general)
    yg_m = 0.5d0*( yg(2:nyg) + yg(1:nyg-1) )

    ! middle points for yg_m (.not. equal to ym in general)
    yg_mm = 0.5d0*( yg_m(2:nyg-1) + yg_m(1:nyg-2) )

    ! local minimum grid size for CFL
    dxmin = Minval ( xg_global(2:nxg_global) - xg_global(1:nxg_global-1) )
    dymin = Minval ( yg_global(2:nyg_global) - yg_global(1:nyg_global-1) )
    dzmin = Minval ( zg_global(2:nzg_global) - zg_global(1:nzg_global-1) )

    ! total domain size
    Lx = x_global(nx_global) - x_global(1)
    Ly = y_global(ny_global) - y_global(1)
    Lz = z_global(nz_global) - z_global(1)

    ! For initial IB implementation only!
    If ( body_type > 0) Then
      Allocate( U_global(nx_global,  nyg_global, nzg_global) )
      Allocate( V_global(nxg_global, ny_global,  nzg_global) )
      Allocate( W_global(nxg_global, nyg_global, nz_global ) )
      Allocate( send_counts_U(nprocs), displs_U(nprocs) )
      Allocate( send_counts_V(nprocs), displs_V(nprocs) )
      Allocate( send_counts_W(nprocs), displs_W(nprocs) )

      local_size_U = nx * nyg * nzm
      local_size_V = nxg * ny * nzm
      local_size_W = nxg * nyg * (nz-2)

      ! Gather send_counts
      Call MPI_Allgather(local_size_U, 1, MPI_INT, send_counts_U, 1, MPI_INT, MPI_COMM_WORLD, ierr)
      Call MPI_Allgather(local_size_V, 1, MPI_INT, send_counts_V, 1, MPI_INT, MPI_COMM_WORLD, ierr)
      Call MPI_Allgather(local_size_W, 1, MPI_INT, send_counts_W, 1, MPI_INT, MPI_COMM_WORLD, ierr)

      displs_U(1) = 0
      displs_V(1) = 0
      displs_W(1) = 0
      Do i = 2, nprocs
        displs_U(i) = displs_U(i-1) + send_counts_U(i-1)
        displs_V(i) = displs_V(i-1) + send_counts_V(i-1)
        displs_W(i) = displs_W(i-1) + send_counts_W(i-1)
      End Do
    End If

    !--------------------------Boundary conditions--------------------------!
    ! local velocity, initial z-planes
    Allocate ( buffer_ui(nx,nyg,2:3), buffer_vi(nxg,ny,2:3), buffer_wi(nxg,nyg), buffer_ci(nxg,nyg,2:3) )
    ! local velocity, ending  z-planes
    Allocate ( buffer_ue(nx,nyg),     buffer_ve(nxg,ny),     buffer_we(nxg,nyg), buffer_ce(nxg,nyg) )
    ! local pressure z-plane
    Allocate ( buffer_p(2:nxg-1,2:nyg-1) ) 

    !------------------------Interior communications------------------------!
    Allocate ( buffer_us(nx ,nyg), buffer_ur(nx ,nyg) )
    Allocate ( buffer_vs(nxg, ny), buffer_vr(nxg, ny) )
    Allocate ( buffer_ws(nxg,nyg), buffer_wr(nxg,nyg) )
    Allocate ( buffer_ps(2:nxg-1,2:nyg-1), buffer_pr(2:nxg-1,2:nyg-1) ) 

    !---------------------------Fourier transform---------------------------!
    If ( myid==0 ) Write(*,*) 'initializing FFT...'
    ! initialize MPI FFTW
    Call fftw_mpi_init()

    ! Fourier constant grid spacing
    dx = dxmin
    dz = dzmin

    ! length for periodic domain
    Lxp = Lx - dx 
    Lzp = Lz - dz

    ! global points for periodic domain in physical space
    nxp_global = nxm_global - 1
    nzp_global = nzm_global - 1

    ! global indices for fourier modes starting from 0
    mx_global = nxp_global - 1
    mz_global = nzp_global - 1

    ! Get local sizes:
    ! local data size in x direction
    nxp = nxp_global
    mx  =  mx_global
    ! local data size in z direction (note dimension reversal)
    alloc_local = fftw_mpi_local_size_2d(nzp_global, nxp_global, MPI_COMM_WORLD, nzp, local_k_offset)
    mz  = nzp - 1

    ! sanity check and restrictions in fftw
    If ( (nzp/=nzm .And. myid/=nprocs-1) .Or. (nzp/=nzm-1 .And. myid==nprocs-1) ) Then 
       Write(*,*) nzp,nzm
       Stop 'Error: something wrong in FFTW size'
    End If

    ! allocate variables
    cplane_fft = fftw_alloc_complex(alloc_local)
    Call c_f_pointer(cplane_fft,plane,[nxp,nzp])
    plane_hat(0:,0:) => plane
    Allocate ( rhs_hat ( 0:mx, 2:nyg-1, 0:mz ) )
    Allocate ( rhs_aux   ( 2:nyg-1 ) )
   
    ! create MPI plan for forward DFT (note dimension reversal and transposed_out/in)
    plan_d = fftw_mpi_plan_dft_2d( nzp_global, nxp_global, plane, plane_hat,           & 
             MPI_COMM_WORLD,  FFTW_FORWARD, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT) ) 
    plan_i = fftw_mpi_plan_dft_2d( nzp_global, nxp_global, plane_hat, plane,           & 
             MPI_COMM_WORLD, FFTW_BACKWARD, ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN)  ) 

    ! global Fourier coeficients with modified wave-number for the second derivative
    Allocate ( kxx(0:mx_global), kzz(0:mz_global) ) 
    kxx = 0d0
    kzz = 0d0
    Do i = 0, Ceiling( Real(nxp_global)/2d0 )
       kxx(i) = 2d0*( dcos(2d0*pi*Real(i,8)/Real(nxp_global,8)) - 1d0 )/dx**2d0  
    End do
    Do i = Ceiling( Real(nxp_global)/2d0 )+1, mx_global
       kxx(i) = 2d0*( dcos(2d0*pi*Real(-nxp_global+i,8)/Real(nxp_global,8)) - 1d0 )/dx**2d0
    End do

    Do k = 0, Ceiling( Real(nzp_global)/2d0 )
       kzz(k) = 2d0*( dcos(2d0*pi*Real(k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0  
    End do
    Do k = Ceiling( Real(nzp_global)/2d0 )+1, mz_global
       kzz(k) = 2d0*( dcos(2d0*pi*Real(-nzp_global+k,8)/Real(nzp_global,8)) - 1d0 )/dz**2d0
    End do
    
    ! MPI mapping for z-modes: from local to global without transposed_out/in
    ! Not used in this version
    Allocate ( kmode_map(0:mz) ) 
    kmode_map = 0
    Do k = 0, mz
       kmode_map(k) = k + myid*nslices_z
    End Do

    ! FFTW+MPI mapping for x and z-modes when using FFTW with transposed_out/in
    ! from local to global
    ! this needs (mz_global+1)*(mx_global+1)/nprocs to be an integer 
    If ( Mod((mz_global+1)*(mx_global+1),nprocs)/=0 ) Stop 'Error: (mz_global+1)*(mx_global+1)/nprocs should be an integer'
    Allocate ( imode_map_fft(0:mx_global,0:mz) ) 
    Allocate ( kmode_map_fft(0:mx_global,0:mz) ) 
    Do i = 0, mx_global
       Do k = 0, mz            
          pos = i + (mx_global+1)*k + (mz_global+1)*(mx_global+1)/nprocs*myid
          imode_map_fft(i,k) = Floor( Real(pos/(mz_global+1)) )
          kmode_map_fft(i,k) = Mod  ( pos, mz_global+1 )
          ! sanity check
       end Do
    End Do

    ! Sanity check for FFTW mapping
    Allocate(A_kmodes      (0:mx_global,0:mz_global))
    Allocate(A_kmodes_local(0:mx_global,0:mz_global))
    A_kmodes       = 0
    A_kmodes_local = 0
    Do i = 0, mx_global
       Do k = 0, mz            
          A_kmodes_local( imode_map_fft(i,k), kmode_map_fft(i,k) ) =  A_kmodes_local(imode_map_fft(i,k), kmode_map_fft(i,k) ) + 1
       end Do
    End Do
    count = (mx_global+1)*(mz_global+1)
    Call MPI_AllReduce(A_kmodes_local,A_kmodes,count,MPI_integer,MPI_sum,MPI_COMM_WORLD,ierr)
    If ( Any(A_kmodes>1) .Or. Any(A_kmodes==0) ) Stop 'Error: wrong combination of nx, nz and processors'
    Deallocate(A_kmodes)
    Deallocate(A_kmodes_local)
     
    !------------------------Tridiagonal linear solver-------------------------!
    If ( myid==0 ) Write(*,*) 'initializing pressure solver...'
    Allocate ( pivot(nyg) )
    Allocate ( Dyy(2:nyg-1,2:nyg-1), M(2:nyg-1,2:nyg-1) )
    Allocate ( D(2:nyg-1), DL(2:nyg-2), DU(2:nyg-2) )
     
    ! second derivative matrix for pressure (full data in y assumed)
    Dyy = 0d0
    Do j=3,nyg-2

       a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
       b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
       c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 

       Dyy(j,j+1) = a
       Dyy(j,j-1) = c 
       Dyy(j,j  ) = b

    End Do

    ! Boundary conditions for pressure (full data in y assumed)
    j = 2
    a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) ) 
    ! Dirichlet in V: p(1)==p(2) 
    Dyy(2,2)   = b + c 
    Dyy(2,3)   = a
    coef_bc_1  = c

    j = nyg-1
    a = 1d0/( y(j)-y(j-1) )/( yg(j+1) - yg(j) )
    b = 1d0/( y(j)-y(j-1) )*( -1d0/( yg(j+1) - yg(j) ) -1d0/( yg(j) - yg(j-1) ) )
    c = 1d0/( y(j)-y(j-1) )/( yg(j) - yg(j-1) )     
    ! Dirichlet in V: p(nyg)==p(nyg-1) 
    Dyy(nyg-1,nyg-1) = a + b
    Dyy(nyg-1,nyg-2) = c
    coef_bc_2        = a

    Allocate ( bc_1(2:nxg-1,2:nzg-1), bc_2(2:nxg-1,2:nzg-1) )
    Allocate ( bc_1_hat(0:mx,0:mz),   bc_2_hat(0:mx,0:mz)   )

    ! some parameters for linear solver
    nr   = nym
    nrhs = 1


    !------------------------statistics---------------------------!
    Allocate (  Umean(nyg),  Vmean(nyg),  Wmean(nyg)              )
    Allocate ( U2mean(nyg), V2mean(nyg), W2mean(nyg), UVmean(nyg) )
    Umean  = 0d0
    Vmean  = 0d0
    Wmean  = 0d0
    U2mean = 0d0
    V2mean = 0d0
    W2mean = 0d0
    UVmean = 0d0

    !------------------------Runge-Kutta 3-------------------------!
    If ( myid==0 ) Write(*,*) 'initializing time integration...'
    Allocate( rk_coef(3,3), rk_t(3) )

    rk_t(1)      =  8d0/15d0
    rk_t(2)      =  2d0/3d0
    rk_t(3)      =  1d0

    rk_coef      =  0d0
    rk_coef(1,1) =  8d0/15d0
    rk_coef(2,1) =  1d0/4d0
    rk_coef(2,2) =  5d0/12d0
    rk_coef(3,1) =  1d0/4d0
    rk_coef(3,2) =  0d0
    rk_coef(3,3) =  3d0/4d0

    Allocate ( Fu1 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    Allocate ( Fu2 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )
    Allocate ( Fu3 ( 2:nx-1,  2:nyg-1, 2:nzg-1 ) )

    Allocate ( Fv1 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    Allocate ( Fv2 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )
    Allocate ( Fv3 ( 2:nxg-1,  2:ny-1, 2:nzg-1 ) )

    Allocate ( Fw1 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    Allocate ( Fw2 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )
    Allocate ( Fw3 ( 2:nxg-1,  2:nyg-1, 2:nz-1 ) )

    !-------------------compute initial mass flow-----------------!    
    Call compute_mean_mass_flow_U(U,Qflow_x_0)
    Call compute_mean_mass_flow_V(V,Qflow_y_0)
    Qflow_y_0 = 0d0
    dPdy      = 0d0


    !-------------------------Done--------------------------------!
    Call Mpi_barrier(MPI_COMM_WORLD,ierr)

    ! Measure time
    time1 = MPI_WTIME()
    
  End Subroutine initialize

  Subroutine initialize_ib_arrays
    !--------------------Initialize main arrays-------------------!    
    If ( myid==0 ) Write(*,*) 'allocating main IB arrays...'
    ! Number of body points
    Call compute_nb

    ! Body coordinates
    Allocate ( xb (  nb), yb (  nb), zb (  nb) )

    ! Wall normal reference coordinate and index
    Allocate ( y_ref_index (  nb) )

    ! Body areas
    Allocate(sb (nb) )

    ! Body forcing
    Allocate (fb(3*nb) )
    fb = 0d0

    ! Body velocity
    Allocate (ub (3 * nb) )
    ub = 0d0

    ! Body normals and tangents
    Allocate (normals (3 * nb) )
    Allocate (tangents_1 (3 * nb) )
    Allocate (tangents_2 (3 * nb) )
    normals= 0d0
    tangents_1= 0d0
    tangents_2= 0d0

    ! Auxiliary surface arrays
    Allocate ( rhs_ib (3 * nb) )
    Allocate ( aux_surface_vector (3 * nb) )
    Allocate ( aux_surface_scalar (nb) )
    rhs_ib = 0d0
    aux_surface_vector = 0d0
    aux_surface_scalar = 0d0
    Allocate ( Fibu (nx,nyg,nzg) )
    Allocate ( Fibv (nxg,ny,nzg) )
    Allocate ( Fibw (nxg,nyg,nz) )

    !--------------------Initialize IB operator variables-------------------!    

    Allocate ( x_pivot_index  (nb) )
    Allocate ( xm_pivot_index (nb) )
    Allocate ( y_pivot_index  (nb) )
    Allocate ( ym_pivot_index (nb) )
    Allocate ( z_pivot_index  (nb) )
    Allocate ( zm_pivot_index (nb) )
    
    Allocate ( u_weights  ( nweights, nb) )
    Allocate ( u_x_indices( nweights, nb) )
    Allocate ( u_y_indices( nweights, nb) )
    Allocate ( u_z_indices( nweights, nb) )
    
    Allocate ( v_weights  ( nweights, nb) )
    Allocate ( v_x_indices( nweights, nb) )
    Allocate ( v_y_indices( nweights, nb) )
    Allocate ( v_z_indices( nweights, nb) )
    
    Allocate ( w_weights  ( nweights, nb) )
    Allocate ( w_x_indices( nweights, nb) )
    Allocate ( w_y_indices( nweights, nb) )
    Allocate ( w_z_indices( nweights, nb) )

    Allocate( send_counts_nb(nprocs), displs_nb(nprocs) )

    !-------------------------Done--------------------------------!
    Call Mpi_barrier(MPI_COMM_WORLD,ierr)

  End Subroutine
  
End Module initialization

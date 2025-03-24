!-----------------------------------------!
! Module with all shared global variables !
!-----------------------------------------!
Module global

  ! General Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use, Intrinsic :: iso_c_binding
  Use mpi

  ! prevent implicit typing
  Implicit None

  ! FFTW
  Include 'fftw3-mpi.f03'  

  !------------------Declarations-----------------!

  ! declarations
  Integer(Int32) :: ierr, myid, nprocs, nslices_z
  Integer        :: istat ( MPI_STATUS_SIZE )

  ! global faces index range for each processor
  Integer(Int32), Dimension(:), Allocatable ::  k1_global,  k2_global

  ! global centers index range for each processor
  Integer(Int32), Dimension(:), Allocatable :: kg1_global, kg2_global

  ! step number
  Integer(Int32) :: istep, rk_step
  Real   (Int64) :: time1, time2

  ! constants
  Real(Int64) :: pi = 4d0*datan(1d0)  

  ! files
  Character(200) :: filein, fileout, fileparams
  Integer(Int32) :: nsave, nmonitor

  ! initial condition
  Integer(Int32) :: init_type

  ! grid type
  Integer(Int32) :: grid_type

  ! domain size
  Real(Int64) :: Lx, Lz, Ly, Lxp, Lzp

  ! steps
  Integer(Int32) :: nsteps, nstep_init
  Real   (Int64) :: dt, t

  ! viscosity
  Real(Int64) :: nu

  ! global face points
  Integer(Int32) :: nx_global, ny_global, nz_global

  ! local face points
  Integer(Int32) :: nx, ny, nz

  ! global center points
  Integer(Int32) :: nxm_global, nym_global, nzm_global

  ! local center points
  Integer(Int32) :: nxm, nym, nzm

  ! global center points + ghost cells
  Integer(Int32) :: nxg_global, nyg_global, nzg_global

  ! local center points + ghost cells
  Integer(Int32) :: nxg, nyg, nzg

  ! global grid at face points
  Real(Int64), Allocatable, Dimension(:) :: x_global, y_global, z_global

  ! local grid at face points
  Real(Int64), Allocatable, Dimension(:) :: x, y, z

  ! global grid at middle points
  Real(Int64), Allocatable, Dimension(:) :: xm_global, ym_global, zm_global

  ! local grid at middle points
  Real(Int64), Allocatable, Dimension(:) :: xm, ym, zm

  ! global grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg_global, yg_global, zg_global

  ! local grid at middle points + ghost cells
  Real(Int64), Allocatable, Dimension(:) :: xg, yg, zg

  ! middle points for yg->yg_m and yg_m->yg_mm
  Real(Int64), Allocatable, Dimension(:) :: yg_m, yg_mm

  ! local velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:,:) :: U,V,W,P
  Real(Int64), Allocatable, Dimension(:,:,:) :: Uo,Vo,Wo,Po
  Real(Int64), Allocatable, Dimension(:,:,:) :: Uoo,Voo,Woo,Poo
  Real(Int64), Allocatable, Dimension(:,:,:) :: Vw 

  ! local auxiliary 
  Real(Int64), Allocatable, Dimension(:,:,:) :: term_1, term_2

  ! local rhs for velocities and pressure
  Real(Int64), Allocatable, Dimension(:,:)   :: px_bottom, px_top
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_p
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_uo, rhs_vo, rhs_wo
  Real(Int64), Allocatable, Dimension(:,:,:) :: rhs_uf, rhs_vf, rhs_wf

  ! local rhs for pressure in Fourier
  Complex(Int64), Dimension(:,:,:), Allocatable :: rhs_p_hat
  Complex(Int64), Dimension(:),     Allocatable :: rhs_aux

  ! local auxiliary arrays for MPI_sendrev boundary conditions
  Real(Int64), Allocatable, Dimension(:,:,:) :: buffer_ui, buffer_vi, buffer_ci
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_ue, buffer_ve, buffer_we, buffer_wi, buffer_ce
  Real(Int64), Allocatable, Dimension(:,:)   :: buffer_p

  ! local auxiliary arrays for MPI_sendrev interior planes
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_us, buffer_ur
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_vs, buffer_vr
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_ws, buffer_wr
  Real(Int64), Allocatable, Dimension(:,:) :: buffer_ps, buffer_pr
  
  ! local auxiliary planes for FFTW
  Type(C_PTR) :: cplane_fft
  Complex(C_DOUBLE_COMPLEX), Pointer, Dimension(:,:) :: plane, plane_hat

  ! Fourier points and wave numbers 
  Integer(C_INTPTR_T) :: nxp_global, nzp_global, local_k_offset
  Integer(C_INTPTR_T) :: nxp, nzp
  Integer(C_INTPTR_T) :: mx_global, mz_global
  Integer(C_INTPTR_T) :: mx, mz
  Real   (Int64)      :: dx, dz
  Real   (Int64), Dimension(:), Allocatable :: kxx, kzz

  ! Mappings for fft modes
  Integer(Int64), Dimension(:),   Allocatable :: kmode_map
  Integer(Int64), Dimension(:,:), Allocatable :: imode_map_fft, kmode_map_fft
  
  ! FFTW plans
  Integer(C_INTPTR_T) :: alloc_local
  Type   (C_PTR)      :: plan_d, plan_i

  ! finite differences (second derivative)
  Real(Int64) :: ddx1, ddx2, ddx3
  Real(Int64) :: ddy1, ddy2, ddy3
  Real(Int64) :: ddz1, ddz2, ddz3

  ! linear solver
  Integer (Int32) :: nr, nrhs
  Integer (Int32), Dimension(:),   Allocatable :: pivot  
  Complex (Int64), Dimension(:),   Allocatable :: D, DL, DU
  Complex (Int64), Dimension(:,:), Allocatable :: M, Dyy

  ! pressure gradients
  Real(Int64) :: dPdx, dPdy, dPdz, dPdx_ref

  ! constant mass flow
  Real   (Int64) :: Qflow_x_0, Qflow_y_0
  Integer(Int32) :: x_mass_cte, y_mass_cte
    
  ! CFL parameters
  Real(Int64) :: CFL, dxmin, dymin, dzmin

  ! actual pressure boundary conditions
  Real   (Int64) :: coef_bc_1, coef_bc_2
  Real   (Int64), Dimension(:,:), Allocatable :: bc_1,     bc_2
  Complex(Int64), Dimension(:,:), Allocatable :: bc_1_hat, bc_2_hat
  Logical(Int32) :: pressure_computed

  ! statistics
  Integer(Int32) :: nstats, Retau_int
  Real   (Int64) :: Retau, utau, Qflow_x, Qflow_y
  Real   (Int64), Dimension(:), Allocatable ::  Umean,  Vmean,  Wmean
  Real   (Int64), Dimension(:), Allocatable :: U2mean, V2mean, W2mean, UVmean

  ! Runge-Kutta 3 coefficients and buffers
  Real(Int64), Dimension(:),     Allocatable :: rk_t
  Real(Int64), Dimension(:,:),   Allocatable :: rk_coef
  Real(Int64), Dimension(:,:,:), Allocatable :: Fu1, Fu2, Fu3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fv1, Fv2, Fv3
  Real(Int64), Dimension(:,:,:), Allocatable :: Fw1, Fw2, Fw3

  ! body mode
  Integer(Int32) :: body_type

  ! number of uniform grid points on each side of the IB
  Integer(Int32) :: nd

  ! body points
  Integer(Int32) :: nb, nxb, nzb
  Real   (Int64) :: dxb, dzb
  Real   (Int64), Dimension(:), Allocatable :: xb, yb, zb

  ! body motion
  Real   (Int64) :: amp, omega
  Integer(Int32) :: wave_nx

  ! body surface areas
  Real(Int64), Dimension(:), Allocatable :: sb ! body face areas interpolated to body nodes

  ! body velocity
  Real(Int64), Dimension(:), Allocatable :: ub

  ! body normals
  Real(Int64), Dimension(:), Allocatable :: normals, tangents_1, tangents_2

  !immersed body forcing
  Real(Int64), Dimension(:), Allocatable :: fb

  ! regularization and interpolation support, weights, and indices
  Integer(Int32) :: suppx, suppy, suppz, nweights
  Real(Int64),    Dimension(:,:),   Allocatable :: u_weights, v_weights, w_weights
  Integer(Int32), Dimension(:,:),   Allocatable :: u_x_indices, u_y_indices, u_z_indices
  Integer(Int32), Dimension(:,:),   Allocatable :: v_x_indices, v_y_indices, v_z_indices
  Integer(Int32), Dimension(:,:),   Allocatable :: w_x_indices, w_y_indices, w_z_indices
  Integer(Int32), Dimension(:),     Allocatable :: x_pivot_index, xm_pivot_index
  Integer(Int32), Dimension(:),     Allocatable :: y_pivot_index, ym_pivot_index
  Integer(Int32), Dimension(:),     Allocatable :: z_pivot_index, zm_pivot_index
  
End Module global

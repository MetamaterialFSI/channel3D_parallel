Module operators_FSI

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use immersed_boundary_operators

  Implicit None

Contains

function Itilde(x)

  REAL(KIND(0.d0)), DIMENSION(nb), INTENT(in) :: x
  real(kind(0.d0)), dimension(3*nb) :: Itilde
  real(kind(0.d0)), dimension(3*nb) :: v_patch_3d
  integer :: i

  v_patch_3d = 0.d0

  do i = 1, nb
    v_patch_3d(nb+i) = x(i)
  end do

  Itilde = v_patch_3d

end function Itilde


function KhatinvQItildeprimeW(x)

  REAL(KIND(0.d0)), DIMENSION(3*nb), INTENT(in) :: x
  real(kind(0.d0)), dimension(nb) :: KhatinvQItildeprimeW
  real(kind(0.d0)), dimension(3*nb) :: v_tp
  real(kind(0.d0)), dimension(nb) :: v_y, v_trunc, v_bg

  v_tp = redistribute(dt*x)
  !v_tp=x
  v_y  = v_tp(nb+1:2*nb)

  v_trunc = v_y
  v_bg = matmul(sol_mat, v_trunc)

  KhatinvQItildeprimeW = v_bg

end function KhatinvQItildeprimeW


function KhatinvQItildeprimeW_subsurface(x)

  REAL(KIND(0.d0)), DIMENSION(3*nb), INTENT(in) :: x
  real(kind(0.d0)), dimension(nblocks) :: KhatinvQItildeprimeW_subsurface
  real(kind(0.d0)), dimension(3*nb) :: v_tp
  real(kind(0.d0)), dimension(nb) :: v_y
  real(kind(0.d0)), dimension(nblocks) :: v_trunc, v_bg
  integer :: i
  real(kind(0.d0)) :: sum

  v_tp = redistribute(dt*x)
  !v_tp=x
  v_y  = v_tp(nb+1:2*nb)

  v_trunc = 0.d0
  sum = 0.d0

  do i = 1, nb
    sum = sum + v_y(i)*dxb*dzb
  end do

  v_trunc(1) = sum/(Lxp*Lzp)
  v_bg = matmul(sol_mat, v_trunc)

  KhatinvQItildeprimeW_subsurface = v_bg

end function KhatinvQItildeprimeW_subsurface


subroutine khatmatrix_testcase(sol_mat, Mmat_testcase, Kmat_testcase, Cmat_testcase)

  real(kind(0.d0)), dimension(nb, nb), intent(inout) :: Kmat_testcase
  real(kind(0.d0)), dimension(nb, nb), intent(inout) :: Mmat_testcase
  real(kind(0.d0)), dimension(nb, nb), intent(inout) :: Cmat_testcase
  real(kind(0.d0)), dimension(nb, nb), intent(inout) :: sol_mat
  integer :: info, neqns, lda, lwork
  real(kind(0.d0)), dimension(:), allocatable :: work
  integer, dimension(nb) :: ipiv_bg

  sol_mat = Kmat_testcase + 4.d0/dt_fsi**2.d0*Mmat_testcase + 2.d0/dt_fsi*Cmat_testcase

  info = 0
  lwork = nb**2
  allocate(work(lwork))

  neqns = nb
  lda = nb
  ipiv_bg = 0

  call dgetrf(neqns, lda, sol_mat, lda, ipiv_bg, info)
  call dgetri(neqns, sol_mat, lda, ipiv_bg, work, lwork, info)

  deallocate(work)

end subroutine khatmatrix_testcase


subroutine khatmatrix_subsurface(sol_mat, Mmat, Kmat)

  real(kind(0.d0)), dimension(nblocks, nblocks), intent(inout) :: Kmat
  real(kind(0.d0)), dimension(nblocks, nblocks), intent(inout) :: Mmat
  real(kind(0.d0)), dimension(nblocks, nblocks), intent(inout) :: sol_mat
  integer :: info, neqns, lda, lwork
  real(kind(0.d0)), dimension(:), allocatable :: work
  integer, dimension(nblocks) :: ipiv_bg

  sol_mat = Kmat + 4.d0/dt_fsi**2.d0*Mmat

  info = 0
  lwork = nblocks**2
  allocate(work(lwork))

  neqns = nblocks
  lda = nblocks
  ipiv_bg = 0

  call dgetrf(neqns, lda, sol_mat, lda, ipiv_bg, info)
  call dgetri(neqns, sol_mat, lda, ipiv_bg, work, lwork, info)

  deallocate(work)

end subroutine khatmatrix_subsurface


function b_times_testcase(x)

  REAL(KIND(0.d0)), DIMENSION(3*nb), INTENT(in) :: x
  real(kind(0.d0)), dimension(3*nb) :: b_times_testcase
  real(kind(0.d0)), dimension(3*nb) :: v_tp, v_patch_3d
  real(kind(0.d0)), dimension(nb) :: v_y, v_patch, v_trunc, v_bg
  integer :: i

  v_tp = redistribute(dt*x)
  !v_tp=x
  v_y = 0.d0

  do i = 1, nb
    v_y(i) = v_tp(nb+i)
  end do

  v_trunc = v_y
  v_bg = matmul(sol_mat, v_trunc)

  do i = 1, nb
    v_patch(i) = v_bg(i)
  end do

  v_patch_3d = 0.d0

  do i = 1, nb
    v_patch_3d(nb+i) = (-2.d0/dt_fsi)*v_patch(i)
  end do

  b_times_testcase = v_patch_3d

end function b_times_testcase


function b_times_subsurface(x)

  REAL(KIND(0.d0)), DIMENSION(3*nb), INTENT(in) :: x
  real(kind(0.d0)), dimension(3*nb) :: b_times_subsurface
  real(kind(0.d0)), dimension(3*nb) :: v_tp, v_patch_3d
  real(kind(0.d0)), dimension(nb) :: v_y, v_patch
  real(kind(0.d0)), dimension(nblocks) :: v_trunc, v_bg
  integer :: i
  real(kind(0.d0)) :: sum

  v_tp = redistribute(dt*x)
  !v_tp=dt*x
  v_y = 0.d0

  do i = 1, nb
    v_y(i) = v_tp(nb+i)
  end do

  v_trunc = 0.d0
  sum = 0.d0

  do i = 1, nb
    sum = sum + v_y(i)*dxb*dzb
  end do

  v_trunc(1) = sum/(Lxp*Lzp)
  v_bg = matmul(sol_mat, v_trunc)

  do i = 1, nb
    v_patch(i) = cos((3.d0*pi/Lxp)*(xb(i)-xtopmass)) * &
                 exp(-(xb(i)-xtopmass)**2.d0/(2.d0*sigma**2.d0)) * v_bg(1)
  end do

  v_patch_3d = 0.d0

  do i = 1, nb
    v_patch_3d(nb+i) = (-2.d0/dt_fsi)*v_patch(i)
  end do

  b_times_subsurface = v_patch_3d

end function b_times_subsurface


function redistribute(f_vector)

  real(kind(0.d0)), dimension(3*nb), intent(in) :: f_vector
  real(kind(0.d0)), dimension(3*nb) :: redistribute
  real(kind(0.d0)), dimension(3*nb) :: one_vect
  real(kind(0.d0)), allocatable :: wghtu(:,:,:), frc_regu(:,:,:)
  real(kind(0.d0)), allocatable :: wghtv(:,:,:), frc_regv(:,:,:)
  real(kind(0.d0)), allocatable :: wghtw(:,:,:), frc_regw(:,:,:)
  integer :: i, j, k

  allocate(wghtu(nx,nyg,nzg), frc_regu(nx,nyg,nzg))
  allocate(wghtv(nxg,ny,nzg), frc_regv(nxg,ny,nzg))
  allocate(wghtw(nxg,nyg,nz), frc_regw(nxg,nyg,nz))

  one_vect = 1.d0
  redistribute = 0.d0

  call regu(wghtu, one_vect)
  call regv(wghtv, one_vect)
  call regw(wghtw, one_vect)

  call regu(frc_regu, f_vector)
  call regv(frc_regv, f_vector)
  call regw(frc_regw, f_vector)

  do i = 1, nx
    do j = 1, nyg
      do k = 1, nzg
        if (wghtu(i,j,k) > 1.d-10) then
          frc_regu(i,j,k) = frc_regu(i,j,k)/wghtu(i,j,k)
        end if
      end do
    end do
  end do

  do i = 1, nxg
    do j = 1, ny
      do k = 1, nzg
        if (wghtv(i,j,k) > 1.d-10) then
          frc_regv(i,j,k) = frc_regv(i,j,k)/wghtv(i,j,k)
        end if
      end do
    end do
  end do

  do i = 1, nxg
    do j = 1, nyg
      do k = 1, nz
        if (wghtw(i,j,k) > 1.d-10) then
          frc_regw(i,j,k) = frc_regw(i,j,k)/wghtw(i,j,k)
        end if
      end do
    end do
  end do

  redistribute = regT(frc_regu, frc_regv, frc_regw)

  deallocate(wghtu, frc_regu)
  deallocate(wghtv, frc_regv)
  deallocate(wghtw, frc_regw)

end function redistribute


function ItildeGprime_subsurface(x)

  REAL(KIND(0.d0)), DIMENSION(nblocks), INTENT(in) :: x
  real(kind(0.d0)), dimension(3*nb) :: ItildeGprime_subsurface
  real(kind(0.d0)), dimension(3*nb) :: v_patch_3d
  real(kind(0.d0)), dimension(nb) :: v_patch
  integer :: i

  do i = 1, nb
    v_patch(i) = cos((3.d0*pi/Lxp)*(xb(i)-xtopmass)) * &
                 exp(-(xb(i)-xtopmass)**2.d0/(2.d0*sigma**2.d0)) * x(1)
  end do

  v_patch_3d = 0.d0

  do i = 1, nb
    v_patch_3d(nb+i) = v_patch(i)
  end do

  ItildeGprime_subsurface = v_patch_3d

end function ItildeGprime_subsurface


subroutine check_slip_FSI

  Real(Int64) :: max_slip
  real(kind(0.d0)), dimension(3*nb) :: slip_fsi

  slip_fsi = regT(U, V, W)


  If ( trim(body_type) ==  'center_wall_deforming_testcase') then
  slip_fsi = slip_fsi - Itilde(zeta)
  end if
  
   If ( trim(body_type) ==  'center_wall_deforming_subsurface') then
  slip_fsi = slip_fsi - ItildeGprime_subsurface(zeta)
  end if 

  max_slip = Maxval(Abs(slip_fsi))
  Write(*,*) 'Maximum slip for fsi         : ', max_slip

end subroutine check_slip_FSI


subroutine writechimax(disp)

  real(kind(0.d0)), dimension(nb), intent(in) :: disp

  open(unit=10001, file="chi_top.dat", form="formatted", &
       status="unknown", position="append")

  write(10001,*) maxval(disp)

  close(10001)

end subroutine writechimax

subroutine writechimax_subsurface(disp)

  real(kind(0.d0)), dimension(nblocks), intent(in) :: disp

  open(unit=10001, file="chi_top.dat", form="formatted", &
       status="unknown", position="append")

  write(10001,*) (disp(1))

  close(10001)

end subroutine writechimax_subsurface


subroutine writestuffchi(disp)

  implicit none
  real(kind(0.d0)), dimension(nb), intent(in) :: disp
  integer, save :: cid = 0
  integer :: k
  character(len=50) :: fname

  cid = cid + 1
  write(fname,'("chi_",I0,".dat")') cid

  open(unit=10001, file=fname, form="formatted", status="replace")

  do k = 1, nb
    write(10001,*) disp(k)
  end do

  close(10001)

end subroutine writestuffchi


subroutine writestuffiterfsi(disp)

  integer, intent(in) :: disp

  open(unit=10001, file="iterfsi.dat", form="formatted", &
       status="unknown", position="append")

  write(10001,*) disp

  close(10001)

end subroutine writestuffiterfsi


subroutine writestuffsurfacestress(disp)

  implicit none
  real(kind(0.d0)), dimension(3*nb), intent(in) :: disp
  integer, save :: file_id = 0
  integer :: k
  character(len=50) :: fname

  file_id = file_id + 1
  write(fname,'("stress_",I0,".dat")') file_id

  open(unit=10001, file=fname, form="formatted", status="replace")

  do k = 1, 3*nb
    write(10001,*) disp(k)
  end do

  close(10001)

end subroutine writestuffsurfacestress


subroutine writestuffyb(yb_in)

  implicit none
  real(kind(0.d0)), dimension(nb), intent(in) :: yb_in
  integer :: k
  integer, save :: yid = 0
  character(len=50) :: fname

  yid = yid + 1
  write(fname,'("yb_",I0,".dat")') yid

  open(unit=10001, file=fname, form="formatted", status="replace")

  do k = 1, nb
    write(10001,*) yb_in(k)
  end do

  close(10001)

end subroutine writestuffyb

End Module operators_FSI
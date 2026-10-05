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

function KhatinvQItildeprimeW_subsurface(x)

  REAL(KIND(0.d0)), DIMENSION(3*nb), INTENT(in) :: x
  real(kind(0.d0)), dimension(nblocks) :: KhatinvQItildeprimeW_subsurface
  real(kind(0.d0)), dimension(3*nb) :: v_tp
  real(kind(0.d0)), dimension(nb) :: v_y
  real(kind(0.d0)), dimension(nblocks) :: v_trunc, v_bg
  integer :: i
  real(kind(0.d0)) :: sum

  v_tp = redistribute(dt_fsi*x)
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
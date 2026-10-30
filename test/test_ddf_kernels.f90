program test_ddf_kernels
  use, intrinsic :: iso_fortran_env, only : int32, int64
  use global, only : use_bao5_ddf, suppx, suppy, suppz, nweights
  use immersed_boundary_operators, only : deltafnc, ddf_support_radius
  implicit none

  real(int64), parameter :: moment = (38d0 - sqrt(69d0)) / 60d0
  real(int64), parameter :: square_sum = ((9d0 - 4d0 * moment)**2 + (4d0 * moment - 1d0)**2) / 128d0
  real(int64), parameter :: spacings(3) = [1d0, 0.125d0, 9.606358701352846d-4]
  real(int64) :: h, marker, grid, r, weight, moments(0:3), squares
  integer(int32) :: spacing, phase, index, pivot

  call require(.not. use_bao5_ddf, 'Yang3 must remain the default')
  call require(suppx == 2 .and. suppy == 2 .and. suppz == 2 .and. nweights == 125, &
    'kernel selection must retain the existing stencil size')
  call require(ddf_support_radius() == 2d0, 'Yang3 support radius')
  call require(deltafnc(2.25d0, 0d0, 1d0) == 0d0, 'Yang3 must vanish outside two cells')

  use_bao5_ddf = .true.
  call require(ddf_support_radius() == 2.5d0, 'Bao5 support radius')
  call require(deltafnc(2.25d0, 0d0, 1d0) > 0d0, 'flag must activate the wider kernel')
  do spacing = 1, size(spacings)
    h = spacings(spacing)
    do phase = 0, 2000
      marker = (real(phase, int64) / 2000d0 - 0.5d0) * h
      pivot = nint(marker / h)
      moments = 0d0
      squares = 0d0
      do index = pivot - 2, pivot + 2
        grid = real(index, int64) * h
        r = (grid - marker) / h
        weight = h * deltafnc(grid, marker, h)
        call require(weight >= -1d-14, 'Bao5 must be nonnegative')
        call require(abs(deltafnc(grid, marker, h) - deltafnc(marker, grid, h)) * h < 1d-14, &
          'Bao5 must be symmetric')
        moments = moments + weight * [1d0, r, r*r, r*r*r]
        squares = squares + weight*weight
      end do
      call require(abs(moments(0) - 1d0) < 2d-13, 'partition of unity at every grid phase')
      call require(abs(moments(1)) < 2d-13, 'zero first moment at every grid phase')
      call require(abs(moments(2) - moment) < 2d-13, 'constant second moment at every grid phase')
      call require(abs(moments(3)) < 2d-13, 'zero third moment at every grid phase')
      call require(abs(squares - square_sum) < 2d-13, 'constant self coupling at every grid phase')
    end do
    call require(deltafnc(2.5d0*h, 0d0, h) == 0d0, 'Bao5 support endpoint')
    call require(deltafnc(3d0*h, 0d0, h) == 0d0, 'Bao5 compact support')
  end do

  use_bao5_ddf = .false.
  call require(ddf_support_radius() == 2d0, 'flag must restore Yang3 support')
  call require(deltafnc(2.25d0, 0d0, 1d0) == 0d0, 'flag must restore Yang3 kernel')

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(A)') 'FAIL: ' // message
      error stop 'DDF kernel tests failed'
    end if
  end subroutine require

end program test_ddf_kernels

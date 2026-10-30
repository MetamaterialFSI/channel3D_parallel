program test_bicgstab_metrics
  use, intrinsic :: iso_fortran_env, only : int32, int64
  use projection, only : compute_surface_block_rms, compute_hybrid_block_tolerances, &
    compute_scaled_block_residual
  implicit none

  integer(int32) :: failures
  real(int64), parameter :: rtol = 1d-2, atol = 1d-4
  real(int64) :: weights2(2), weights4(4), rhs(6), residual(6), rotated(6)
  real(int64) :: scaled_rhs(6), scaled_residual_vector(6)
  real(int64) :: vector2(6), vector4(12)
  real(int64) :: component_rms(3), grouped_rms(2), rhs_rms(3), tolerances(2), eta(2)
  real(int64) :: component_rms4(3), grouped_rms4(2), rotated_rms(3), rotated_eta(2)
  real(int64) :: scaled_rhs_rms(3), scaled_tolerances(2), scaled_component_rms(3), scaled_eta(2)

  failures = 0
  weights2 = [1d0, 2d0]
  rhs = 0d0
  call compute_hybrid_block_tolerances(rhs, weights2, rtol, atol, rhs_rms, tolerances)
  call require_close(tolerances(1), atol, 'zero-RHS tangent tolerance', failures)
  call require_close(tolerances(2), atol, 'zero-RHS pressure tolerance', failures)

  residual = 0d0
  residual(5:6) = 0.5d0 * atol
  call compute_scaled_block_residual(residual, weights2, tolerances, component_rms, eta)
  call require_close(eta(2), 0.5d0, 'absolute-floor passing residual', failures)
  residual(5:6) = 2d0 * atol
  call compute_scaled_block_residual(residual, weights2, tolerances, component_rms, eta)
  call require(eta(2) > 1d0 .and. eta(1) == 0d0, 'one failing block must reject', failures)

  rhs = [3d0, 3d0, 4d0, 4d0, 5d0, 5d0]
  call compute_hybrid_block_tolerances(rhs, weights2, rtol, atol, rhs_rms, tolerances)
  call require_close(rhs_rms(1), 3d0, 'relative branch first tangent RHS RMS', failures)
  call require_close(rhs_rms(2), 4d0, 'relative branch second tangent RHS RMS', failures)
  call require_close(rhs_rms(3), 5d0, 'relative branch pressure RHS RMS', failures)
  call require_close(tolerances(1), rtol * sqrt(12.5d0), &
    'relative branch tangent tolerance', failures)
  call require_close(tolerances(2), 5d0 * rtol, 'relative branch pressure tolerance', failures)
  call require(all(tolerances > atol), 'relative tolerance must dominate both blocks', failures)

  residual = 0.5d0 * rtol * rhs
  call compute_scaled_block_residual(residual, weights2, tolerances, component_rms, eta)
  call require_close(eta(1), 0.5d0, 'relative tangent residual ratio', failures)
  call require_close(eta(2), 0.5d0, 'relative pressure residual ratio', failures)

  scaled_rhs = rhs
  scaled_residual_vector = residual
  scaled_rhs(1:4) = 100d0 * scaled_rhs(1:4)
  scaled_residual_vector(1:4) = 100d0 * scaled_residual_vector(1:4)
  scaled_rhs(5:6) = 0.1d0 * scaled_rhs(5:6)
  scaled_residual_vector(5:6) = 0.1d0 * scaled_residual_vector(5:6)
  call compute_hybrid_block_tolerances(scaled_rhs, weights2, rtol, atol, &
    scaled_rhs_rms, scaled_tolerances)
  call require(all(scaled_tolerances > atol), &
    'rescaled blocks must remain relative-tolerance dominated', failures)
  call compute_scaled_block_residual(scaled_residual_vector, weights2, scaled_tolerances, &
    scaled_component_rms, scaled_eta)
  call require(maxval(abs(eta - scaled_eta)) < 1d-14, &
    'relative residual ratios must be independently block-rescaling invariant', failures)

  residual = 0d0
  residual(1:2) = sqrt(2d0) * atol
  call compute_scaled_block_residual(residual, weights2, tolerances, component_rms, eta)
  rotated = 0d0
  rotated(1:4) = atol
  call compute_scaled_block_residual(rotated, weights2, tolerances, rotated_rms, rotated_eta)
  call require_close(eta(1), rotated_eta(1), 'tangent-basis rotation invariance', failures)

  vector2 = [1d0, 3d0, 2d0, 4d0, 5d0, 7d0]
  weights4 = [0.5d0, 0.5d0, 1d0, 1d0]
  vector4 = [1d0, 1d0, 3d0, 3d0, 2d0, 2d0, 4d0, 4d0, 5d0, 5d0, 7d0, 7d0]
  call compute_surface_block_rms(vector2, weights2, component_rms, grouped_rms)
  call compute_surface_block_rms(vector4, weights4, component_rms4, grouped_rms4)
  call require(maxval(abs(component_rms - component_rms4)) < 1d-14, &
    'surface RMS must be marker-count invariant', failures)
  call require(maxval(abs(grouped_rms - grouped_rms4)) < 1d-14, &
    'grouped RMS must be marker-count invariant', failures)

  if (failures /= 0) error stop 'BiCGSTAB metric tests failed'

contains

  subroutine require(condition, message, failure_count)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer(int32), intent(inout) :: failure_count

    if (.not. condition) then
      write(*, '(A)') 'FAIL: ' // message
      failure_count = failure_count + 1
    end if
  end subroutine require

  subroutine require_close(actual, expected, message, failure_count)
    real(int64), intent(in) :: actual, expected
    character(len=*), intent(in) :: message
    integer(int32), intent(inout) :: failure_count

    call require(abs(actual - expected) <= 1d-13 * max(1d0, abs(expected)), message, failure_count)
  end subroutine require_close

end program test_bicgstab_metrics

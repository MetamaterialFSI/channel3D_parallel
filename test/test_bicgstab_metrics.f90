program test_bicgstab_metrics
  use, intrinsic :: iso_fortran_env, only : int32, int64
  use projection, only : compute_surface_block_rms, compute_hybrid_block_tolerances, &
    compute_scaled_block_residual
  implicit none

  integer(int32) :: failures
  real(int64), parameter :: rtol = 1d-2, atol = 1d-4
  real(int64), parameter :: inv_sqrt_two = 1d0 / sqrt(2d0)
  real(int64), parameter :: inv_sqrt_three = 1d0 / sqrt(3d0)
  real(int64), parameter :: inv_sqrt_six = 1d0 / sqrt(6d0)
  real(int64) :: weights2(2), weights4(4), rhs(6), residual(6)
  real(int64) :: scaled_rhs(6), scaled_residual_vector(6), vector2(6), vector4(12)
  real(int64) :: normals2(6), tangent1_2(6), tangent2_2(6)
  real(int64) :: normals4(12), tangent1_4(12), tangent2_4(12)
  real(int64) :: rotated_tangent1(6), rotated_tangent2(6)
  real(int64) :: oblique_vector(6), oblique_normals(6), oblique_tangent1(6), oblique_tangent2(6)
  real(int64) :: component_rms(3), grouped_rms(2), rhs_rms(3), tolerances(2), eta(2)
  real(int64) :: component_rms4(3), grouped_rms4(2), rotated_rms(3), rotated_eta(2)
  real(int64) :: scaled_rhs_rms(3), scaled_tolerances(2), scaled_component_rms(3), scaled_eta(2)

  failures = 0
  weights2 = [1d0, 2d0]
  call set_cartesian_frames(normals2, tangent1_2, tangent2_2)

  rhs = 0d0
  call compute_hybrid_block_tolerances(rhs, weights2, normals2, tangent1_2, tangent2_2, &
    rtol, atol, rhs_rms, tolerances)
  call require_close(tolerances(1), atol, 'zero-RHS tangent tolerance', failures)
  call require_close(tolerances(2), atol, 'zero-RHS normal tolerance', failures)

  residual = 0d0
  residual(5:6) = 0.5d0 * atol
  call compute_scaled_block_residual(residual, weights2, normals2, tangent1_2, tangent2_2, &
    tolerances, component_rms, eta)
  call require_close(eta(2), 0.5d0, 'absolute-floor passing residual', failures)
  residual(5:6) = 2d0 * atol
  call compute_scaled_block_residual(residual, weights2, normals2, tangent1_2, tangent2_2, &
    tolerances, component_rms, eta)
  call require(eta(2) > 1d0 .and. eta(1) == 0d0, 'one failing block must reject', failures)

  rhs = [3d0, 3d0, 4d0, 4d0, 5d0, 5d0]
  call compute_hybrid_block_tolerances(rhs, weights2, normals2, tangent1_2, tangent2_2, &
    rtol, atol, rhs_rms, tolerances)
  call require_close(rhs_rms(1), 3d0, 'relative branch first tangent RHS RMS', failures)
  call require_close(rhs_rms(2), 4d0, 'relative branch second tangent RHS RMS', failures)
  call require_close(rhs_rms(3), 5d0, 'relative branch normal RHS RMS', failures)
  call require_close(tolerances(1), rtol * sqrt(12.5d0), &
    'relative branch tangent tolerance', failures)
  call require_close(tolerances(2), 5d0 * rtol, 'relative branch normal tolerance', failures)
  call require(all(tolerances > atol), 'relative tolerance must dominate both blocks', failures)

  call set_oblique_frames(oblique_normals, oblique_tangent1, oblique_tangent2)
  oblique_vector(1:2) = 3d0 * inv_sqrt_two + 4d0 * inv_sqrt_six + &
    5d0 * inv_sqrt_three
  oblique_vector(3:4) = -3d0 * inv_sqrt_two + 4d0 * inv_sqrt_six + &
    5d0 * inv_sqrt_three
  oblique_vector(5:6) = -8d0 * inv_sqrt_six + 5d0 * inv_sqrt_three
  call compute_surface_block_rms(oblique_vector, weights2, oblique_normals, oblique_tangent1, &
    oblique_tangent2, component_rms, grouped_rms)
  call require_close(component_rms(1), 3d0, 'oblique-frame first tangent projection', failures)
  call require_close(component_rms(2), 4d0, 'oblique-frame second tangent projection', failures)
  call require_close(component_rms(3), 5d0, 'oblique-frame normal projection', failures)

  residual = 0.5d0 * rtol * rhs
  call compute_scaled_block_residual(residual, weights2, normals2, tangent1_2, tangent2_2, &
    tolerances, component_rms, eta)
  call require_close(eta(1), 0.5d0, 'relative tangent residual ratio', failures)
  call require_close(eta(2), 0.5d0, 'relative normal residual ratio', failures)

  scaled_rhs = rhs
  scaled_residual_vector = residual
  scaled_rhs(1:4) = 100d0 * scaled_rhs(1:4)
  scaled_residual_vector(1:4) = 100d0 * scaled_residual_vector(1:4)
  scaled_rhs(5:6) = 0.1d0 * scaled_rhs(5:6)
  scaled_residual_vector(5:6) = 0.1d0 * scaled_residual_vector(5:6)
  call compute_hybrid_block_tolerances(scaled_rhs, weights2, normals2, tangent1_2, tangent2_2, &
    rtol, atol, scaled_rhs_rms, scaled_tolerances)
  call require(all(scaled_tolerances > atol), &
    'rescaled blocks must remain relative-tolerance dominated', failures)
  call compute_scaled_block_residual(scaled_residual_vector, weights2, normals2, tangent1_2, &
    tangent2_2, scaled_tolerances, scaled_component_rms, scaled_eta)
  call require(maxval(abs(eta - scaled_eta)) < 1d-14, &
    'relative residual ratios must be independently block-rescaling invariant', failures)

  residual = 0d0
  residual(1:2) = sqrt(2d0) * atol
  call compute_scaled_block_residual(residual, weights2, normals2, tangent1_2, tangent2_2, &
    tolerances, component_rms, eta)
  rotated_tangent1 = 0d0
  rotated_tangent2 = 0d0
  rotated_tangent1(1:2) = inv_sqrt_two
  rotated_tangent1(3:4) = inv_sqrt_two
  rotated_tangent2(1:2) = -inv_sqrt_two
  rotated_tangent2(3:4) = inv_sqrt_two
  call compute_scaled_block_residual(residual, weights2, normals2, rotated_tangent1, &
    rotated_tangent2, tolerances, rotated_rms, rotated_eta)
  call require_close(eta(1), rotated_eta(1), 'tangent-basis rotation invariance', failures)

  vector2 = [1d0, 3d0, 2d0, 4d0, 5d0, 7d0]
  weights4 = [0.5d0, 0.5d0, 1d0, 1d0]
  vector4 = [1d0, 1d0, 3d0, 3d0, 2d0, 2d0, 4d0, 4d0, 5d0, 5d0, 7d0, 7d0]
  call set_cartesian_frames(normals4, tangent1_4, tangent2_4)
  call compute_surface_block_rms(vector2, weights2, normals2, tangent1_2, tangent2_2, &
    component_rms, grouped_rms)
  call compute_surface_block_rms(vector4, weights4, normals4, tangent1_4, tangent2_4, &
    component_rms4, grouped_rms4)
  call require(maxval(abs(component_rms - component_rms4)) < 1d-14, &
    'surface RMS must be marker-count invariant', failures)
  call require(maxval(abs(grouped_rms - grouped_rms4)) < 1d-14, &
    'grouped RMS must be marker-count invariant', failures)

  if (failures /= 0) error stop 'BiCGSTAB metric tests failed'

contains

  subroutine set_cartesian_frames(normals, tangent1, tangent2)
    real(int64), intent(out) :: normals(:), tangent1(:), tangent2(:)
    integer(int32) :: marker_count

    marker_count = size(normals) / 3
    normals = 0d0
    tangent1 = 0d0
    tangent2 = 0d0
    tangent1(1:marker_count) = 1d0
    tangent2(marker_count + 1:2 * marker_count) = 1d0
    normals(2 * marker_count + 1:3 * marker_count) = 1d0
  end subroutine set_cartesian_frames

  subroutine set_oblique_frames(normals, tangent1, tangent2)
    real(int64), intent(out) :: normals(:), tangent1(:), tangent2(:)
    integer(int32) :: marker_count

    marker_count = size(normals) / 3
    tangent1(1:marker_count) = inv_sqrt_two
    tangent1(marker_count + 1:2 * marker_count) = -inv_sqrt_two
    tangent1(2 * marker_count + 1:3 * marker_count) = 0d0
    tangent2(1:marker_count) = inv_sqrt_six
    tangent2(marker_count + 1:2 * marker_count) = inv_sqrt_six
    tangent2(2 * marker_count + 1:3 * marker_count) = -2d0 * inv_sqrt_six
    normals = inv_sqrt_three
  end subroutine set_oblique_frames

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

program channel3d_verification
  use, intrinsic :: iso_fortran_env, only : int32, int64
  use global
  use input_output, only : read_input_parameters
  use initialization, only : initialize, initialize_ib_arrays
  use immersed_boundary_geometry, only : setup_IB_geometry
  use immersed_boundary_operators, only : setup_IB_operators
  use heaviside, only : compute_heaviside
  use boundary_conditions, only : apply_boundary_conditions
  use mass_flow, only : compute_mean_mass_flow_U, compute_mean_mass_flow_V, compute_mean_mass_flow_W
  use time_integration, only : compute_dt, compute_time_step_Euler, compute_time_step_RK3
  use finalization, only : finalize
  use mpi
  implicit none

  character(len=200) :: verification_case
  character(len=200) :: bulk_clearance_argument
  real(int64) :: h, rel_l1, rel_l2, rel_linf, rel_linf_bulk, rel_linf_support
  real(int64) :: bulk_clearance, inactive_linf
  real(int64) :: requested_bulk_clearance

  call MPI_init(ierr)
  call MPI_comm_size(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_comm_rank(MPI_COMM_WORLD, myid, ierr)

  call get_command_argument(1, fileparams)
  call get_command_argument(2, verification_case)
  call get_command_argument(3, bulk_clearance_argument)
  if (len_trim(fileparams) == 0 .or. len_trim(verification_case) == 0) then
    if (myid == 0) write(*, '(A)') 'usage: channel3d_verification INPUT CASE'
    call MPI_abort(MPI_COMM_WORLD, 2, ierr)
  end if
  requested_bulk_clearance = -1d0
  if (len_trim(bulk_clearance_argument) > 0) read(bulk_clearance_argument, *) requested_bulk_clearance

  call read_input_parameters
  call initialize
  call initialize_ib_arrays
  call setup_IB_geometry
  call setup_IB_operators
  call compute_heaviside
  call set_exact_velocity(trim(verification_case))
  call apply_boundary_conditions(U, V, W)

  P = 0d0
  call compute_mean_mass_flow_U(U, Qflow_x_0)
  call compute_mean_mass_flow_V(V, Qflow_y_0)
  call compute_mean_mass_flow_W(W, Qflow_z_0)
  Qflow_y_0 = 0d0
  dPdy = 0d0

  do istep = 1, nsteps
    call compute_dt
    if (use_rk3) then
      call compute_time_step_RK3
    else
      call compute_time_step_Euler
    end if
  end do

  call compute_velocity_errors(trim(verification_case), h, rel_l1, rel_l2, rel_linf, rel_linf_bulk, &
    rel_linf_support, bulk_clearance, inactive_linf)
  if (myid == 0) then
    write(*, '(A,1X,A,9(1X,ES24.16))') 'CHANNEL3D_METRICS', trim(verification_case), &
      h, rel_l1, rel_l2, rel_linf, rel_linf_bulk, rel_linf_support, bulk_clearance, inactive_linf, t
  end if

  call finalize

contains

  real(int64) function center_wall_profile(y_coord, forcing)
    real(int64), intent(in) :: y_coord, forcing
    real(int64) :: wall

    wall = 0.5d0 * Ly_channel
    if (y_coord < wall) then
      center_wall_profile = forcing / (2d0 * nu) * y_coord * (wall - y_coord)
    else
      center_wall_profile = forcing / (2d0 * nu) * (y_coord - wall) * (Ly_channel - y_coord)
    end if
  end function center_wall_profile

  subroutine cylinder_velocity(x_coord, y_coord, u_exact, v_exact)
    real(int64), intent(in) :: x_coord, y_coord
    real(int64), intent(out) :: u_exact, v_exact
    real(int64) :: r1, r2, xc, yc, x_rel, y_rel, radius, omega, azimuthal_velocity

    r1 = body_param_1
    r2 = body_param_2
    xc = 0.5d0 * Lxp
    yc = 0.5d0 * Ly_channel
    x_rel = x_coord - xc
    y_rel = y_coord - yc
    radius = sqrt(x_rel * x_rel + y_rel * y_rel)
    omega = body_param_3 / r1

    if (radius < r1) then
      azimuthal_velocity = -omega * radius
    else if (radius <= r2) then
      azimuthal_velocity = -omega * r1 * r1 / (r2 * r2 - r1 * r1) * &
        (r2 * r2 / radius - radius)
    else
      azimuthal_velocity = 0d0
    end if

    if (radius > 0d0) then
      u_exact = azimuthal_velocity * y_rel / radius
      v_exact = -azimuthal_velocity * x_rel / radius
    else
      u_exact = 0d0
      v_exact = 0d0
    end if
  end subroutine cylinder_velocity

  subroutine set_exact_velocity(case_name)
    character(len=*), intent(in) :: case_name
    integer(int32) :: i, j, k
    real(int64) :: u_exact, v_exact

    U = 0d0
    V = 0d0
    W = 0d0

    select case (case_name)
      case ('center_wall_x')
        if (trim(body_type) /= 'center_wall') error stop 'center_wall_x requires body_type=center_wall'
        do k = 1, nzg
          do j = 2, nyg - 1
            U(:, j, k) = center_wall_profile(yg(j), dPdx)
          end do
        end do

      case ('center_wall_z')
        if (trim(body_type) /= 'center_wall') error stop 'center_wall_z requires body_type=center_wall'
        do k = 1, nz
          do j = 2, nyg - 1
            W(:, j, k) = center_wall_profile(yg(j), dPdz)
          end do
        end do

      case ('double_cylinder')
        if (trim(body_type) /= 'double_cylinders_z') &
          error stop 'double_cylinder requires body_type=double_cylinders_z'
        do k = 1, nzg
          do j = 2, nyg - 1
            do i = 2, nx - 1
              call cylinder_velocity(x(i), yg(j), u_exact, v_exact)
              U(i, j, k) = u_exact
            end do
          end do
        end do
        do k = 1, nzg
          do j = 2, ny - 1
            do i = 2, nxg - 2
              call cylinder_velocity(xg(i), y(j), u_exact, v_exact)
              V(i, j, k) = v_exact
            end do
          end do
        end do

      case default
        error stop 'unknown verification case'
    end select
  end subroutine set_exact_velocity

  subroutine accumulate_error(value, exact_value, in_bulk, outside_support, sums, maxima)
    real(int64), intent(in) :: value, exact_value
    logical, intent(in) :: in_bulk, outside_support
    real(int64), intent(inout) :: sums(4), maxima(4)
    real(int64) :: error_value

    error_value = abs(value - exact_value)
    sums(1) = sums(1) + error_value
    sums(2) = sums(2) + error_value * error_value
    sums(3) = sums(3) + 1d0
    maxima(1) = max(maxima(1), error_value)
    if (outside_support) maxima(2) = max(maxima(2), error_value)
    if (in_bulk) maxima(3) = max(maxima(3), error_value)
  end subroutine accumulate_error

  subroutine compute_velocity_errors(case_name, grid_h, l1_error, l2_error, linf_error, linf_bulk_error, &
    linf_support_error, bulk_clearance, inactive_error)
    character(len=*), intent(in) :: case_name
    real(int64), intent(out) :: grid_h, l1_error, l2_error, linf_error, linf_bulk_error, linf_support_error
    real(int64), intent(out) :: bulk_clearance, inactive_error
    integer(int32) :: i, j, k, center_k_end
    real(int64) :: exact_value, u_exact, v_exact, velocity_scale, native_support_clearance
    real(int64) :: local_sums(4), global_sums(4), local_maxima(4), global_maxima(4)
    logical :: in_bulk, outside_support

    local_sums = 0d0
    local_maxima = 0d0
    if (myid < nprocs - 1) then
      center_k_end = nzg - 1
    else
      center_k_end = nzg - 2
    end if

    select case (case_name)
      case ('center_wall_x')
        native_support_clearance = real(suppy + 1, int64) * dymin
        bulk_clearance = native_support_clearance
        if (requested_bulk_clearance > 0d0) bulk_clearance = requested_bulk_clearance
        do k = 2, center_k_end
          do j = 2, nyg - 1
            exact_value = center_wall_profile(yg(j), dPdx)
            in_bulk = abs(yg(j) - 0.5d0 * Ly_channel) > bulk_clearance
            outside_support = abs(yg(j) - 0.5d0 * Ly_channel) > native_support_clearance
            do i = 2, nx - 1
              call accumulate_error(U(i, j, k), exact_value, in_bulk, outside_support, local_sums, local_maxima)
            end do
          end do
        end do
        local_maxima(4) = max(maxval(abs(V)), maxval(abs(W)))
        grid_h = dymin
        velocity_scale = abs(dPdx) * Ly_channel * Ly_channel / (32d0 * nu)

      case ('center_wall_z')
        native_support_clearance = real(suppy + 1, int64) * dymin
        bulk_clearance = native_support_clearance
        if (requested_bulk_clearance > 0d0) bulk_clearance = requested_bulk_clearance
        do k = 2, nz - 1
          do j = 2, nyg - 1
            exact_value = center_wall_profile(yg(j), dPdz)
            in_bulk = abs(yg(j) - 0.5d0 * Ly_channel) > bulk_clearance
            outside_support = abs(yg(j) - 0.5d0 * Ly_channel) > native_support_clearance
            do i = 2, nxg - 2
              call accumulate_error(W(i, j, k), exact_value, in_bulk, outside_support, local_sums, local_maxima)
            end do
          end do
        end do
        local_maxima(4) = max(maxval(abs(U)), maxval(abs(V)))
        grid_h = dymin
        velocity_scale = abs(dPdz) * Ly_channel * Ly_channel / (32d0 * nu)

      case ('double_cylinder')
        bulk_clearance = 0d0
        do k = 2, center_k_end
          do j = 2, nyg - 1
            do i = 2, nx - 1
              call cylinder_velocity(x(i), yg(j), u_exact, v_exact)
              call accumulate_error(U(i, j, k), u_exact, .true., .true., local_sums, local_maxima)
            end do
          end do
        end do
        do k = 2, center_k_end
          do j = 2, ny - 1
            do i = 2, nxg - 2
              call cylinder_velocity(xg(i), y(j), u_exact, v_exact)
              call accumulate_error(V(i, j, k), v_exact, .true., .true., local_sums, local_maxima)
            end do
          end do
        end do
        local_maxima(4) = maxval(abs(W))
        grid_h = max(dxmin, dymin)
        velocity_scale = abs(body_param_3)

      case default
        error stop 'unknown verification case'
    end select

    call MPI_reduce(local_sums, global_sums, 4, MPI_real8, MPI_sum, 0, MPI_COMM_WORLD, ierr)
    call MPI_reduce(local_maxima, global_maxima, 4, MPI_real8, MPI_max, 0, MPI_COMM_WORLD, ierr)

    if (myid == 0) then
      l1_error = global_sums(1) / global_sums(3) / velocity_scale
      l2_error = sqrt(global_sums(2) / global_sums(3)) / velocity_scale
      linf_error = global_maxima(1) / velocity_scale
      linf_bulk_error = global_maxima(3) / velocity_scale
      linf_support_error = global_maxima(2) / velocity_scale
      inactive_error = global_maxima(4)
    else
      l1_error = 0d0
      l2_error = 0d0
      linf_error = 0d0
      linf_bulk_error = 0d0
      linf_support_error = 0d0
      inactive_error = 0d0
    end if
  end subroutine compute_velocity_errors

end program channel3d_verification

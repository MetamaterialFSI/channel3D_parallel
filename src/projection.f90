Module projection

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use, Intrinsic :: ieee_arithmetic, Only : ieee_is_finite
  Use global
  Use immersed_boundary_operators
  Use immersed_boundary_geometry
  Use equations
  Use poisson
  Use boundary_conditions
  Use mass_flow
  Use differential_operators

  ! prevent implicit typing
  Implicit None

Contains

  !--------------------------------------------!
  !    Compute incompressible velocity with    !
  !         fractional step method             !
  !--------------------------------------------!
  Subroutine compute_non_IB_projection

    ! (U,V,W)_interim store the intermediate velocity fields from solving the NS terms without pressure gradient or forcing
    ! save for use in the subsequent projection steps
    U_interim = U
    V_interim = V
    W_interim = W
  
    Call divergence(rhs_p, U, V, W)

    Call solve_poisson_equation(rhs_p)
    rhs_p = rhs_p / dt

    ! save for use in the subsequent projection steps
    P_interim = rhs_p

    rhs_p = dt * rhs_p
    Call gradient(U, V, W, rhs_p)

    ! U* = U** - Gp*
    U = U_interim - U
    V = V_interim - V
    W = W_interim - W

  End Subroutine compute_non_IB_projection

  Subroutine compute_IB_projection

    ! - E u* + ub 
    Call regT(Eu, U, V, W)
    ib_surface_vector_work = ub - Eu
    Call cart_to_frame_t(ib_surface_vector_work, rhs_ib(1 : 2 * nb))
    rhs_ib(2 * nb + 1 : 3 * nb) = -(ib_surface_vector_work(1 : nb) * normals(1 : nb) &
      + ib_surface_vector_work(nb + 1 : 2 * nb) * normals(nb + 1 : 2 * nb) &
      + ib_surface_vector_work(2 * nb + 1 : 3 * nb) * normals(2 * nb + 1 : 3 * nb))
    aux_surface_scalar = 0d0
    Call dpmean_ds(P_interim, 1d0, aux_surface_scalar, dpds_known, ring_tang)
    Call regTc_1n(E1np, P_interim)
    rhs_ib(2 * nb + 1 : 3 * nb) = rhs_ib(2 * nb + 1 : 3 * nb) &
      - E1np + E2n1c * pn_mean + Ecross * dpds_known

    ! Remove mean from rhs pressure

    ! solve for IB forcing
    Call bicgstab(fb, rhs_ib)

    Call frame_to_cart_t(fb, dudn_jump)
    p_jump = fb(2 * nb + 1 : 3 * nb)
    Call remove_mean_per_body(p_jump)

    ! U_reg = R f
    ib_surface_vector_work = kdt * dudn_jump
    aux_surface_scalar = kdt * p_jump
    Call apply_momentum_forcing(U_reg, V_reg, W_reg, ib_surface_vector_work, aux_surface_scalar)

    ! rhs_p = D R f
    Call divergence(rhs_p, U_reg, V_reg, W_reg)
    Call apply_continuity_forcing(con_F, dudn_jump)
    rhs_p = rhs_p + con_F

    ! rhs_p = Linv D R f
    Call solve_poisson_equation(rhs_p, .False.)
    rhs_p = rhs_p / dt

    ! Pnp1 = P* - Linv D R f
    rhs_p = P_interim - rhs_p
    P = 0d0
    P(2:nxg-1, 2:nyg-1, 2:nzg-1) = rhs_p(2:nxg-1, 2:nyg-1, 2:nzg-1)
    Call apply_periodic_bc_x(P, 2)
    P(:, 1, :) = P(:, 2, :)
    P(:, nyg, :) = P(:, nyg-1, :)
    If (nprocs == 1) Then
      P(:, :, 1) = P(:, :, nzg-2)
      P(:, :, nzg-1) = P(:, :, 2)
      P(:, :, nzg) = P(:, :, 3)
    Else
      Call update_ghost_interior_planes(P, 4)
      Call apply_periodic_bc_z(P, 4)
    End If
    pressure_computed = .True.

    ! U, V, W = G Pnp1
    rhs_p = dt * rhs_p
    Call gradient(U, V, W, rhs_p)

    ! Unp1 = U** - R f - G Pnp1
    U = U_interim - U_reg - U
    V = V_interim - V_reg - V
    W = W_interim - W_reg - W

    Call apply_boundary_conditions(U, V, W)

  End Subroutine compute_IB_projection

  Subroutine schur(Sf_, f_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: Sf_(:)
    Real(Int64), Dimension(nbodies) :: pavg

    Sf_ = 0.d0

    Call frame_to_cart_t(f_, dudn_jump)
    p_jump = f_(2 * nb + 1 : 3 * nb)
    Call remove_mean_per_body(p_jump)

    ! U_reg = R f
    ib_surface_vector_work = kdt * dudn_jump
    aux_surface_scalar = kdt * p_jump
    Call apply_momentum_forcing(U_reg, V_reg, W_reg, ib_surface_vector_work, aux_surface_scalar)
    Call apply_continuity_forcing(con_F, dudn_jump)

    ! rhs_p = D R f
    Call divergence(rhs_p, U_reg, V_reg, W_reg)
    rhs_p = rhs_p + con_F

    ! rhs_p = Linv D R f
    call solve_poisson_equation(rhs_p, .False.)
    rhs_p = rhs_p / dt

    Call regTc_1n(E1np, rhs_p)
    Call dpmean_ds(rhs_p, -1d0, p_jump, dpds_known, ring_tang)

    ! U, V, W = G Linv D R f
    rhs_p = dt * rhs_p
    Call gradient(U, V, W, rhs_p)

    ! U, V, W = -R f + G Linv D R f
    U = U - U_reg 
    V = V - V_reg 
    W = W - W_reg 
    Call apply_boundary_conditions(U, V, W)

    ! Sf_ = -E (I -  G Linv D) R f
    Call regT(Eu, U, V, W)
    ib_surface_vector_work = Eu - E1nH_exterior * dudn_jump
    Call cart_to_frame_t(ib_surface_vector_work, Sf_(1 : 2 * nb))
    Sf_(2 * nb + 1 : 3 * nb) = -(ib_surface_vector_work(1 : nb) * normals(1 : nb) &
      + ib_surface_vector_work(nb + 1 : 2 * nb) * normals(nb + 1 : 2 * nb) &
      + ib_surface_vector_work(2 * nb + 1 : 3 * nb) * normals(2 * nb + 1 : 3 * nb))
    Call compute_average_pressure_per_body(pavg, f_(2 * nb + 1 : 3 * nb))
    Call add_per_body_constant(Sf_(2 * nb + 1 : 3 * nb), -pjv_gscale * pavg)
    Sf_(2 * nb + 1 : 3 * nb) = Sf_(2 * nb + 1 : 3 * nb) &
      - E1np - E1nHc_exterior * p_jump - Ecross * dpds_known

  End Subroutine schur

  Subroutine compute_surface_block_rms(vector, weights, component_rms, grouped_rms)
    Real(Int64), Contiguous, Intent(In) :: vector(:), weights(:)
    Real(Int64), Intent(Out) :: component_rms(3), grouped_rms(2)
    Integer(Int32) :: k, n
    Real(Int64) :: component_sum(3), weight_sum

    n = Size(weights)
    If (Size(vector) /= 3 * n) Error Stop 'Surface block vector has the wrong size'
    weight_sum = Sum(weights)
    If (weight_sum <= 0d0) Error Stop 'Surface block weights must have a positive sum'

    component_sum = 0d0
    Do k = 1, n
      component_sum(1) = component_sum(1) + weights(k) * vector(k)**2
      component_sum(2) = component_sum(2) + weights(k) * vector(n + k)**2
      component_sum(3) = component_sum(3) + weights(k) * vector(2 * n + k)**2
    End Do
    component_rms = Sqrt(component_sum / weight_sum)
    grouped_rms(1) = Sqrt(0.5d0 * (component_rms(1)**2 + component_rms(2)**2))
    grouped_rms(2) = component_rms(3)
  End Subroutine compute_surface_block_rms

  Subroutine compute_hybrid_block_tolerances(rhs, weights, relative_tolerance, absolute_tolerance, &
      component_rms, tolerances)
    Real(Int64), Contiguous, Intent(In) :: rhs(:), weights(:)
    Real(Int64), Intent(In) :: relative_tolerance, absolute_tolerance
    Real(Int64), Intent(Out) :: component_rms(3), tolerances(2)
    Real(Int64) :: grouped_rms(2)

    Call compute_surface_block_rms(rhs, weights, component_rms, grouped_rms)
    tolerances = Max(absolute_tolerance, relative_tolerance * grouped_rms)
  End Subroutine compute_hybrid_block_tolerances

  Subroutine compute_scaled_block_residual(residual, weights, tolerances, component_rms, scaled_residual)
    Real(Int64), Contiguous, Intent(In) :: residual(:), weights(:)
    Real(Int64), Intent(In) :: tolerances(2)
    Real(Int64), Intent(Out) :: component_rms(3), scaled_residual(2)
    Real(Int64) :: grouped_rms(2)

    If (Any(tolerances <= 0d0)) Error Stop 'Surface block tolerances must be positive'
    Call compute_surface_block_rms(residual, weights, component_rms, grouped_rms)
    scaled_residual = grouped_rms / tolerances
  End Subroutine compute_scaled_block_residual

  Subroutine bicgstab(bcg_x, bcg_b)
    Integer(Int32) :: iter, report_slot
    Real(Int64), Contiguous, Intent(In) :: bcg_b(:)
    Real(Int64), Contiguous, Intent(InOut) :: bcg_x(:)
    Real(Int64) :: rho_o, rho_n, alpha, om, error, bta, denominator
    Real(Int64) :: component_rms(3), rhs_rms(3), tolerances(2), scaled_residual(2)
    Logical :: record_report

    Call compute_hybrid_block_tolerances(bcg_b, sb, cg_tol, cg_abs_tol, rhs_rms, tolerances)
    record_report = .False.
    If (nmonitor > 0) record_report = Mod(istep, nmonitor) == 0
    report_slot = 0
    If (record_report) Then
      If (cg_report_step /= istep) Then
        cg_report_count = 0
        cg_report_step = istep
      End If
      If (cg_report_count >= cg_report_max_solves) Error Stop 'Too many BiCGSTAB solves in one monitored step'
      cg_report_count = cg_report_count + 1
      report_slot = cg_report_count
      cg_report_stages(report_slot) = rk_step
      cg_report_rhs_rms(:, report_slot) = rhs_rms
      cg_report_tolerances(:, report_slot) = tolerances
    End If

    iter = 0
    Call schur(bcg_r, bcg_x)
    bcg_r = bcg_b - bcg_r
    Call compute_scaled_block_residual(bcg_r, sb, tolerances, component_rms, scaled_residual)
    error = Maxval(scaled_residual)
    If (.Not. ieee_is_finite(error)) Error Stop 'BiCGSTAB produced a non-finite initial residual'
    Call Mpi_bcast (error, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    If (record_report) Then
      cg_report_history(1:3, 0, report_slot) = component_rms
      cg_report_history(4:5, 0, report_slot) = scaled_residual
    End If
    bcg_rhat = bcg_r
    rho_o = 1.d0
    alpha = 1.d0
    om = 1.d0
    bcg_nu = 0.d0
    bcg_p = 0.d0
    Do While ((iter .lt. cg_max_iter) .And. (error > 1d0))
      rho_n = dot_product(bcg_rhat, bcg_r)
      If (.Not. ieee_is_finite(rho_n) .Or. rho_n == 0d0) &
        Error Stop 'BiCGSTAB rho breakdown'
      bta = (rho_n / rho_o) * (alpha / om)
      If (.Not. ieee_is_finite(bta)) Error Stop 'BiCGSTAB beta breakdown'
      rho_o = rho_n
      bcg_p = bcg_r + bta * (bcg_p - om * bcg_nu)
      Call schur(bcg_nu, bcg_p)
      denominator = dot_product(bcg_rhat, bcg_nu)
      If (.Not. ieee_is_finite(denominator) .Or. denominator == 0d0) &
        Error Stop 'BiCGSTAB alpha breakdown'
      alpha = rho_n / denominator
      bcg_h = bcg_x + alpha * bcg_p
      bcg_sv = bcg_r - alpha * bcg_nu
      Call schur(bcg_tv, bcg_sv)
      denominator = dot_product(bcg_tv, bcg_tv)
      If (.Not. ieee_is_finite(denominator)) Error Stop 'BiCGSTAB omega breakdown'
      If (denominator == 0d0) Then
        If (dot_product(bcg_sv, bcg_sv) /= 0d0) Error Stop 'BiCGSTAB omega breakdown'
        om = 0d0
      Else
        om = dot_product(bcg_tv, bcg_sv) / denominator
      End If
      If (.Not. ieee_is_finite(om)) Error Stop 'BiCGSTAB omega breakdown'
      bcg_x = bcg_h + om * bcg_sv
      bcg_r = bcg_sv - om * bcg_tv
      iter = iter + 1
      Call compute_scaled_block_residual(bcg_r, sb, tolerances, component_rms, scaled_residual)
      error = Maxval(scaled_residual)
      If (.Not. ieee_is_finite(error)) Error Stop 'BiCGSTAB produced a non-finite recursive residual'
      Call Mpi_bcast (error, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
      If (record_report) Then
        cg_report_history(1:3, iter, report_slot) = component_rms
        cg_report_history(4:5, iter, report_slot) = scaled_residual
      End If
    End Do

    cg_accum_iter = cg_accum_iter + iter
    cg_nsolves    = cg_nsolves + 1

    If (record_report) cg_report_iterations(report_slot) = iter

    If (iter == cg_max_iter .And. error > 1d0 .And. myid == 0) Then
      Write(*,'(A,I4,A,ES10.3,A)') 'WARNING: BiCGSTAB reached max iterations (', &
        iter, ') but the block-scaled residual is ', error, '!'
    End If
  End Subroutine
End Module projection

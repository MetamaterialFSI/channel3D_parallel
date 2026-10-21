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
    rhs_ib = -Eu + ub

    ! solve for IB forcing
    Call bicgstab(fb, rhs_ib)

    ! U_reg = R f
    aux_surface_scalar = dt * fb(1 : nb)
    Call regu(U_reg, aux_surface_scalar)
    aux_surface_scalar = dt * fb(nb + 1 : 2 * nb)
    Call regv(V_reg, aux_surface_scalar)
    aux_surface_scalar = dt * fb(2 * nb + 1 : 3 * nb)
    Call regw(W_reg, aux_surface_scalar)

    Call apply_boundary_conditions(U_reg, V_reg, W_reg)

    ! rhs_p = D R f
    Call divergence(rhs_p, U_reg, V_reg, W_reg)

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

    Sf_ = 0.d0

    ! U_reg = R f
    aux_surface_scalar = dt * f_(1 : nb)
    Call regu(U_reg, aux_surface_scalar)
    aux_surface_scalar = dt * f_(nb + 1 : 2 * nb)
    Call regv(V_reg, aux_surface_scalar)
    aux_surface_scalar = dt * f_(2 * nb + 1 : 3 * nb)
    Call regw(W_reg, aux_surface_scalar)

    Call apply_boundary_conditions(U_reg, V_reg, W_reg)

    ! rhs_p = D R f
    Call divergence(rhs_p, U_reg, V_reg, W_reg)

    ! rhs_p = Linv D R f
    call solve_poisson_equation(rhs_p, .False.)
    rhs_p = rhs_p / dt



    ! U, V, W = G Linv D R f
    rhs_p = dt * rhs_p
    Call gradient(U, V, W, rhs_p)
    Call apply_boundary_conditions(U, V, W)

    ! U, V, W = -R f + G Linv D R f
    U = U - U_reg 
    V = V - V_reg 
    W = W - W_reg 

    ! Sf_ = -E (I -  G Linv D) R f
    Call regT(Eu, U, V, W)
    Sf_(1 : 3 * nb) = Eu

  End Subroutine schur

  Subroutine compute_surface_block_rms(vector, weights, normal_vectors, tangent_vectors_1, &
      tangent_vectors_2, component_rms, grouped_rms)
    Real(Int64), Contiguous, Intent(In) :: vector(:), weights(:)
    Real(Int64), Contiguous, Intent(In) :: normal_vectors(:), tangent_vectors_1(:), tangent_vectors_2(:)
    Real(Int64), Intent(Out) :: component_rms(3), grouped_rms(2)
    Integer(Int32) :: k, n
    Real(Int64) :: component_sum(3), weight_sum
    Real(Int64) :: vector_x, vector_y, vector_z
    Real(Int64) :: tangent_1_component, tangent_2_component, normal_component

    n = Size(weights)
    If (Size(vector) /= 3 * n) Error Stop 'Surface block vector has the wrong size'
    If (Size(normal_vectors) /= 3 * n .Or. Size(tangent_vectors_1) /= 3 * n .Or. &
        Size(tangent_vectors_2) /= 3 * n) Error Stop 'Surface frame vector has the wrong size'
    weight_sum = Sum(weights)
    If (weight_sum <= 0d0) Error Stop 'Surface block weights must have a positive sum'

    component_sum = 0d0
    Do k = 1, n
      vector_x = vector(k)
      vector_y = vector(n + k)
      vector_z = vector(2 * n + k)
      tangent_1_component = vector_x * tangent_vectors_1(k) + &
        vector_y * tangent_vectors_1(n + k) + vector_z * tangent_vectors_1(2 * n + k)
      tangent_2_component = vector_x * tangent_vectors_2(k) + &
        vector_y * tangent_vectors_2(n + k) + vector_z * tangent_vectors_2(2 * n + k)
      normal_component = vector_x * normal_vectors(k) + vector_y * normal_vectors(n + k) + &
        vector_z * normal_vectors(2 * n + k)
      component_sum(1) = component_sum(1) + weights(k) * tangent_1_component**2
      component_sum(2) = component_sum(2) + weights(k) * tangent_2_component**2
      component_sum(3) = component_sum(3) + weights(k) * normal_component**2
    End Do
    component_rms = Sqrt(component_sum / weight_sum)
    grouped_rms(1) = Sqrt(0.5d0 * (component_rms(1)**2 + component_rms(2)**2))
    grouped_rms(2) = component_rms(3)
  End Subroutine compute_surface_block_rms

  Subroutine compute_hybrid_block_tolerances(rhs, weights, normal_vectors, tangent_vectors_1, &
      tangent_vectors_2, relative_tolerance, absolute_tolerance, component_rms, tolerances)
    Real(Int64), Contiguous, Intent(In) :: rhs(:), weights(:)
    Real(Int64), Contiguous, Intent(In) :: normal_vectors(:), tangent_vectors_1(:), tangent_vectors_2(:)
    Real(Int64), Intent(In) :: relative_tolerance, absolute_tolerance
    Real(Int64), Intent(Out) :: component_rms(3), tolerances(2)
    Real(Int64) :: grouped_rms(2)

    Call compute_surface_block_rms(rhs, weights, normal_vectors, tangent_vectors_1, tangent_vectors_2, &
      component_rms, grouped_rms)
    tolerances = Max(absolute_tolerance, relative_tolerance * grouped_rms)
  End Subroutine compute_hybrid_block_tolerances

  Subroutine compute_scaled_block_residual(residual, weights, normal_vectors, tangent_vectors_1, &
      tangent_vectors_2, tolerances, component_rms, scaled_residual)
    Real(Int64), Contiguous, Intent(In) :: residual(:), weights(:)
    Real(Int64), Contiguous, Intent(In) :: normal_vectors(:), tangent_vectors_1(:), tangent_vectors_2(:)
    Real(Int64), Intent(In) :: tolerances(2)
    Real(Int64), Intent(Out) :: component_rms(3), scaled_residual(2)
    Real(Int64) :: grouped_rms(2)

    If (Any(tolerances <= 0d0)) Error Stop 'Surface block tolerances must be positive'
    Call compute_surface_block_rms(residual, weights, normal_vectors, tangent_vectors_1, tangent_vectors_2, &
      component_rms, grouped_rms)
    scaled_residual = grouped_rms / tolerances
  End Subroutine compute_scaled_block_residual

  Subroutine bicgstab(bcg_x, bcg_b)
    Integer(Int32) :: iter, report_slot
    Real(Int64), Contiguous, Intent(In) :: bcg_b(:)
    Real(Int64), Contiguous, Intent(InOut) :: bcg_x(:)
    Real(Int64) :: rho_o, rho_n, alpha, om, error, bta, denominator
    Real(Int64) :: component_rms(3), rhs_rms(3), tolerances(2), scaled_residual(2)
    Logical :: record_report

    Call compute_hybrid_block_tolerances(bcg_b, sb, normals, tangents_1, tangents_2, cg_tol, &
      cg_abs_tol, rhs_rms, tolerances)
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
    Call compute_scaled_block_residual(bcg_r, sb, normals, tangents_1, tangents_2, tolerances, &
      component_rms, scaled_residual)
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
      Call compute_scaled_block_residual(bcg_r, sb, normals, tangents_1, tangents_2, tolerances, &
        component_rms, scaled_residual)
      error = Maxval(scaled_residual)
      If (.Not. ieee_is_finite(error)) Error Stop 'BiCGSTAB produced a non-finite recursive residual'
      If (om == 0d0 .And. error > 1d0) Error Stop 'BiCGSTAB omega breakdown'
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

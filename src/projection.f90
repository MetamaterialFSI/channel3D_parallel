Module projection

  Use iso_fortran_env, Only : error_unit, Int32, Int64
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
    Call solve_poisson_equation(rhs_p)
    rhs_p = rhs_p / dt

    ! Pnp1 = P* - Linv D R f
    rhs_p = P_interim - rhs_p
    P(2:nxg-1, 2:nyg-1, 2:nzg-1) = rhs_p(2:nxg-1, 2:nyg-1, 2:nzg-1)
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

  Subroutine bicgstab(bcg_x, bcg_b)
    Integer :: j, iter
    Real(Int64), Contiguous, Intent(In) :: bcg_b(:)
    Real(Int64), Contiguous, Intent(InOut) :: bcg_x(:)
    Real(Int64) :: rho_o, rho_n, alpha, om, eps, error, bta

    !initialize
    eps = dot_product(bcg_b, bcg_b)
    If (eps > 0d0) Then
      eps = cg_tol * cg_tol * eps
    Else
      eps = cg_tol * cg_tol
    End If
    iter = 0
    Call schur(bcg_r, bcg_x)
    bcg_r = bcg_b - bcg_r
    error = dot_product( bcg_r, bcg_r)
    Call Mpi_bcast (error, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    bcg_rhat = bcg_r
    rho_o = 1.d0
    alpha = 1.d0
    om = 1.d0
    bcg_nu = 0.d0
    bcg_p = 0.d0
    Do While ((iter .lt. cg_max_iter) .and. (error .ge. eps))
      rho_n = dot_product(bcg_rhat, bcg_r)
      bta = (rho_n / rho_o) * (alpha / om)
      rho_o = rho_n
      bcg_p = bcg_r + bta * (bcg_p - om * bcg_nu)
      Call schur(bcg_nu, bcg_p)
      alpha = rho_n / dot_product(bcg_rhat, bcg_nu)
      bcg_h = bcg_x + alpha * bcg_p
      bcg_sv = bcg_r - alpha * bcg_nu
      Call schur(bcg_tv, bcg_sv)
      om = dot_product( bcg_tv, bcg_sv) / dot_product(bcg_tv, bcg_tv)
      bcg_x = bcg_h + om * bcg_sv
      bcg_r = bcg_sv - om * bcg_tv
      error = dot_product( bcg_r, bcg_r)
      iter = iter + 1
      Call Mpi_bcast (error, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    End Do
    cg_accum_iter = cg_accum_iter + iter
    cg_nsolves    = cg_nsolves + 1
    If ((iter .eq. cg_max_iter) .and. (error .ge. eps) .and. (myid == 0)) Then
      Write(*,'(A,I4,A,ES10.3,A,ES10.3,A)') "WARNING: BiCGSTAB reached max iterations (", &
      iter, ") but residual norm (", sqrt(error), ") exceeds requested tolerance (", cg_tol, ")!"
    End If
  End Subroutine
End Module projection

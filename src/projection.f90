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
    Call solve_poisson_equation(rhs_p)
    rhs_p = rhs_p / dt

    ! Pnp1 = P* - Linv D R f
    rhs_p = P_interim - rhs_p

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
    call solve_poisson_equation(rhs_p)
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

  Subroutine bicgstab(bcg_x, bcg_b)
    Integer :: j, iter
    Real(Int64), Contiguous, Intent(In) :: bcg_b(:)
    Real(Int64), Contiguous, Intent(InOut) :: bcg_x(:)
    Real(Int64) :: rho_o, rho_n, alpha, om, eps, error, bta

    !initialize
    eps = cg_tol * cg_tol
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

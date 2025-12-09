Module projection

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use immersed_boundary_operators
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

    ! save for use in the subsequent projection steps
    P_interim = rhs_p

    Call gradient(U, V, W, rhs_p)

    ! U* = U** - Gp*
    U = U_interim - U
    V = V_interim - V
    W = W_interim - W

  End Subroutine compute_non_IB_projection

  Subroutine compute_IB_projection

    ! - E u* + ub 
    rhs_ib = -regT(U, V, W) + ub

    ! solve for IB forcing
    call bicgstab(fb, rhs_ib)

    ! U_reg = R f
    Call regu(U_reg, fb)
    Call regv(V_reg, fb)
    Call regw(W_reg, fb)
    Call apply_boundary_conditions(U_reg, V_reg, W_reg)

    ! rhs_p = D R f
    Call divergence(rhs_p, U_reg, V_reg, W_reg)
    Call solve_poisson_equation(rhs_p)

    ! Pnp1 = P* - Linv D R f
    rhs_p = P_interim - rhs_p

    ! U, V, W = G Pnp1
    call gradient(U, V, W, rhs_p)

    ! Unp1 = U** - R f - G Pnp1
    U = U_interim - U_reg - U
    V = V_interim - V_reg - V
    W = W_interim - W_reg - W

    Call apply_boundary_conditions(U, V, W)

  End Subroutine compute_IB_projection

  Function schur(f_)
    Implicit None
    Real(Int64), Dimension(3*nb), Intent(In) :: f_
    Real(Int64), Dimension(3*nb) :: schur

    schur = 0.d0

    ! Fibu, Fibv, Fibw = R f
    Call regu(Fibu, f_)
    Call regv(Fibv, f_)
    Call regw(Fibw, f_)
    Call apply_boundary_conditions(Fibu, Fibv, Fibw)

    ! rhs_p = D R f
    Call divergence(rhs_p, Fibu, Fibv, Fibw)

    ! rhs_p = Linv D R f
    call solve_poisson_equation(rhs_p)

    ! U, V, W = G Linv D R f
    Call gradient(U, V, W, rhs_p)
    Call apply_boundary_conditions(U, V, W)

    ! U, V, W = -R f + G Linv D R f
    U = U - Fibu 
    V = V - Fibv 
    W = W - Fibw 

    ! schur = -E (I -  G Linv D) R f
    schur = regT(U, V, W)

  End Function schur

  Subroutine bicgstab( bcg_x, bcg_b)
    Integer :: j, iter
    Real(Int64), Dimension(3*nb), Intent(In) :: bcg_b
    Real(Int64), Dimension(3*nb), Intent(Inout) :: bcg_x
    Real(Int64) :: rho_o, rho_n, alpha, om, eps, error, bta

    !initialize
    error = 1.d0
    eps = cg_tol * cg_tol
    iter = 0
    bcg_r = bcg_b - schur(bcg_x)
    bcg_rhat = bcg_r
    rho_o = 1.d0
    alpha = 1.d0
    om = 1.d0
    bcg_nu = 0.d0
    bcg_p = 0.d0
    Do While ((iter .le. cg_max_iter) .and. (error .ge. eps))
      rho_n = dot_product(bcg_rhat, bcg_r)
      bta = (rho_n / rho_o) * (alpha / om)
      rho_o = rho_n
      bcg_p = bcg_r + bta * (bcg_p - om * bcg_nu)
      bcg_nu = schur(bcg_p)
      alpha = rho_n / dot_product(bcg_rhat, bcg_nu)
      bcg_h = bcg_x + alpha * bcg_p
      bcg_sv = bcg_r - alpha * bcg_nu
      bcg_tv = schur(bcg_sv )
      om = dot_product( bcg_tv, bcg_sv) / dot_product(bcg_tv, bcg_tv)
      bcg_x = bcg_h + om * bcg_sv
      bcg_r = bcg_sv - om * bcg_tv
      error = dot_product( bcg_r, bcg_r)
      iter = iter + 1
      Call Mpi_bcast (error, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    End Do
    cg_accum_iter = cg_accum_iter + (iter - 1)
    If (iter .gt. cg_max_iter .and. myid == 0) Then
      Write(*,*)  "......WARNING, bicgstab used maximum number of iterations (", cg_max_iter, ")"
      Write(*,*)  "......max |residual| = ", Maxval(Abs(bcg_r))
    End If
  End Subroutine
End Module projection

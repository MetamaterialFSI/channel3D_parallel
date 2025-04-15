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
    prev_time = MPI_WTIME()
    rhs_ib = -regT(U, V, W) + ub
    last_time = MPI_WTIME()
    E_1st = E_1st+last_time-prev_time

    ! solve for IB forcing
    prev_time = MPI_WTIME()
    call bicgstab(fb, rhs_ib)
    last_time = MPI_WTIME()
    IB_force = IB_force +last_time-prev_time

    ! U_reg = R f
    prev_time = MPI_WTIME()
    Call regu(U_reg, fb)
    Call regv(V_reg, fb)
    Call regw(W_reg, fb)
    Call apply_boundary_conditions(U_reg, V_reg, W_reg)
    last_time = MPI_WTIME()
    R_1st = R_1st +last_time-prev_time


    ! rhs_p = D R f
    prev_time = MPI_WTIME()
    Call divergence(rhs_p, U_reg, V_reg, W_reg)
    last_time = MPI_WTIME()
    D_1st = D_1st +last_time-prev_time

    prev_time = last_time
    Call solve_poisson_equation(rhs_p)
    last_time = MPI_WTIME()
    IB_possion = IB_possion +last_time-prev_time

    ! Pnp1 = P* - Linv D R f
    prev_time = MPI_WTIME()
    rhs_p = P_interim - rhs_p
    last_time = MPI_WTIME()
    proj_1st = proj_1st +last_time-prev_time

    ! U, V, W = G Pnp1
    prev_time = MPI_WTIME()
    call gradient(U, V, W, rhs_p)
    last_time = MPI_WTIME()
    grad_1st = grad_1st +last_time-prev_time

    ! Unp1 = U** - R f - G Pnp1
    prev_time = MPI_WTIME()
    U = U_interim - U_reg - U
    V = V_interim - V_reg - V
    W = W_interim - W_reg - W

    Call apply_boundary_conditions(U, V, W)
    last_time = MPI_WTIME()
    proj_2nd = proj_2nd +last_time-prev_time

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

  Subroutine bicgstab( x, b)
    Integer :: j, iter, cg_max_iter
    Real(Int64), Dimension(3*nb), Intent(In) :: b
    Real(Int64), Dimension(3*nb), Intent(Inout) :: x
    Real(Int64), Dimension(3*nb) :: r, rhat, p, nu, h, sv, tv
    Real(Int64) :: rho_o, rho_n, alpha, om, eps, error, bta, cgtol

    !initialize
    cgtol=1.E-12
    cg_max_iter=50
    error = 1.d0
    eps = cgtol * cgtol
    iter = 0
    r = b - schur( x)
    rhat = r
    rho_o = 1.d0
    alpha = 1.d0
    om = 1.d0
    nu = 0.d0
    p = 0.d0
    Do While ((iter .le. cg_max_iter) .and. (error .ge. eps))
      rho_n = dot_product(rhat, r)
      bta = (rho_n/rho_o) * (alpha/om)
      rho_o = rho_n
      p = r + bta * (p - om * nu)
      nu = schur(p)
      alpha = rho_n/dot_product(rhat, nu)
      h = x + alpha * p
      sv = r - alpha * nu
      tv = schur(sv )
      om = dot_product( tv, sv)/ dot_product(tv, tv)
      x = h + om * sv
      r = sv - om * tv
      error = dot_product( r, r)
      iter = iter + 1
      Call Mpi_bcast (error, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    End Do
    If (iter .eq. cg_max_iter) Then
      Write(*,*)  "......WARNING, bicgstab used maximum number of iterations"
      Write(*,*)  "......Iterations = ",iter,", residual = ", error
    End If
    ! output iteration
    if ( rk_step==1 )then
      RK1_iter=RK1_iter+iter
    elseif (rk_step ==2) then
      RK2_iter=RK2_iter+iter
    elseif (rk_step==3) then
      RK3_iter=RK3_iter+iter
    end if 
  End Subroutine

End Module projection

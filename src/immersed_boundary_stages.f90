Module immersed_boundary_stages

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use immersed_boundary_operators
  Use equations
  Use projection
  Use boundary_conditions
  Use mass_flow

  ! prevent implicit typing
  Implicit None

Contains

  Subroutine IB_trial_velocity
    Integer(Int32) :: k

    ! Ustar, Vstar, Wstar will have a trial velocity that will be used for the surface stress solve
    ! U,V,W now has the guess and also obeys BCs
    ! save the guess field to Uor1 for later use
    Uor1=U
    Vor1=V
    Wor1=W
  
    Call divergence(rhs_p, U, V, W)

    If ( no_exterior_velocity ) Then
      ! ILM continuity term. ub is the velocity on the minus side of the interface
      Do k=1,nb
        aux_surface_scalar(k) = -normals(k) * ub(k) - &
                                normals(nb + k) * ub(nb + k) - &
                                normals(2 * nb + k) * ub(2 * nb + k)
      End Do
      term = 0.d0
      term(2:nxm,2:nym+1,2:nzm) = regp(aux_surface_scalar)
      Call apply_periodic_bc_x(term, 2) 
      Call apply_periodic_bc_z(term, 2) 
      ! DEBUG output
      ! output_1 = term
      rhs_p = rhs_p - term(2:nxg-1, 2:nyg-1, 2:nzg-1)
    End if

    Call solve_poisson_equation(rhs_p)
    ! rhs_p now has the solution. Assign to Pstar
    Pstar = rhs_p
    Call gradient(U, V, W, Pstar)

    ! U* = U** - Gp*
    Ustar = Uor1 - U
    Vstar = Vor1 - V
    Wstar = Wor1 - W
    Call apply_boundary_conditions(Ustar, Vstar, Wstar)

  end subroutine IB_trial_velocity


  subroutine IB_surface_stress

    rhsib = -regT(Ustar,Vstar,Wstar)

    ! - E u* + 1/2 (ubplus + ubminus) 
    If ( no_exterior_velocity ) Then
      rhsib = rhsib + 0.5 * ub
    Else
      rhsib = rhsib + ub
    End if

    If ( static_body ) Then
      ! use inverted matrix to obtain surface stress if the body is static
      fb = matmul( lumat, rhsib )
    Else
      ! solve for surface stress iteratively if the body is moving
      call bicgstab( fb, rhsib )
    End If

  end subroutine IB_surface_stress

  subroutine IB_pressure
    U=regu(fb)
    V=regv(fb)
    W=regw(fb)
    Call apply_boundary_conditions(U, V, W)

    ! Uor = R f
    Uor=U
    Vor=V
    Wor=W

    ! rhs_p = D R f
    Call divergence(rhs_p, U, V, W)
    Call solve_poisson_equation(rhs_p)

    ! rhs_p = P* - Linv D R f
    rhs_p = Pstar - rhs_p

  end subroutine IB_pressure


  subroutine IB_projection
    U=0.D0
    V=0.D0
    W=0.D0

    ! U, V, W = G Pnp1
    call gradient(U, V, W, rhs_p)

    ! Unp1 = U** - R f - G Pnp1
    U = Uor1 - Uor - U
    V = Vor1 - Vor - V
    W = Wor1 - Wor - W

    Call apply_boundary_conditions(U, V, W)

  end subroutine IB_projection

  subroutine preprocess

    !try to compute the inverse of the surface stress solve operator B
    !output: lumat
    Real(Int64), DIMENSION(3*nb) :: z
    INTEGER :: i, ii,jj
    Real(Int64) ::  Qflow_ref
    Integer :: info, neqns, lda, lwork
    Real(Int64), dimension(:), allocatable :: work
    Integer, dimension(3*nb) :: ipiv_bg

    WRITE(*,*) 'precomputing body matrix for stationary geometry...please be patient'

    DO i=1,3*nb   ! build matrix one column at a time
    z = 0.d0
    z(i) = 1.d0
    lumat(1:3*nb,i) = a_times( z )
    END DO


    !take the inverse of lumat
    info = 0
    lwork = (3*nb)** 2
    allocate( work( lwork ) )
    neqns = 3*nb
    lda = 3*nb
    ipiv_bg = 0
    call dgetrf( neqns, lda, lumat, lda, ipiv_bg, info)
    call dgetri( neqns, lumat, lda, ipiv_bg, work, lwork, info)
    deallocate( work )

    !write out the B matrix of the surface stress solve
    call writestufflumat(lumat)
    !reasoning for the portion below :
    !preprocess ends up changing the velovity variables U,V,W initialized in "init_flow" when it uses a_times
    !so here we reinstate the velocity guesses at t=0

    U=0d0
    V=0d0
    W=0d0

  end subroutine preprocess

  FUNCTION a_times(f_)
    !INPUT: surface stress (f_)
    !OUTPUT: a_times
    !variables changing in the routine: U,V,W,rhs_p
    IMPLICIT NONE
    REAL(Int64), DIMENSION(3*nb), INTENT(IN) :: f_
    REAL(Int64), DIMENSION(3*nb) :: a_times

    a_times = 0.d0

    ! Fibu, Fibv, Fibw = R f
    Fibu=regu(f_)
    Fibv=regv(f_)
    Fibw=regw(f_)
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

    ! a_times = -E (I -  G Linv D) R f
    a_times = regT(U, V, W)

  END FUNCTION a_times

  subroutine bicgstab( x, b)
    !iterative solver
    integer :: j, iter, cg_max_iter
    real(Int64), dimension(3*nb), intent(in) :: b
    real(Int64), dimension(3*nb), intent(inout) :: x
    real(Int64), dimension(3*nb) :: r, rhat, p, nu, h, sv, tv
    real(Int64) :: rho_o, rho_n, alpha, om, eps, err, bta, cgtol
    !initialize
    cgtol=1.E-8
    cg_max_iter=3000
    err = 1.d0
    eps = cgtol * cgtol
    iter = 0
    r = b - a_times( x)
    rhat = r
    rho_o = 1.d0
    alpha = 1.d0
    om = 1.d0
    nu = 0.d0
    p = 0.d0
    do while ((iter .le. cg_max_iter) .and. (err .ge. eps))
      rho_n = dot_product(rhat, r)
      bta = (rho_n/rho_o) * (alpha/om)
      rho_o = rho_n
      p = r + bta * (p - om * nu)
      nu = a_times(p)
      alpha = rho_n/dot_product(rhat, nu)
      h = x + alpha * p
      sv = r - alpha * nu
      tv = a_times(sv )
      om = dot_product( tv, sv)/ dot_product(tv, tv)
      x = h + om * sv
      r = sv - om * tv
      err = dot_product( r, r)
      iter = iter + 1
    end do
    IF (iter.eq.cg_max_iter) THEN
      WRITE(*,*)  "......WARNING, bicg used max iterations"
      WRITE(*,*)    "......itmax = ",iter,", res = ", err
    end if
  end subroutine

  subroutine check_slip
    Real   (Int64) :: max_slip
    ! aux_surface_vector=(regt(Uint,Vint,Wint)+regTHXnf(surface))
    aux_surface_vector=regt(U,V,W)
    aux_surface_vector = aux_surface_vector - ub
    max_slip = Maxval( Abs(aux_surface_vector) )
    Write(*,*) 'Maximum slip          : ',max_slip
    ! aux_surface_vector = regTHXnf(surface)
    ! max_slip = Maxval( Abs(aux_surface_vector) )
    ! Write(*,*) 'regTHXnf          : ',max_slip
  end subroutine check_slip

  subroutine writestufflumat(disp)
    REAL(KIND(0.D0)), dimension(:,:) :: disp
    integer::i,j,k
    OPEN(unit=10001,file="lumat.dat",form="formatted",status="replace")


    DO k=1,3*nb
      DO i=1,3*nb
     WRITE(10001,*)  disp(i,k)
      END DO
    END DO
  
    CLOSE(10001)

  end subroutine writestufflumat

End Module immersed_boundary_stages

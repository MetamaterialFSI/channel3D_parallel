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
  Use operators_FSI
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

    Call gradient(U, V, W, dt * rhs_p)

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
    Call regu(U_reg, dt * fb(1 : nb)             )
    Call regv(V_reg, dt * fb(nb + 1 : 2 * nb)    )
    Call regw(W_reg, dt * fb(2 * nb + 1 : 3 * nb))

    Call apply_boundary_conditions(U_reg, V_reg, W_reg)

    ! rhs_p = D R f
    Call divergence(rhs_p, U_reg, V_reg, W_reg)

    ! rhs_p = Linv D R f
    Call solve_poisson_equation(rhs_p)
    rhs_p = rhs_p / dt

    ! Pnp1 = P* - Linv D R f
    rhs_p = P_interim - rhs_p

    ! U, V, W = G Pnp1
    call gradient(U, V, W, dt * rhs_p)

    ! Unp1 = U** - R f - G Pnp1
    U = U_interim - U_reg - U
    V = V_interim - V_reg - V
    W = W_interim - W_reg - W

    Call apply_boundary_conditions(U, V, W)

  End Subroutine compute_IB_projection
  subroutine compute_IB_FSI_projection

     If ( trim(body_type) ==  'center_wall_deforming_testcase' ) Then
      Call IB_FSI_loop_testcase
     end if
      If ( trim(body_type) ==  'center_wall_deforming_subsurface' ) Then
       Call IB_FSI_loop_subsurface
     end if

    ! U_reg = R f
    Call regu(U_reg, dt * fb(1 : nb)             )
    Call regv(V_reg, dt * fb(nb + 1 : 2 * nb)    )
    Call regw(W_reg, dt * fb(2 * nb + 1 : 3 * nb))

    Call apply_boundary_conditions(U_reg, V_reg, W_reg)

    ! rhs_p = D R f
    Call divergence(rhs_p, U_reg, V_reg, W_reg)

    ! rhs_p = Linv D R f
    Call solve_poisson_equation(rhs_p)
    rhs_p = rhs_p / dt

    ! Pnp1 = P* - Linv D R f
    rhs_p = P_interim - rhs_p

    ! U, V, W = G Pnp1
    call gradient(U, V, W, dt * rhs_p)

    ! Unp1 = U** - R f - G Pnp1
    U = U_interim - U_reg - U
    V = V_interim - V_reg - V
    W = W_interim - W_reg - W

    Call apply_boundary_conditions(U, V, W)

  end subroutine compute_IB_FSI_projection
subroutine IB_FSI_loop_testcase
   real(kind(0.d0)), dimension(3*nb) :: r_chirhs, r_zetarhs, r_c_rhs, zeta_k2, chi_k2, v_tp, r_crhs1
REAL(KIND(0.D0)), Dimension(:), Allocatable :: v_y
real(kind(0.d0)), dimension(nb) :: r_chi, r_zeta, dchi, dchi_1, v_trunc, v_bg, Fint
integer :: j
real(kind(0.d0)) :: sum

Allocate(v_y(nb))
Fint=0.d0
r_crhs1=0.d0
err_FSI = 10.0
tol_FSI = 1.e-5
iter_FSI=0

call khatmatrix_testcase(sol_mat,Mmat_testcase,Kmat_testcase,Cmat_testcase)

Ustarfsi=U
Vstarfsi=V
Wstarfsi=W

do while ( (err_FSI .ge. tol_FSI) .and. (iter_FSI .le. 1000) )

      call setup_IB_operators

      ! get rhs contribution from the guess
       rhsib = -regT(Ustarfsi, Vstarfsi, Wstarfsi)
      ! get rhs contribution from structural motion
      r_chi=zeta+((2/dt_fsi)*chi)-((2/dt_fsi)*chi_k)+zeta_k
      r_chirhs=-Itilde(r_chi)
      r_zeta = matmul(Mmat_testcase, ( zetadot+ (4.0d0/dt_fsi)*zeta + (4.0d0/(dt_fsi*dt_fsi))*(chi - chi_k) ))+ matmul(Cmat_testcase, ( zeta+ (2.0d0/dt_fsi)*(chi - chi_k) ))+ F_bf- matmul(Kmat_testcase, chi_k)
      r_zetarhs = (2/dt_fsi)*Itilde(matmul(sol_mat,r_zeta))
      r_crhs1= (Itilde(zeta_k))
      r_c_rhs=r_crhs1
      rhsib=rhsib+r_c_rhs+r_zetarhs+r_chirhs


      call bicgstab_fsi_testcase( fb, rhsib)

           


      dchi_1=KhatinvQItildeprimeW(fb)
      dchi= (matmul(sol_mat,r_zeta))+dchi_1
      chi_k = chi_k + dchi
      zeta_k = -zeta + ((2.0/dt_fsi) * (chi_k - chi))
      zetadot_k = ((4.0/(dt_fsi*dt_fsi))* (chi_k - chi)) - ((4.0/dt_fsi)*zeta) - zetadot

      if (maxval(abs( chi_k)) .ge. 1.d-13) then
            err_FSI = maxval(abs(dchi))/maxval(abs(chi_k))
      else
            err_FSI = maxval(abs(dchi))
      end if

      chi_k2=Itilde(chi_k)
      xb = xbref +  chi_k2(1:nb)
      yb=  ybref +  chi_k2(nb+1:2*nb)
      zb=  zbref +  chi_k2((2*nb)+1:3*nb)
      zeta_k2=Itilde(zeta_k)
      iter_FSI=iter_FSI+1

end do

xb=xb
yb=yb
zb=zb

  end subroutine IB_FSI_loop_testcase

subroutine IB_FSI_loop_subsurface
    real(kind(0.d0)), dimension(3*nb) :: r_chirhs, r_zetarhs, r_c_rhs,  zeta_k2,chi_k2,v_tp,r_crhs1
    REAL(KIND(0.D0)), Dimension(:), Allocatable ::v_y
    real(kind(0.d0)), dimension(nblocks) ::r_chi,r_zeta,dchi,dchi_1,v_trunc,v_bg, Fint
    integer :: i,j
    real(kind(0.d0))::sum

    Allocate(v_y(nb))
    Fint=0.d0
    r_crhs1=0.d0
    err_FSI = 10.0
    tol_FSI = 1.e-5
    iter_FSI=0

    call khatmatrix_subsurface(sol_mat,Mmat,Kmat)

              Ustarfsi=U
              Vstarfsi=V
              Wstarfsi=W
      do while ( (err_FSI .ge. tol_FSI) .and. (iter_FSI .le. 1000) )

              call setup_IB_operators
              !get rhs contribution from the guess
              rhsib = -regT(Ustarfsi, Vstarfsi, Wstarfsi)

              !get rhs contribution from structural motion
              r_chi=zeta+((2/dt_fsi)*chi)-((2/dt_fsi)*chi_k)+zeta_k

              r_chirhs=-ItildeGprime_subsurface(r_chi)

              r_zeta = matmul(Mmat, ( zetadot+ ((4.0/dt_fsi)* zeta) + ((4.0/(dt_fsi*dt_fsi))*(chi - chi_k)) )) +F_bf - matmul(Kmat,chi_k)
              r_zetarhs = (2/dt_fsi)*ItildeGprime_subsurface(matmul(sol_mat,r_zeta))

              r_crhs1= (ItildeGprime_subsurface(zeta_k))
              r_c_rhs=r_crhs1
              rhsib=rhsib+r_c_rhs+r_zetarhs+r_chirhs
             
              call bicgstab_fsi_subsurface( fb, rhsib)
              
              dchi_1=KhatinvQItildeprimeW_subsurface(fb)
          
              dchi= (matmul(sol_mat,r_zeta))+dchi_1
              chi_k = chi_k + dchi
              zeta_k = -zeta + ((2.0/dt_fsi) * (chi_k - chi))
              zetadot_k = ((4.0/(dt_fsi*dt_fsi))* (chi_k - chi)) - ((4.0/dt_fsi)*zeta) - zetadot

              if (maxval(abs( chi_k)) .ge. 1.d-13) then
                    err_FSI = maxval(abs(dchi))/maxval(abs(chi_k))
                else
                    err_FSI = maxval(abs(dchi))
                end if
            chi_k2=ItildeGprime_subsurface(chi_k)

            xb = xbref +  chi_k2(1:nb)
            yb=  ybref +  chi_k2(nb+1:2*nb)
            zb=  zbref +  chi_k2((2*nb)+1:3*nb)

            zeta_k2=ItildeGprime_subsurface(zeta_k)
            iter_FSI=iter_FSI+1;


      end do
      xb=xb
      yb=yb
      zb=zb


  end subroutine IB_FSI_loop_subsurface

  Function schur(f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(3 * nb) :: schur

    schur = 0.d0

    ! U_reg = R f
    Call regu(U_reg, dt * f_(1 : nb)             )
    Call regv(V_reg, dt * f_(nb + 1 : 2 * nb)    )
    Call regw(W_reg, dt * f_(2 * nb + 1 : 3 * nb))

    Call apply_boundary_conditions(U_reg, V_reg, W_reg)

    ! rhs_p = D R f
    Call divergence(rhs_p, U_reg, V_reg, W_reg)

    ! rhs_p = Linv D R f
    call solve_poisson_equation(rhs_p)
    rhs_p = rhs_p / dt

    ! U, V, W = G Linv D R f
    Call gradient(U, V, W, dt * rhs_p)
    Call apply_boundary_conditions(U, V, W)

    ! U, V, W = -R f + G Linv D R f
    U = U - U_reg
    V = V - V_reg
    W = W - W_reg

    ! schur = -E (I -  G Linv D) R f
    schur = regT(U, V, W)

  End Function schur

  Subroutine bicgstab( bcg_x, bcg_b)
    Integer :: j, iter
    Real(Int64), Dimension(3 * nb), Intent(In) :: bcg_b
    Real(Int64), Dimension(3 * nb), Intent(Inout) :: bcg_x
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
    cg_accum_iter = cg_accum_iter + iter
    If (iter .gt. cg_max_iter .and. myid == 0) Then
      Write(*,*)  "......WARNING, bicgstab used maximum number of iterations (", cg_max_iter, ")"
      Write(*,*)  "......max |residual| = ", Maxval(Abs(bcg_r))
    End If
  End Subroutine

Subroutine bicgstab_fsi_testcase( bcg_x, bcg_b)
    Integer :: j, iter
    Real(Int64), Dimension(3*nb), Intent(In) :: bcg_b
    Real(Int64), Dimension(3*nb), Intent(Inout) :: bcg_x
    Real(Int64) :: rho_o, rho_n, alpha, om, eps, error, bta

    !initialize
    !bcg_x=0.d0

    error = 1.d0
    eps = cg_tol * cg_tol
    iter = 0
    bcg_r = bcg_b - schur(bcg_x)-b_times_testcase(bcg_x)
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
      bcg_nu = schur(bcg_p)+b_times_testcase(bcg_p)
      alpha = rho_n / dot_product(bcg_rhat, bcg_nu)
      bcg_h = bcg_x + alpha * bcg_p
      bcg_sv = bcg_r - alpha * bcg_nu
      bcg_tv = schur(bcg_sv )+b_times_testcase(bcg_sv)
      om = dot_product( bcg_tv, bcg_sv) / dot_product(bcg_tv, bcg_tv)
      bcg_x = bcg_h + om * bcg_sv
      bcg_r = bcg_sv - om * bcg_tv
      error = dot_product( bcg_r, bcg_r)
      iter = iter + 1
      Call Mpi_bcast (error, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    End Do
    If (iter .gt. cg_max_iter .and. myid == 0) Then
      Write(*,*)  "......WARNING, bicgstab used maximum number of iterations (", cg_max_iter, ")"
      Write(*,*)  "......max |residual| = ", Maxval(Abs(bcg_r))
    End If
    
  End Subroutine

Subroutine bicgstab_fsi_subsurface( bcg_x, bcg_b)
    Integer :: j, iter
    Real(Int64), Dimension(3*nb), Intent(In) :: bcg_b
    Real(Int64), Dimension(3*nb), Intent(Inout) :: bcg_x
    Real(Int64) :: rho_o, rho_n, alpha, om, eps, error, bta

    !initialize
    !bcg_x=0.d0

    error = 1.d0
    eps = cg_tol * cg_tol
    iter = 0
    bcg_r = bcg_b - schur(bcg_x)-b_times_subsurface(bcg_x)
    
  
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
      bcg_nu = schur(bcg_p)+b_times_subsurface(bcg_p)
      alpha = rho_n / dot_product(bcg_rhat, bcg_nu)
      bcg_h = bcg_x + alpha * bcg_p
      bcg_sv = bcg_r - alpha * bcg_nu
      bcg_tv = schur(bcg_sv )+b_times_subsurface(bcg_sv)
      om = dot_product( bcg_tv, bcg_sv) / dot_product(bcg_tv, bcg_tv)
      bcg_x = bcg_h + om * bcg_sv
      bcg_r = bcg_sv - om * bcg_tv
      error = dot_product( bcg_r, bcg_r)
      iter = iter + 1
      Call Mpi_bcast (error, 1, MPI_real8, 0, MPI_COMM_WORLD, ierr)
    End Do
    If (iter .gt. cg_max_iter .and. myid == 0) Then
      Write(*,*)  "......WARNING, bicgstab used maximum number of iterations (", cg_max_iter, ")"
      Write(*,*)  "......max |residual| = ", Maxval(Abs(bcg_r))
    End If
    
  End Subroutine





End Module projection

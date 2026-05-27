!---------------------------------------------!
!     Module for temporal integration         !
!---------------------------------------------!
Module time_integration

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use equations
  Use projection
  Use boundary_conditions
  Use mass_flow
  Use immersed_boundary_geometry
  Use immersed_boundary_operators
  Use heaviside

  ! prevent implicit typing
  Implicit None

Contains

  !-----------------------------------------------!
  !                Explicit Euler                 !
  !-----------------------------------------------!
  Subroutine compute_time_step_Euler

    ! equivalent to last rk steps
    rk_step = 3

    ! save current step
    Uo = U
    Vo = V
    Wo = W

    ! compute rhs for U
    Call compute_rhs_u(Uo,Vo,Wo,rhs_uo)

    ! Advance U interior points
    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*rhs_uo

    ! compute rhs for V
    Call compute_rhs_v(Uo,Vo,Wo,rhs_vo)

    ! Advance V interior points
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*rhs_vo

    ! compute rhs for W
    Call compute_rhs_w(Uo,Vo,Wo,rhs_wo)

    ! Advance W interior points
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*rhs_wo

    ! Advance time
    t = t + dt

    ! boundary conditions
    Call apply_boundary_conditions(U, V, W)

    ! projection step
    Call compute_non_IB_projection

    ! boundary conditions
    Call apply_boundary_conditions(U, V, W)

    ! compute mean pressure gradient for constant mass flow in x
    If ( x_mass_cte == 1 ) Then
       Call compute_dPx_for_constant_mass_flow(U,dPdx)
       U(2:nx-1,2:nyg-1,2:nzg-1) = U(2:nx-1,2:nyg-1,2:nzg-1) + dPdx
       dPdx = dPdx/dt ! to be used later by rhs_*
       Call apply_boundary_conditions(U, V, W)
    End If

    ! compute mean pressure gradient for constant mass flow in y
    If ( y_mass_cte == 1 ) Then
       Call compute_dPy_for_constant_mass_flow(V,dPdy)
       V(2:nxg-1,1:ny,2:nzg-1) = V(2:nxg-1,1:ny,2:nzg-1) + dPdy
       dPdy = dPdy/dt ! to be used later by rhs_*
       Call apply_boundary_conditions(U, V, W)
    End If

  End Subroutine compute_time_step_Euler

  !-----------------------------------------------!
  !          Explicit Runge-Kutta 3 steps         !
  !-----------------------------------------------!
  Subroutine compute_time_step_RK3

    Real(Int64) :: to
 real(kind(0.d0)) ::  pi
  pi = 4.0d0*atan(1.0d0)
    ! save previous state
    to = t
    Uo = U
    Vo = V
    Wo = W

    ! Vw(:,1,:) = -V(:,24,:)
    ! Vw(:,2,:) = -V(:,ny-23,:)

    ! step 1
    rk_step = 1
    Call compute_rhs_u(U,V,W,Fu1)
    Call compute_rhs_v(U,V,W,Fv1)
    Call compute_rhs_w(U,V,W,Fw1)

    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*rk_coef(1,1)*Fu1
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*rk_coef(1,1)*Fv1
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*rk_coef(1,1)*Fw1
    t = to + rk_t(rk_step)*dt

    ! update body point positions and velocities if body is moving
    If ( moving_body ) Then
      call setup_IB_geometry
      call setup_IB_operators
      call compute_heaviside
    End If

    Call apply_boundary_conditions(U, V, W)
    Call compute_non_IB_projection
    If ( trim(body_type) /= 'none' .and. functionality==0) Then
      Call apply_boundary_conditions(U, V, W)
      Call compute_IB_projection
    End If
    If ( trim(body_type) ==  'center_wall_deforming_testcase' .and. functionality==1 ) Then
              call  compute_initial_acceleration_testcase
              Call apply_boundary_conditions(U, V, W)
              Call compute_IB_FSI_projection
    End If

    If ( trim(body_type) ==  'center_wall_deforming_subsurface' .and. functionality==1 ) Then
              call compute_initial_acceleration_subsurface
              Call apply_boundary_conditions(U, V, W)
              Call compute_IB_FSI_projection

    End If
    Call apply_boundary_conditions(U, V, W)

    ! Vw(:,1,:) = -V(:,24,:)
    ! Vw(:,2,:) = -V(:,ny-23,:)

    ! step 2
    rk_step = 2
    Call compute_rhs_u(U,V,W,Fu2)
    Call compute_rhs_v(U,V,W,Fv2)
    Call compute_rhs_w(U,V,W,Fw2)

    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + dt*( rk_coef(2,1)*Fu1 + rk_coef(2,2)*Fu2 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + dt*( rk_coef(2,1)*Fv1 + rk_coef(2,2)*Fv2 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + dt*( rk_coef(2,1)*Fw1 + rk_coef(2,2)*Fw2 )
    t = to + rk_t(rk_step)*dt

    ! update body point positions and velocities if body is moving
    If ( moving_body ) Then
      call setup_IB_geometry
      call setup_IB_operators
      call compute_heaviside
    End If

    Call apply_boundary_conditions(U, V, W)
    Call compute_non_IB_projection
    If ( trim(body_type) /= 'none'.and. functionality==0  ) Then
      Call apply_boundary_conditions(U, V, W)
      Call compute_IB_projection
    End If
       If ( trim(body_type) /= 'none' .and. functionality==1 ) Then
      dt_fsi= rk_t(2)*dt
      Call compute_IB_FSI_projection
    End If
    Call apply_boundary_conditions(U, V, W)

    ! Vw(:,1,:) = -V(:,24,:)
    ! Vw(:,2,:) = -V(:,ny-23,:)

    ! step 3
    rk_step = 3
    Call compute_rhs_u(U,V,W,Fu3)
    Call compute_rhs_v(U,V,W,Fv3)
    Call compute_rhs_w(U,V,W,Fw3)

    U(2:nx-1,2:nyg-1,2:nzg-1) = Uo(2:nx-1,2:nyg-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fu1 + rk_coef(3,2)*Fu2 + rk_coef(3,3)*Fu3 )
    V(2:nxg-1,2:ny-1,2:nzg-1) = Vo(2:nxg-1,2:ny-1,2:nzg-1) + &
         dt*( rk_coef(3,1)*Fv1 + rk_coef(3,2)*Fv2 + rk_coef(3,3)*Fv3 )
    W(2:nxg-1,2:nyg-1,2:nz-1) = Wo(2:nxg-1,2:nyg-1,2:nz-1) + &
         dt*( rk_coef(3,1)*Fw1 + rk_coef(3,2)*Fw2 + rk_coef(3,3)*Fw3 )
    t = to + rk_t(rk_step)*dt

    ! update body point positions and velocities if body is moving
    If ( moving_body ) Then
      call setup_IB_geometry
      call setup_IB_operators
      call compute_heaviside
    End If

    Call apply_boundary_conditions(U, V, W)
    Call compute_non_IB_projection
    If ( trim(body_type) /= 'none'.and. functionality==0 ) Then
      Call apply_boundary_conditions(U, V, W)
      Call compute_IB_projection
    End If
    If ( trim(body_type) /= 'none' .and. functionality==1 ) Then
       dt_fsi= rk_t(3)*dt
      Call compute_IB_FSI_projection
         chi=chi_k
         zeta=zeta_k
         zetadot=zetadot_k
    End If
    Call apply_boundary_conditions(U, V, W)

    ! compute mean pressure gradient for constant mass flow in x
    If ( x_mass_cte == 1 ) Then
       Call compute_dPx_for_constant_mass_flow(U,dPdx)
       U(2:nx-1,2:nyg-1,2:nzg-1) = U(2:nx-1,2:nyg-1,2:nzg-1) + dPdx
       dPdx = dPdx/dt ! to be used later by rhs_*
       Call apply_boundary_conditions(U, V, W)
    End If

    ! compute mean pressure gradient for constant mass flow in y
    If ( y_mass_cte == 1 ) Then
       Call compute_dPy_for_constant_mass_flow(V,dPdy)
       V(2:nxg-1,2:ny-1,2:nzg-1) = V(2:nxg-1,2:ny-1,2:nzg-1) + dPdy
       dPdy = dPdy/dt ! to be used later by rhs_*
       Call apply_boundary_conditions(U, V, W)
    End If

    ! compute mean pressure gradient for constant mass flow in z
    If ( z_mass_cte == 1 ) Then
       Call compute_dPz_for_constant_mass_flow(W,dPdz)
       W(2:nxg-1,2:nyg-1,2:nz-1) = W(2:nxg-1,2:nyg-1,2:nz-1) + dPdz
       dPdz = dPdz/dt ! to be used later by rhs_*
       Call apply_boundary_conditions(U, V, W)
    End If

  End Subroutine compute_time_step_RK3

  !-----------------------------------------------!
  !            compute dt based on CFL            !
  !-----------------------------------------------!
  ! NOTE: add rotating and eddy viscosity CFL
  Subroutine compute_dt

    Integer(Int32) :: i, j, k
    Real   (Int64) :: lUmax, lVmax, lWmax, dt_local
    Real   (Int64) :: dt_conv_u, dt_conv_v, dt_conv_w, dt_conv
    Real   (Int64) :: dt_vis_u, dt_vis_v, dt_vis_w, dt_vis
    Real   (Int64) :: dt_max

    ! convective time step
    lUmax = 0d0
    lVmax = 0d0
    lWmax = 0d0
    Do i=2,nxg-1
       Do j=2,nyg-1
          Do k=2,nzg-1
             lUmax = Max( lUmax,(xg(i+1)-xg(i))/Abs(U(i,j,k)) )
             lVmax = Max( lVmax,(yg(j+1)-yg(j))/Abs(V(i,j,k)) )
             lWmax = Max( lWmax,(zg(k+1)-zg(k))/Abs(W(i,j,k)) )
          End Do
       End Do
    End Do

    dt_conv_u = CFL*lUmax
    dt_conv_v = CFL*lVmax
    dt_conv_w = CFL*lWmax

    dt_conv = Minval( (/dt_conv_u,dt_conv_v,dt_conv_w/) )

    ! viscous time step
    dt_vis_u = CFL*dxmin**2d0/nu
    dt_vis_v = CFL*dymin**2d0/nu
    dt_vis_w = CFL*dzmin**2d0/nu

    dt_vis = Minval( (/dt_vis_u,dt_vis_v,dt_vis_w/) )

    ! time step
    dt_local = Min ( dt_conv,dt_vis )

    ! compute global minimum and communicate results to all processors
    Call MPI_Allreduce(dt_local,dt,1,MPI_real8,MPI_min,MPI_COMM_WORLD,ierr)

    ! time step limiter
    dt_max = 1d-1
    dt     = Min( dt, dt_max )
    If ( CFL<0 ) Then
       dt = -CFL
    End If

   End Subroutine compute_dt



Subroutine compute_initial_acceleration_subsurface

    integer :: info, neqns, lda, lwork
    real(kind(0.d0)), dimension(:), allocatable :: work
    integer, dimension(nblocks) :: ipiv_bg
    real(kind(0.d0)), dimension(nblocks) ::Fint

       If ( trim(body_type) ==  'center_wall_deforming_subsurface' .and. functionality==1 ) Then
           if (istep .eq. 1) then
            print *, "computing consistent initial body acceleration..."
            info = 0
            lwork = (nblocks)** 2
            allocate( work( lwork ) )
            neqns = nblocks
            lda = nblocks
            ipiv_bg = 0

            call dgetrf( neqns, lda, Mmat, lda, ipiv_bg, info)
            call dgetri( neqns, Mmat, lda, ipiv_bg, work, lwork, info)
            deallocate( work )
             Fint=0.0
             F_bf=0.0
             F_bf(1)=0.1
             Fint=matmul(Kmat,chi_k)
            zetadot_k = matmul( Mmat, -Fint +F_bf )
            Call mass_matrix(Mmat)

            print *, "initial acceleration computed"

        end if
              dt_fsi= rk_t(1)*dt
              chi=chi_k
              zeta=zeta_k
              zetadot=zetadot_k
        end if

end subroutine compute_initial_acceleration_subsurface


Subroutine compute_initial_acceleration_testcase

    integer :: info, neqns, lda, lwork
    real(kind(0.d0)), dimension(:), allocatable :: work
    integer, dimension(nb) :: ipiv_bg
    real(kind(0.d0)), dimension(nb) ::Fint
    integer :: i, j, flag


     If ( trim(body_type) ==  'center_wall_deforming_testcase' .and. functionality==1 ) Then
            if (istep .eq. 1) then
            print *, "computing consistent initial body acceleration..."
            info = 0
            lwork = (nb)**2
            allocate(work(lwork))
            neqns = nb
            lda = nb
            ipiv_bg = 0
            call dgetrf(neqns, lda, Mmat_testcase, lda, ipiv_bg, info)
            call dgetri(neqns, Mmat_testcase, lda, ipiv_bg, work, lwork, info)
            deallocate(work)
            Fint = 0.0d0
            F_bf = 0.0d0
            do i = 1, nb
                F_bf(i) =  0.1*sin((2.0d0*pi*xb(i))/Lxp) * sin((2.0d0*pi*zb(i))/Lzp)
                 !F_bf(i) =  ((xb(i)))
            end do
            Fint = matmul(Kmat_testcase, chi_k)
            zeta_k=0 ! intial non-zero velocity due to forcing
            zetadot_k = matmul(Mmat_testcase, -Fint + F_bf)
            call mass_matrix_testcase(Mmat_testcase)
            print *, "initial acceleration computed"
            end if
              dt_fsi= rk_t(1)*dt
              chi=chi_k
              zeta=zeta_k
              zetadot=zetadot_k
      end if
  
end subroutine compute_initial_acceleration_testcase



End Module time_integration

!--------------------------------------------------------!
! Module to compute right-hand side of Navier-Stokes eq. !
!--------------------------------------------------------!
Module equations

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global,          Only : x, xm, xg, y, ym, yg, z, zm, zg, term_1, & 
                              term_2, nx, nxg, ny, nyg, nz, nzg, nu,   & 
                              dPdx, dPdy, dPdz, yg_m, exterior_pressure_gradient, &
                              Hu_interior, Hv_interior, Hw_interior
  Use interpolation
  
  ! prevent implicit typing
  Implicit None
  
Contains

  !--------------------------------------------------------------------!
  !                      Compute RHS for du/dt                         !
  !                                                                    !
  ! du/dt = -du^2/dx - duv/dy - duw/dz + div(nu grad(u))               ! 
  !--------------------------------------------------------------------!
  Subroutine compute_rhs_u(U_,V_,W_,rhs_u)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nx-1,2:nyg-1,2:nzg-1), Intent(Out) :: rhs_u

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dx_1, dx_2, dx_3, maxerr
    Real   (Int64) :: dy_1, dy_2, dy_3
    Real   (Int64) :: dz_1, dz_2, dz_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2 

    !-------------compute convective term_2s----------------!

    ! 1)---------compute -du^2/dx--------------

    ! interpolate u in x (faces to centers)
    Call interpolate_x(U_,term_1(2:nxg-1,1:nyg,1:nzg)) 

    ! u^2
    term_1 = term_1**2d0

    ! -du^2/dx
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             term_2(i,j,k) = -(term_1(i+1,j,k)-term_1(i,j,k))/( xg(i+1) - xg(i) )
          End Do
       End Do
    End Do

    rhs_u = term_2(2:nx-1,2:nyg-1,2:nzg-1)

    ! 2)-----------compute -duv/dy--------------

    ! interpolate u in y (centers to faces)
    Call interpolate_y(U_,term_1(1:nx,1:ny,1:nzg)) 

    ! interpolate v in x (centers to faces)
    Call interpolate_x(V_,term_2(1:nx,1:ny,1:nzg)) 

    ! uv at faces
    term_1 = term_1*term_2

    ! -duv/dy (this derivative goes to ym)
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1
             term_2(i,j,k) = -( term_1(i,j,k) - term_1(i,j-1,k) )/( y(j) - y(j-1) )
          End Do
       End Do
    End Do

    rhs_u = rhs_u + term_2(2:nx-1,2:nyg-1,2:nzg-1)

    ! 3)--------------compute -duw/dz--------------

    ! interpolate u in z (centers to faces)
    Call interpolate_z(U_,term_1(1:nx,1:nyg,1:nz)) 
 
    ! interpolate w in x (centers to faces)
    Call interpolate_x(W_,term_2(1:nx,1:nyg,1:nz)) 

    ! uw at faces
    term_1 = term_1*term_2

    ! -duw/dz
    Do k=2,nzg-1
       Do j=2,nyg-1
          Do i=2,nx-1       
             term_2(i,j,k) = -( term_1(i,j,k) - term_1(i,j,k-1) )/( z(k) - z(k-1) )
          End Do
       End Do
    End Do

    rhs_u = rhs_u + term_2(2:nx-1,2:nyg-1,2:nzg-1)

    !--------------compute viscous term_2s------------!

    ! 4)--------compute div( nu grad(u) )------------

    ! first derivation in y: du/dy goes to yg_m(1:ny)
    Do i = 1, nx
       Do k = 1, nzg
          term_1(i,1:nyg-1,k) = ( U_(i,2:nyg,k) - U_(i,1:nyg-1,k) )/( yg(2:nyg) - yg(1:nyg-1) ) 
       End Do
    End Do

    Do k=2,nzg-1
       
       ! first derivation in z (U in centers)
       dz_1 = zg(  k) - zg(k-1) 
       dz_2 = zg(k+1) - zg(k  )
       ! second derivation in z (U at faces)
       dz_3 = z(k) - z(k-1)
       
       Do j=2,nyg-1

          ! second derivation in y (U at faces)
          dy_3 = y(j) - y(j-1)

          Do i=2,nx-1
             
             ! first derivation in x (U at faces)
             dx_1 = x(  i) - x(i-1) 
             dx_2 = x(i+1) - x(i  )
             ! second derivation in x (U at centers)
             dx_3 = xg(i+1) - xg(i)   

             ! viscous term_2
             term_2(i,j,k) = 1d0/dx_3*(nu*1d0/dx_2*(U_(i+1,j,k)-U_(i,j,k)) - nu*1d0/dx_1*(U_(i,j,k)-U_(i-1,j,k)) ) + & ! d^2u/dx^2
                             1d0/dy_3*(nu*term_1(i,j,k)                    - nu*term_1(i,j-1,k) )                  + & ! d^2u/dy^2
                             1d0/dz_3*(nu*1d0/dz_2*(U_(i,j,k+1)-U_(i,j,k)) - nu*1d0/dz_1*(U_(i,j,k)-U_(i,j,k-1)) )     ! d^2u/dz^2

          End Do
       End Do
    End Do

    rhs_u = rhs_u + term_2(2:nx-1,2:nyg-1,2:nzg-1)

    !--------------Constant pressure gradient-------------!
    If (exterior_pressure_gradient == .True.) Then
      rhs_u = rhs_u + dPdx * Hu_interior(2:nx-1,2:nyg-1,2:nzg-1)
    Else
      rhs_u = rhs_u + dPdx
    End If

  End Subroutine compute_rhs_u

  !--------------------------------------------------------------------!
  !                       Compute RHS for dv/dt                        !
  !                                                                    !
  ! dv/dt = -duv/dx - dv^2/dy - dvw/dz + div(nu grad(v))               ! 
  !--------------------------------------------------------------------!
  Subroutine compute_rhs_v(U_,V_,W_,rhs_v)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nxg-1,2:ny-1,2:nzg-1), Intent(Out) :: rhs_v

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: dx_1, dx_2, dx_3, maxerr
    Real   (Int64) :: dy_1, dy_2, dy_3
    Real   (Int64) :: dz_1, dz_2, dz_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2 

    !-------------compute convective term_2s----------------!

    ! 1)---------compute -dv^2/dy--------------

    ! interpolate v in y (faces to centers)
    Call interpolate_y(V_,term_1(1:nxg,2:nyg-1,1:nzg)) 

    ! compute v^2 at centers
    term_1 = term_1**2d0

    ! -dv^2/dy ( this derivative goes to yg_m(2:ny-1) ) 
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             term_2(i,j,k) = -(term_1(i,j+1,k)-term_1(i,j,k))/( yg(j+1) - yg(j) )
          End Do
       End Do
    End Do
    
    rhs_v = term_2(2:nxg-1,2:ny-1,2:nzg-1)

    ! 2)-----------compute -duv/dx--------------

    ! interpolate u in y (centers to faces)
    Call interpolate_y(U_,term_1(1:nx,1:ny,1:nzg)) 

    ! interpolate v in x (centers to faces)
    Call interpolate_x(V_,term_2(1:nx,1:ny,1:nzg)) 

    ! uv
    term_1 = term_1*term_2

    ! -duv/dx
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             term_2(i,j,k) = -( term_1(i,j,k) - term_1(i-1,j,k) )/( x(i) - x(i-1) )
          End Do
       End Do
    End Do

    rhs_v = rhs_v + term_2(2:nxg-1,2:ny-1,2:nzg-1)

    ! 3)--------------compute -dvw/dz--------------

    ! interpolate v in z (centers to faces)
    Call interpolate_z(V_,term_1(1:nxg,1:ny,1:nz)) 

    ! interpolate w in y (centers to faces)
    Call interpolate_y(W_,term_2(1:nxg,1:ny,1:nz)) 

    ! vw
    term_1 = term_1*term_2

    ! -dvw/dz
    Do k=2,nzg-1
       Do j=2,ny-1
          Do i=2,nxg-1
             term_2(i,j,k) = -( term_1(i,j,k) - term_1(i,j,k-1) )/( z(k) - z(k-1) )
          End Do
       End Do
    End Do

    rhs_v = rhs_v + term_2(2:nxg-1,2:ny-1,2:nzg-1)

    !--------------compute viscous term_2s------------!
    
    ! 4)-----------compute div( nu grad(v) )----------

    ! interpolate eddy viscosity to faces

    ! second order remain, no need to interpolate
    Do k=2,nzg-1
       
       ! first derivation in z (V at centers)
       dz_1 = zg(  k) - zg(k-1) 
       dz_2 = zg(k+1) - zg(k  )
       ! second derivation in z (V at faces)
       dz_3 = z(k) - z(k-1)
       
       Do j=2,ny-1

          ! first derivation in y (V at faces)
          dy_1 = y(  j) - y(j-1) 
          dy_2 = y(j+1) - y(j  )
          ! second derivation in y (V at centers)
          dy_3 = yg(j+1) - yg(j)

          Do i=2,nxg-1
             
             ! first derivation in x (V at centers)
             dx_1 = xg(  i) - xg(i-1) 
             dx_2 = xg(i+1) - xg(i  )
             ! second derivation in x (V at faces)
             dx_3 = x(i) - x(i-1)   

             ! viscous term_2
             term_2(i,j,k) = 1d0/dx_3*(nu*1d0/dx_2*(V_(i+1,j,k) - V_(i,j,k)) - nu*1d0/dx_1*(V_(i,j,k)-V_(i-1,j,k)) ) + & ! d^2v/dx^2
                             1d0/dy_3*(nu*1d0/dy_2*(V_(i,j+1,k) - V_(i,j,k)) - nu*1d0/dy_1*(V_(i,j,k)-V_(i,j-1,k)) ) + & ! d^2v/dy^2
                             1d0/dz_3*(nu*1d0/dz_2*(V_(i,j,k+1) - V_(i,j,k)) - nu*1d0/dz_1*(V_(i,j,k)-V_(i,j,k-1)) )     ! d^2v/dz^2

          End Do
       End Do
    End Do

    rhs_v = rhs_v + term_2(2:nxg-1,2:ny-1,2:nzg-1)

    !--------------Constant pressure gradient-------------!
    If (exterior_pressure_gradient == .True.) Then
      rhs_v = rhs_v + dPdy * Hv_interior(2:nxg-1,2:ny-1,2:nzg-1)
    Else
      rhs_v = rhs_v + dPdy
    End If

  End Subroutine compute_rhs_v

  !-------------------------------------------------------!
  !                Compute RHS for dw/dt                  !
  !                                                       !
  !  dw/dt = -duw/dx - dvw/dy - dw^2/dz + div(nu grad(w)) !
  !-------------------------------------------------------!
  Subroutine compute_rhs_w(U_,V_,W_,rhs_w)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg,nyg,nz), Intent(In) :: W_
    Real(Int64), Dimension(2:nxg-1,2:nyg-1,2:nz-1), Intent(Out) :: rhs_w

    ! local variables
    Integer(Int32) :: i, j, k    
    Real   (Int64) :: dx_1, dx_2, dx_3, maxerr, w0, w1, yd
    Real   (Int64) :: dy_1, dy_2, dy_3
    Real   (Int64) :: dz_1, dz_2, dz_3
    Real   (Int64) :: nu_x1, nu_x2, nu_y1, nu_y2, nu_z1, nu_z2 

    !-------------compute convective term_2s---------------!

    ! 1)---------compute -dw^2/dz--------------

    ! interpolate w in z (faces to centers)
    Call interpolate_z(W_,term_1(1:nxg,1:nyg,2:nz))

    ! compute w^2
    term_1 = term_1**2d0

    ! -dw^2/dz
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             term_2(i,j,k) = -(term_1(i,j,k+1)-term_1(i,j,k))/( zg(k+1) - zg(k) )
          End Do
       End Do
    End Do

    rhs_w = term_2(2:nxg-1,2:nyg-1,2:nz-1)

    ! 2)-----------compute -duw/dx--------------

    ! interpolate u in z (centers to faces)
    Call interpolate_z(U_,term_1(1:nx,1:nyg,1:nz)) 

    ! interpolate w in x (centers to faces)
    Call interpolate_x(W_,term_2(1:nx,1:nyg,1:nz)) 

    ! uw
    term_1 = term_1*term_2

    ! -duw/dx
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             term_2(i,j,k) = -( term_1(i,j,k) - term_1(i-1,j,k) )/( x(i) - x(i-1) )
          End Do
       End Do
    End Do

    rhs_w = rhs_w + term_2(2:nxg-1,2:nyg-1,2:nz-1)

    ! 3)--------------compute -dvw/dy--------------

    ! interpolate v in z (centers to faces)
    Call interpolate_z(V_,term_1(1:nxg,1:ny,1:nz)) 

    ! interpolate w in y (centers to faces)
    Call interpolate_y(W_,term_2(1:nxg,1:ny,1:nz)) 

    ! vw at faces
    term_1 = term_1*term_2

    ! -dvw/dy (this derivative goes to ym)
    Do k=2,nz-1
       Do j=2,nyg-1
          Do i=2,nxg-1
             term_2(i,j,k) = -( term_1(i,j,k) - term_1(i,j-1,k) )/( y(j) - y(j-1) )
          End Do
       End Do
    End Do

    rhs_w = rhs_w + term_2(2:nxg-1,2:nyg-1,2:nz-1)

    !--------------compute viscous term_2s------------!

    ! 4)----------compute div( nu grad(w) )-----------
    ! first derivation in y: dw/dy goes to yg_m(1:ny)
    Do i = 1, nxg
       Do k = 1, nz
          term_1(i,1:nyg-1,k) = ( W_(i,2:nyg,k) - W_(i,1:nyg-1,k) )/( yg(2:nyg) - yg(1:nyg-1) ) 
       End Do
    End Do

    Do k=2,nz-1
       
       ! first derivation in z (W at faces)
       dz_1 = z(  k) - z(k-1) 
       dz_2 = z(k+1) - z(k  )
       ! second derivation in z (W at centers)
       dz_3 = zg(k+1) - zg(k)
       
       Do j=2,nyg-1

          ! second derivation in y (W at faces)
          dy_3 = y(j) - y(j-1)

          Do i=2,nxg-1
             
             ! first derivation in x (W at centers)
             dx_1 = xg(  i) - xg(i-1) 
             dx_2 = xg(i+1) - xg(i  )
             ! second derivation in x (W at faces)
             dx_3 = x(i) - x(i-1)   

             ! viscous term_2
             term_2(i,j,k) = 1d0/dx_3*(nu*1d0/dx_2*(W_(i+1,j,k) - W_(i,j,k)) - nu*1d0/dx_1*(W_(i,j,k) - W_(i-1,j,k)) ) + & ! d^2w/dx^2
                           1d0/dy_3*(nu*term_1(i,j,k)                      - nu*term_1(i,j-1,k) )                    + & ! d^2w/dy^2
                           1d0/dz_3*(nu*1d0/dz_2*(W_(i,j,k+1) - W_(i,j,k)) - nu*1d0/dz_1*(W_(i,j,k) - W_(i,j,k-1)) )     ! d^2w/dz^2

          End Do
       End Do
    End Do

    rhs_w = rhs_w + term_2(2:nxg-1,2:nyg-1,2:nz-1)

    !--------------Constant pressure gradient-------------!
    If (exterior_pressure_gradient == .True.) Then
      rhs_w = rhs_w + dPdz * Hw_interior(2:nxg-1,2:nyg-1,2:nz-1)
    Else
      rhs_w = rhs_w + dPdz
    End If

  End Subroutine compute_rhs_w

End Module equations

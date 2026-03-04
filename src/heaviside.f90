!------------------------------------------!
!      Module for Heaviside fields         !
!------------------------------------------!
Module heaviside

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use immersed_boundary_operators
  Use boundary_conditions
  Use poisson
  Use interpolation
  Use differential_operators

  ! prevent implicit typing
  Implicit None

Contains

  Subroutine compute_heaviside

    Integer(Int32) :: i, j, k
    Hc_exterior = 0d0
    Hc_interior = 0d0
    debug_rhs_p=0d0

    Fibu = 0d0
    Fibv = 0d0
    Fibw = 0d0

    Call regu(Fibu, normals)
    Call regv(Fibv, normals)
    Call regw(Fibw, normals)
    Call apply_boundary_conditions(Fibu, Fibv, Fibw)

    !Call divergence(Hc_interior(2:nxg-1, 2:nyg-1, 2:nzg-1), Fibu, Fibv, Fibw)
    Do k = 2, nzg-1
      Do j = 2, nyg-1
         Do i = 2, nxg-1
          Hc_interior(i,j,k) = ( Fibu(i,j,k) - Fibu(i-1,j,k) ) / ( x(i)-x(i-1) ) + & ! d u/dx
                          ( Fibv(i,j,k) - Fibv(i,j-1,k) ) / ( y(j)-y(j-1) ) + & ! d v/dy
                          ( Fibw(i,j,k) - Fibw(i,j,k-1) ) / ( z(k)-z(k-1) )     ! d w/dz
         End Do
      End Do
    End Do
    Hc_interior = -Hc_interior
    !debug_rhs_p=Hc_interior
    Call apply_periodic_xz_centers(Hc_interior)
    Call update_ghost_interior_planes_centers(Hc_interior)
    !debug_rhs_p=Hc_interior
    !debug_rhs_p=Hc_interior
    Call solve_poisson_equation(Hc_interior)
    !debug_rhs_p=Hc_interior
    Hc_exterior = 1 - Hc_interior

    Call interpolate_x(Hc_exterior, Hu_exterior(2:nx-1 , 2:nyg-1, 2:nzg)) 
    Call interpolate_y(Hc_exterior, Hv_exterior(2:nxg-1, 2:ny-1 , 2:nzg)) 
    Call interpolate_z(Hc_exterior, Hw_exterior(2:nxg-1, 2:nyg-1, 2:nz))
    Call apply_boundary_conditions(Hu_exterior, Hv_exterior, Hw_exterior)

    Hu_interior = 1 - Hu_exterior
    Hv_interior = 1 - Hv_exterior
    Hw_interior = 1 - Hw_exterior

    Call apply_boundary_conditions(Hu_interior, Hv_interior, Hw_interior)
    
    E1nHc_exterior = regTc_1n(Hc_exterior)
    E1nH_exterior = regT_1n(Hu_exterior, Hv_exterior, Hw_exterior)

    debug_surface_scalar = E1nHc_exterior
    !debug_rhs_p=Hc_exterior
    ! debug_u=Fibu
    ! debug_v=Fibv
    ! debug_w=Fibw

  End Subroutine compute_heaviside

End module heaviside

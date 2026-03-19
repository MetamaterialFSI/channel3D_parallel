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

    U_reg = 0d0
    V_reg = 0d0
    W_reg = 0d0

    Call regu(U_reg, normals(1 : nb))             
    Call regv(V_reg, normals(nb + 1 : 2 * nb))    
    Call regw(W_reg, normals(2 * nb + 1 : 3 * nb))
    Call apply_boundary_conditions(U_reg, V_reg, W_reg)

    Call divergence(Hc_interior, U_reg, V_reg, W_reg)
    Hc_interior = -Hc_interior
    Call solve_poisson_equation(Hc_interior)
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

  End Subroutine compute_heaviside

End module heaviside

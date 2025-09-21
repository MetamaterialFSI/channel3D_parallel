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
    Real(Int64), Dimension(3 * nb) :: f_
    Integer(Int32) :: k

    Hc_exterior = 0d0
    Hc_interior = 0d0

    ! test for v
    ! Call regu(Fibu, normals)
    ! Call regv(Fibv, normals)
    ! Call regw(Fibw, normals)
    ! test for w
    
    Call regu(Fibu, tangents_2)
    Call regv(Fibv, tangents_2)
    Call regw(Fibw, tangents_2)
    Call apply_boundary_conditions(Fibu, Fibv, Fibw)

    ! f_ = 0d0
    ! Do k=1,nb
    !   f_(nb + k) = sin(2*pi*zb(k)/lzp)
    ! End Do

    ! Hc_exterior = 0d0
    ! Hc_interior = 0d0

    ! Call regu(Fibu, f_)
    ! Call regv(Fibv, f_)
    ! Call regw(Fibw, f_)
    ! Call apply_boundary_conditions(Fibu, Fibv, Fibw)

    ! temporarily for debugging:
    Hu_exterior = Fibu
    Hv_exterior = Fibv
    Hw_exterior = Fibw

    ! Call divergence(Hc_interior, Fibu, Fibv, Fibw)
    ! Hc_interior =  -Hc_interior
    ! Call solve_poisson_equation(Hc_interior)
    ! Hc_exterior = 1 - Hc_interior

    ! Call interpolate_x(Hc_exterior, Hu_exterior(2:nx-1 , 2:nyg-1, 2:nzg-1)) 
    ! Call interpolate_y(Hc_exterior, Hv_exterior(2:nxg-1, 2:ny-1 , 2:nzg-1)) 
    ! Call interpolate_z(Hc_exterior, Hw_exterior(2:nxg-1, 2:nyg-1, 2:nz-1))

    ! Call apply_boundary_conditions(Hu_exterior, Hv_exterior, Hw_exterior)

    ! Hu_interior = 1 - Hu_exterior
    ! Hv_interior = 1 - Hv_exterior
    ! Hw_interior = 1 - Hw_exterior

  End Subroutine compute_heaviside

End module heaviside

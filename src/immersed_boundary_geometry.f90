Module immersed_boundary_geometry

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global

  ! prevent implicit typing
  Implicit None

  ! declarations
Contains

  Subroutine compute_nb
    
    nb = nxb * nzb
    dxb = real(Lxp / nxb, 8)
    dzb = real(Lzp / nzb, 8)

  End Subroutine compute_nb

  Subroutine setup_IB_geometry
    Integer(Int32) :: i, j, k, l

    Select Case (body_type)
      Case (0) ! No IB

      Case (1) ! Static wall
        ! Scalar arrays
        Do j=1,nzb
           Do i=1,nxb
              xb(i + (j-1) * nxb) = (real(i,8) - 0.75d0) * dxb
              yb(i + (j-1) * nxb) = 1.0d0
              zb(i + (j-1) * nxb) = (real(j,8) - 0.75d0) * dzb
           End Do
        End Do
        sb = dxb * dzb 

      Case (2) ! Standing wave motion

      Case (3) ! Travelling wave motion

    End Select

  End Subroutine setup_IB_geometry

End Module immersed_boundary_geometry

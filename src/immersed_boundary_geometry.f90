Module immersed_boundary_geometry

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global

  ! prevent implicit typing
  Implicit None

  ! declarations
Contains

  Subroutine compute_nb

    If (body_type == 0) Then
      nxb = 0
      nzb = 0
    End If

    nb = nxb * nzb
    dxb = real(Lxp / nxb, 8)
    dzb = real(Lzp / nzb, 8)

  End Subroutine compute_nb

  Subroutine setup_IB_geometry
    Integer(Int32) :: i, j, k, l, idx

    Select Case (body_type)
      Case (0) ! No IB
        nb_start = 0
        nb_end = -1
      Case (1) ! Static wall
        moving_body = .False.
        ub = 0d0
        ! Scalar arrays. Arrange such that the points treated by one partition are contiguous
        ! (i.e., fall between an nb_start and nb_end)
        nb_start = nb + 1  ! Initialize to an invalid value (beyond the max index)
        nb_end = 0         ! Initialize to the lowest possible index
        Do j=1,nzb
          Do i=1,nxb
            idx = i + (j-1) * nxb
            xb(idx) = (real(i,8) - 0.75d0) * dxb
            yb(idx) = 1.0d0
            zb(idx) = (real(j,8) - 0.75d0) * dzb
            If (zb(idx) >= z(1) .and. nb_start > idx) then
              nb_start = idx
            End If
            If (zb(idx) <= z(nz-1) .and. nb_end < idx) then
              nb_end = idx
            End If  
          End Do
        End Do
        sb = dxb * dzb 

      Case (2) ! Standing wave motion

      Case (3) ! Travelling wave motion

    End Select

  End Subroutine setup_IB_geometry

End Module immersed_boundary_geometry

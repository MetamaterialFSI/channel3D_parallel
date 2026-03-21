!------------------------------------------------!
!      Module for differential operators         !
!------------------------------------------------!
Module differential_operators

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global, Only : x, y, z, xg, yg, zg, nx, ny, nz, nxg, nyg, nzg
  !Use initialization

  ! prevent implicit typing
  Implicit None

Contains

  Subroutine gradient(Gpx_, Gpy_, Gpz_, p_)

    Integer(Int32) :: i, j, k
    Real   (Int64), Contiguous, Intent(In)  :: p_(2:, 2:, 2:)
    Real   (Int64), Contiguous, Intent(Out) :: Gpx_(:, :, :)
    Real   (Int64), Contiguous, Intent(Out) :: Gpy_(:, :, :)
    Real   (Int64), Contiguous, Intent(Out) :: Gpz_(:, :, :)

    Gpx_ = 0d0
    Gpy_ = 0d0
    Gpz_ = 0d0

    Do k = 2, nzg - 1
      Do j = 2, nyg - 1
        Do i = 2, nx - 1
          Gpx_(i, j, k) = (p_(i+1, j, k) - p_(i, j, k)) / (xg(i + 1) - xg(i))
        End Do
      End Do
    End Do

    Do k = 2, nzg - 1
      Do j = 2, ny - 1
        Do i = 2, nx - 1
          Gpy_(i, j, k) = (p_(i, j + 1, k) - p_(i, j, k)) / (yg(j + 1) - yg(j))
        End Do
      End Do
    End Do

    Do k = 2, nz - 1
      Do j = 2, nyg - 1
        Do i = 2, nx - 1
          Gpz_(i, j, k) = (p_(i, j, k + 1) - p_(i, j, k)) / (zg(k + 1) - zg(k))
        End Do
      End Do
    End Do

  End Subroutine gradient

  Subroutine divergence(div_, U_, V_, W_)

    Integer(Int32) :: i, j, k
    Real   (Int64), Contiguous, Intent(Out) :: div_(2:, 2:, 2:)
    Real   (Int64), Contiguous, Intent(In)  :: U_(:, :, :)
    Real   (Int64), Contiguous, Intent(In)  :: V_(:, :, :)
    Real   (Int64), Contiguous, Intent(In)  :: W_(:, :, :)

    div_ = 0d0

    Do k = 2, nzg - 1
      Do j = 2, nyg - 1
        Do i = 2, nxg - 1
          div_(i,j,k) = (U_(i, j, k) - U_(i-1, j,   k))   / (x(i) - x(i - 1)) + & ! du/dx
                        (V_(i, j, k) - V_(i,   j-1, k))   / (y(j) - y(j - 1)) + & ! dv/dy
                        (W_(i, j, k) - W_(i,   j,   k-1)) / (z(k) - z(k - 1))     ! dw/dz
        End Do
      End Do
    End Do

  End Subroutine divergence

End Module differential_operators

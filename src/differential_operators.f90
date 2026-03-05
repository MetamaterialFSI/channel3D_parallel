!------------------------------------------------!
!      Module for differential operators         !
!------------------------------------------------!
Module differential_operators

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global, Only : x, y, z, xg, yg, zg, nx, ny, nz, nxg, nyg, nzg, nxm, nym, nzm
  !Use initialization

  ! prevent implicit typing
  Implicit None

Contains

  Subroutine gradient(Gpx_, Gpy_, Gpz_, p_)

    Integer(Int32) :: i, j, k
    Real   (Int64), DIMENSION( 2:nxg-1, 2:nyg-1, 2:nzg   ), INTENT(IN)  :: p_
    Real   (Int64), DIMENSION(      nx,   nym+2,   nzm+2 ), INTENT(OUT) :: Gpx_
    Real   (Int64), DIMENSION(   nxm+2,      ny,   nzm+2 ), INTENT(OUT) :: Gpy_
    Real   (Int64), DIMENSION(   nxm+2,   nym+2,      nz ), INTENT(OUT) :: Gpz_

    Gpx_ = 0d0
    Gpy_ = 0d0
    Gpz_ = 0d0

    Do i=2,nx-1
       Gpx_(i,2:nyg-1,2:nzg-1) = ( p_(i+1,2:nyg-1,2:nzg-1) - p_(i,2:nyg-1,2:nzg-1) ) / ( xg(i+1) - xg(i) )
    End Do

    Do j=2,ny-1
       Gpy_(2:nxg-1,j,2:nzg-1) = ( p_(2:nxg-1,j+1,2:nzg-1) - p_(2:nxg-1,j,2:nzg-1) ) / ( yg(j+1) - yg(j) )
    End Do

    Do k=2,nz-1
       Gpz_(2:nxg-1,2:nyg-1,k) = ( p_(2:nxg-1,2:nyg-1,k+1) - p_(2:nxg-1,2:nyg-1,k) ) / ( zg(k+1) - zg(k) )
    End Do

  End Subroutine gradient

  Subroutine velocity_gradient_u(Gux_, Guy_, Guz_, U_)

    Integer(Int32) :: i, j, k
    Real   (Int64), DIMENSION(   nx,      nyg,      nzg   ), INTENT(IN)  :: U_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Gux_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Guy_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Guz_

    Gux_ = 0d0
    Guy_ = 0d0
    Guz_ = 0d0

    Do i=1,nx-1
       Gux_(i+1,1:nyg,1:nzg) = ( U_(i+1,1:nyg,1:nzg) - U_(i,1:nyg,1:nzg) ) / ( x(i+1) - x(i) )
    End Do

    Do j=1,nyg-1
       Guy_(1:nx,j,1:nzg)    = ( U_(1:nx,j+1,1:nzg)  - U_(1:nx,j,1:nzg) ) / ( yg(j+1) - yg(j) )
    End Do

    Do k=1,nzg-1
       Guz_(1:nx,1:nyg,k)    = ( U_(1:nx,1:nyg,k+1)  - U_(1:nx,1:nyg,k) ) / ( zg(k+1) - zg(k) )
    End Do

  End Subroutine velocity_gradient_u

  Subroutine velocity_gradient_v(Gvx_, Gvy_, Gvz_, V_)

    Integer(Int32) :: i, j, k
    Real   (Int64), DIMENSION(   nxg,      ny,     nzg   ), INTENT(IN)  :: V_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Gvx_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Gvy_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Gvz_

    Gvx_ = 0d0
    Gvy_ = 0d0
    Gvz_ = 0d0

    Do i=1,nxg-1
       Gvx_(i,1:ny,1:nzg)    = ( V_(i+1,1:ny,1:nzg) - V_(i,1:ny,1:nzg) ) / ( xg(i+1) - xg(i) )
    End Do

    Do j=1,ny-1
       Gvy_(1:nxg,j+1,1:nzg) = ( V_(1:nxg,j+1,1:nzg) - V_(1:nxg,j,1:nzg) ) / ( y(j+1) - y(j) )
    End Do

    Do k=1,nzg-1
       Gvz_(1:nxg,1:ny,k)    = ( V_(1:nxg,1:ny,k+1) - V_(1:nxg,1:ny,k) ) / ( zg(k+1) - zg(k) )
    End Do

  End Subroutine velocity_gradient_v

  Subroutine velocity_gradient_w(Gwx_, Gwy_, Gwz_, W_)
    Integer(Int32) :: i, j, k
    Real   (Int64), DIMENSION(   nxg,      nyg,     nz   ), INTENT(IN)  :: W_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Gwx_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Gwy_
    Real   (Int64), DIMENSION(   nxg, nyg, nzg   ), INTENT(OUT) :: Gwz_

    Gwx_ = 0d0
    Gwy_ = 0d0
    Gwz_ = 0d0

    Do i=1,nxg-1
       Gwx_(i,1:nyg,1:nz)    = ( W_(i+1,1:nyg,1:nz) - W_(i,1:nyg,1:nz) ) / ( xg(i+1) - xg(i) )
    End Do

    Do j=1,nyg-1
       Gwy_(1:nxg,j,1:nz)    = ( W_(1:nxg,j+1,1:nz) - W_(1:nxg,j,1:nz) ) / ( yg(j+1) - yg(j) )
    End Do

    Do k=1,nz-1
       Gwz_(1:nxg,1:nyg,k+1) = ( W_(1:nxg,1:nyg,k+1) - W_(1:nxg,1:nyg,k) ) / ( z(k+1) - z(k) )
    End Do

  End Subroutine velocity_gradient_w

  Subroutine divergence(div_, U_, V_, W_)

    Integer(Int32) :: i, j, k
    Real   (Int64), DIMENSION( 2:nxg-1, 2:nyg-1, 2:nzg   ), INTENT(OUT) :: div_
    Real   (Int64), DIMENSION(      nx,   nym+2,   nzm+2 ), INTENT(IN)  :: U_
    Real   (Int64), DIMENSION(   nxm+2,      ny,   nzm+2 ), INTENT(IN)  :: V_
    Real   (Int64), DIMENSION(   nxm+2,   nym+2,      nz ), INTENT(IN)  :: W_

    div_ = 0d0

    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             div_(i,j,k) = ( U_(i,j,k) - U_(i-1,j,k) ) / ( x(i)-x(i-1) ) + & ! d u/dx
                           ( V_(i,j,k) - V_(i,j-1,k) ) / ( y(j)-y(j-1) ) + & ! d v/dy
                           ( W_(i,j,k) - W_(i,j,k-1) ) / ( z(k)-z(k-1) )     ! d w/dz
          End Do
       End Do
    End Do

  End Subroutine divergence

End Module differential_operators

!------------------------------------------!
!        Module for interpolation          !
!------------------------------------------!
Module interpolation

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64

  ! prevent implicit typing
  Implicit None

Contains

  !-------------------------------------------------------!
  !            Linear interpolation of u in x             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  !-------------------------------------------------------!
  ! uniform mesh assumed
  Subroutine interpolate_x(u,ui)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)

    Integer(Int32) :: n(3), n1, n2, n3

    n  = Shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ui = 0d0
    ui( 1:n1-1, 1:n2, 1:n3) = 0.5d0*( u(1:n1-1,:,:) + u(2:n1,:,:) )

  End Subroutine interpolate_x

  !-------------------------------------------------------!
  !            Linear interpolation of u in y             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  !-------------------------------------------------------!
  Subroutine interpolate_y(u,ui)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)

    Integer (Int32) :: n(3), n1, n2, n3, i1, i3

    n  = Shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ui = 0d0
    ui(1:n1, 1:n2-1, 1:n3) = 0.5d0*( u(:,1:n2-1,:) + u(:,2:n2,:) ) 
   
  End Subroutine interpolate_y

  !-------------------------------------------------------!
  !            Linear interpolation of u in z             !
  !                                                       !
  ! Input : u                                             !
  ! Output: ui                                            !
  !-------------------------------------------------------!
  ! uniform mesh assumed
  Subroutine interpolate_z(u,ui)

    Real    (Int64), Intent(In)  :: u(:,:,:)
    Real    (Int64), Intent(Out) :: ui(:,:,:)

    Integer(Int32) :: n(3), n1, n2, n3

    n  = Shape(u)
    n1 = n(1)
    n2 = n(2)
    n3 = n(3)

    ui = 0d0
    ui(1:n1, 1:n2, 1:n3-1) = 0.5d0*( u(:,:,1:n3-1) + u(:,:,2:n3) )
    
  End Subroutine interpolate_z

End Module interpolation

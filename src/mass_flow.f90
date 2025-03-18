!----------------------------------------!
!    Module for constant mass flow       !
!----------------------------------------!
Module mass_flow

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global,          Only : nx, nxg, nyg, ny, nzg, yg, y, & 
                              Qflow_x_0, Qflow_y_0, ierr
  Use mpi 

  ! prevent implicit typing
  Implicit None

Contains

  !-----------------------------------------------!
  !             Compute mean mass flow            !
  !            for U between yg(2:nyg-1)          !
  !              with trapezoidal rule            !
  ! Input:   U                                    !
  ! Output:  Qflow                                !
  !-----------------------------------------------!
  Subroutine compute_mean_mass_flow_U(U,Qflow)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U    
    Real(Int64), Intent(Out) :: Qflow

    Real   (Int64) :: Qflow_local, norm_local
    Real   (Int64) :: norm
    Integer(Int32) :: i, j, k

    ! compute local mass flow
    Qflow_local = 0d0
    norm_local  = 0d0
    Do k=2,nzg-1
       Do j=3,nyg-1
          Do i=2,nx-1
             Qflow_local = Qflow_local + ( U(i,j,k) + U(i,j-1,k) )*0.5d0*( yg(j)-yg(j-1) )
             norm_local  = norm_local  + ( 1d0 + 1d0 )*0.5d0*( yg(j)-yg(j-1) )
          End Do
       End Do
    End Do

    ! compute total mass flow
    Call MPI_AllReduce(Qflow_local,Qflow,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_AllReduce(norm_local,  norm,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    
    Qflow = Qflow/norm

  End Subroutine compute_mean_mass_flow_U

  !-----------------------------------------------!
  !             Compute mean mass flow            !
  !              for V between y(1:ny)            !
  !              with trapezoidal rule            !
  ! Input:   V                                    !
  ! Output:  Qflow                                !
  !-----------------------------------------------!
  Subroutine compute_mean_mass_flow_V(V,Qflow)

    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V
    Real(Int64), Intent(Out) :: Qflow

    Real   (Int64) :: Qflow_local, norm_local
    Real   (Int64) :: norm
    Integer(Int32) :: i, j, k

    ! compute local mass flow
    Qflow_local = 0d0
    norm_local  = 0d0
    Do k=2,nzg-1
       Do j=2,ny
          Do i=2,nxg-1
             Qflow_local = Qflow_local + ( V(i,j,k) + V(i,j-1,k) )*0.5d0*( y(j)-y(j-1) )
             norm_local  = norm_local  + ( 1d0 + 1d0 )*0.5d0*( y(j)-y(j-1) )
          End Do
       End Do
    End Do

    ! compute total mass flow
    Call MPI_AllReduce(Qflow_local,Qflow,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_AllReduce(norm_local,  norm,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    
    Qflow = Qflow/norm

  End Subroutine compute_mean_mass_flow_V
  
  !------------------------------------------------------------!
  !           Compute dPx for constant mass flow in x          !
  !        The mass flow is conserved in U at ym points        !
  !                                                            !
  !     dPdx      = Qflow_0 - trapz( U(tn+1) )/Volume          !
  !     Qflow_x_0 = initial mass flow                          !
  !     trapz     = trapezoidal rule                           !
  !     ym        = yg(2:nyg-1)                                !
  !     y-weights for trapz = 0.5d0*(yg(2:nyg-2)-yg(3:nyg-1))  !
  !                                                            !
  ! Input:  U, Qflow_x_0                                       !
  ! Output: dPdx                                               !
  !------------------------------------------------------------!
  Subroutine compute_dPx_for_constant_mass_flow(U,dPdx)

    Real(Int64), Dimension(nx,nyg,nzg), Intent(In) :: U
    Real(Int64), Intent(Out) :: dPdx

    ! local variables
    Real   (Int64) :: norm_local, Int_U_local
    Real   (Int64) :: norm,       Int_U
    Integer(Int32) :: i, k, j

    ! compute integrals with trapezoidal rule
    Int_U_local = 0d0
    norm_local  = 0d0
    Do k=2,nzg-1
       Do j=3,nyg-1
          Do i=2,nx-1
             Int_U_local = Int_U_local + ( U (i,j,k) + U (i,j-1,k) )*0.5d0*( yg(j)-yg(j-1) )
             norm_local  = norm_local  + ( 1d0 + 1d0 )*0.5d0*( yg(j)-yg(j-1) )
          End Do
       End Do
    End Do

    ! compute total mass flow
    Call MPI_AllReduce(Int_U_local,Int_U,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_AllReduce(norm_local,  norm,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)

    ! compute pressure gradient
    dPdx = Qflow_x_0 - Int_U/norm

  End Subroutine compute_dPx_for_constant_mass_flow

  !------------------------------------------------------------!
  !           Compute dPy for constant mass flow in y          !
  !        The mass flow is conserved in V at y points         !
  !                                                            !
  !     dPdy      = Qflow_y_0 - trapz( V(tn+1) )/Volume        !
  !     Qflow_y_0 = initial mass flow                          !
  !     trapz     = trapezoidal rule                           !
  !     y         = y(1:ny)                                    !
  !     y-weights for trapz = 0.5d0*(y(1:ny-1)-y(2:ny))        !
  !                                                            !
  ! Input:  V, Qflow_y_0                                       !
  ! Output: dPdy                                               !
  !------------------------------------------------------------!
  Subroutine compute_dPy_for_constant_mass_flow(V,dPdy)

    Real(Int64), Dimension(nxg,ny,nzg), Intent(In) :: V
    Real(Int64), Intent(Out) :: dPdy

    ! local variables
    Real   (Int64) :: norm_local, Int_U_local
    Real   (Int64) :: norm,       Int_U
    Integer(Int32) :: i, k, j

    ! compute integrals with trapezoidal rule
    Int_U_local = 0d0
    norm_local  = 0d0
    Do k=2,nzg-1
       Do j=2,ny
          Do i=2,nxg-1
             Int_U_local = Int_U_local + ( V (i,j,k) + V (i,j-1,k) )*0.5d0*( y(j)-y(j-1) )
             norm_local  = norm_local  + ( 1d0 + 1d0 )*0.5d0*( y(j)-y(j-1) )
          End Do
       End Do
    End Do

    ! compute total mass flow
    Call MPI_AllReduce(Int_U_local,Int_U,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)
    Call MPI_AllReduce(norm_local,  norm,1,MPI_real8,MPI_sum,MPI_COMM_WORLD,ierr)

    ! compute pressure gradient
    dPdy = Qflow_y_0 - Int_U/norm

  End Subroutine compute_dPy_for_constant_mass_flow

End Module mass_flow

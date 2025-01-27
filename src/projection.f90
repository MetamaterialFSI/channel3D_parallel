!------------------------------------------------!
!      Module for fractional step method         !
!------------------------------------------------!
Module projection

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi

  ! prevent implicit typing
  Implicit None

Contains

  !--------------------------------------------!
  !    Compute incompressible velocity with    !
  !         fractional step method             !
  !--------------------------------------------!
  Subroutine compute_projection_step

    Call compute_pseudo_pressure_rhs

    Call solve_poisson_equation

    Call project_velocity

  End Subroutine compute_projection_step

  !------------------------------------------------------!
  ! Compute right-hand size for pseudo-pressure equation !
  !                                                      !
  ! Equation:                                            !
  !   Laplacian(p) = rhs_p = div(vel*)                   !
  !                                                      !
  ! Input:  U, V, W                                      !
  ! Output: rhs_p                                        !
  !                                                      !
  !------------------------------------------------------!
  Subroutine compute_pseudo_pressure_rhs

    Integer(Int32) :: i, j, k

    ! rhs_p locate at cell centers
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             rhs_p(i,j,k)  = ( U(i,j,k) - U(i-1,j,k) )/( x(i)-x(i-1) ) + & ! du*/dx 
                             ( V(i,j,k) - V(i,j-1,k) )/( y(j)-y(j-1) ) + & ! dv*/dy
                             ( W(i,j,k) - W(i,j,k-1) )/( z(k)-z(k-1) )     ! dw*/dz
          End Do
       End Do
    End Do

  End Subroutine compute_pseudo_pressure_rhs

  !----------------------------------------------------!
  !             Solve pseudo-pressure equation         !
  !                                                    !
  !   Equation:                                        !
  !     Laplacian(p) = rhs_p = div(vel*)               !
  !                                                    !
  !   Boundary conditions:                             !
  !     dP/dy = 0 at bottom wall                       !
  !     dP/dy = 0 at top    wall                       !
  !                                                    !
  !         This is particular for channels:           !
  !               Fast Poisson solver                  !
  !                                                    !
  ! Input:  rhs_p (rhs for the pseudo-pressure)        !
  ! Output: rhs_p (pseudo-pressure)                    !
  !                                                    !
  !----------------------------------------------------!
  Subroutine solve_poisson_equation

    Integer(Int32) :: i, j, k, k_global, i_global, info
    Real   (Int64) :: dum, dumref, maxerr

    ! Lapack function for solving tridiagonal systems
    External :: zgesv

    ! 2D Fourier transform interior points 
    Do j = 2, nyg-1
       plane = dcmplx( rhs_p ( 2:nxp+1, j, 2:nzp+1 ) )
       Call fftw_mpi_execute_dft(plan_d, plane, plane_hat)
       rhs_p_hat (:, j, :) = plane_hat
    End Do

    ! solve for each mode
    Do k = 0, mz
       Do i = 0, mx
          ! mapping to x-mode 
          i_global = imode_map_fft( i, k )
          ! mapping to z-mode 
          k_global = kmode_map_fft( i, k )
          ! form matrix          
          Do j = 2, nyg-1 ! diagonal
             D(j) = Dyy(j,j) + kxx(i_global) + kzz(k_global) 
          End Do
          Do j = 2, nyg-2 ! lower diagonal
             DL(j) = Dyy(j+1,j)
          End Do
          Do j = 2, nyg-2 ! upper diagonal
             DU(j) = Dyy(j,j+1)
          End Do
          ! remove singularity 00 mode (set a reference pseudo-pressure)
          If ( i_global==0 .And. k_global==0 ) D(2) = 3d0/2d0*D(2) 
          ! solve M*u = rhs (solution stored in rhs_p_hat)
          rhs_aux = rhs_p_hat (i, :, k)
          Call Zgtsv( nr, nrhs, DL, D, DU, rhs_aux, nr, info) 
          rhs_p_hat (i, :, k) = rhs_aux
       End Do
    End Do

    ! 2D inverse Fourier transform 
    Do j = 2, nyg-1
       plane_hat = rhs_p_hat (:, j, :)
       Call fftw_mpi_execute_dft(plan_i, plane_hat, plane)
       rhs_p ( 2:nxp+1, j, 2:nzp+1 ) = plane/Real( nxp_global*nzp_global, 8)
    End Do

    ! periodic boundary conditions in x and z
    Call apply_periodic_xz_pressure

    ! update ghost interior planes
    Call update_ghost_interior_planes_pressure

    ! save pressure for statistics
    If ( rk_step == 1 ) Then
       P( 2:nxg-1, 2:nyg-1, 2:nzg-1 ) = rhs_p/(dt*rk_coef(1,1)) ! to be checked
       P(  1,:,:) = P(    2,:,:)
       P(nxg,:,:) = P(nxg-1,:,:)
       P(:,  1,:) = P(:,    2,:)
       P(:,nyg,:) = P(:,nyg-1,:)
       P(:,:,  1) = P(:,:,    2)
       P(:,:,nzg) = P(:,:,nzg-1)
    End If

  End Subroutine solve_poisson_equation

  !--------------------------------------------------!
  ! Periodic boundary conditions for pseudo-pressure !
  !--------------------------------------------------!
  Subroutine apply_periodic_xz_pressure

    ! apply periodicity in x (All processors, no MPI needed)
    rhs_p ( nxg-1, :, : ) = rhs_p ( 2, :, : )

    ! apply periodicity in z (Only first and last processor, MPI needed) 
    If     ( myid==0 ) Then
       buffer_p = rhs_p ( 2:nxg-1, :, 2 ) 
       ! send data to nprocs-1
       Call Mpi_send(buffer_p, (nxg-2)*(nyg-2), MPI_real8, nprocs-1, 0, &
             MPI_COMM_WORLD,ierr)
    Elseif ( myid==nprocs-1 ) Then
       ! receive data from 0
       Call Mpi_recv(buffer_p, (nxg-2)*(nyg-2), MPI_real8, 0, 0, &
            MPI_COMM_WORLD,istat,ierr)
       rhs_p ( 2:nxg-1, :, nzp+1+1 ) = buffer_p
    End If

  End Subroutine apply_periodic_xz_pressure

  !------------------------------------------------!
  !    Update ghost interior planes for pressure   !
  !------------------------------------------------!
  Subroutine update_ghost_interior_planes_pressure

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom
    
    !----------------------update P-----------------------!    
    ! send to bottom processor, receive from top one
    sendto   = myid - 1
    tagto    = myid - 1
    recvfrom = myid + 1
    tagfrom  = myid 
    If ( myid==0 ) Then
       sendto = MPI_PROC_NULL
       tagto  = 0
    End If
    If ( myid==nprocs-1 ) Then
       recvfrom = MPI_PROC_NULL
       tagfrom  = MPI_ANY_TAG
    End If
    buffer_ps = rhs_p(:,:,2)  ! send buffer
    Call Mpi_sendrecv(buffer_ps, (nxg-2)*(nyg-2), Mpi_real8, sendto, tagto,        &
         buffer_pr, (nxg-2)*(nyg-2), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
         istat, ierr)   
    If ( myid/=nprocs-1 ) rhs_p(:,:,nzg) = buffer_pr ! received buffer 

  End Subroutine update_ghost_interior_planes_pressure

  !------------------------------------------------!
  !             Project velocity field             !
  !                                                !
  !             vel = vel* - gradient(p)           !
  !                                                !
  ! Input:  U*,V*,W*,P* (-> U,V,W,rhs_p)           !
  ! Output: U, V, W     (projected)                !
  !                                                !
  !------------------------------------------------!
  Subroutine project_velocity

    Integer(Int32) :: i, j, k

    ! project U (interior points)
    Do i = 2, nx-1
       U(i,2:nyg-1,2:nzg-1) = U(i,2:nyg-1,2:nzg-1) - & 
            ( rhs_p(i+1,2:nyg-1,2:nzg-1) - rhs_p(i,2:nyg-1,2:nzg-1) )/( xg(i+1) - xg(i) )
    End Do

    ! project V (interior points)
    Do j = 2, ny-1
       V(2:nxg-1,j,2:nzg-1) = V(2:nxg-1,j,2:nzg-1) - & 
            ( rhs_p(2:nxg-1,j+1,2:nzg-1) - rhs_p(2:nxg-1,j,2:nzg-1) )/( yg(j+1) - yg(j) )
    End Do

    ! project W (interior points)
    Do k = 2, nz-1
       W(2:nxg-1,2:nyg-1,k) = W(2:nxg-1,2:nyg-1,k) - & 
            ( rhs_p(2:nxg-1,2:nyg-1,k+1) - rhs_p(2:nxg-1,2:nyg-1,k) )/( zg(k+1) - zg(k) )
    End Do

  End Subroutine project_velocity

  !----------------------------------------------------!
  !   Compute maximum divergence for interior points   !
  !                                                    !
  ! Input:  U, V, W                                    !
  ! Output: div, max_divergence                        !
  !                                                    !
  !----------------------------------------------------!
  Subroutine check_divergence(max_divergence)

    Real(Int64), Intent(Out) :: max_divergence

    Real   (Int64) :: max_divergence_local, div
    Integer(Int32) :: i, j, k    

    div = 0d0
    Do k = 2, nzg-1     
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             div  = Max(div, ( U(i,j,k) - U(i-1,j,k) )/( x(i)-x(i-1) ) + &  ! du/dx 
                             ( V(i,j,k) - V(i,j-1,k) )/( y(j)-y(j-1) ) + &  ! dv/dy
                             ( W(i,j,k) - W(i,j,k-1) )/( z(k)-z(k-1) )     )! dw/dz
          End Do
       End Do
    End Do

    max_divergence_local = div

    Call MPI_Reduce(max_divergence_local,max_divergence,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

  End Subroutine check_divergence

End module projection

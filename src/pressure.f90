!------------------------------------------------!
!      Module for computing actual pressure      !
!------------------------------------------------!
Module pressure

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use interpolation
  Use equations
  Use boundary_conditions

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------------!
  !              Compute actual pressure               !
  !----------------------------------------------------! 
  Subroutine compute_pressure

    Call compute_pressure_BC
  
    Call compute_pressure_rhs

    Call solve_poisson_equation_for_pressure

  End Subroutine compute_pressure

  !-------------------------------------------------------!
  !     Compute right-hand size for pressure equation     !
  !                                                       !
  !   Equation:                                           !
  !    dvel/dt      = - grad(P) + rhs_vel                 !
  !    Laplacian(p) = rhs_p = div(rhs_vel)                !
  !                                                       !
  ! Input:  U, V, W                                       !
  ! Output: rhs_p                                         !
  !                                                       !
  !-------------------------------------------------------!
  Subroutine compute_pressure_rhs

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: maxerr

    ! compute rhs terms 
    Call compute_rhs_u( U, V, W, rhs_uf(2:nx-1, 2:nyg-1,2:nzg-1) )
    Call compute_rhs_v( U, V, W, rhs_vf(2:nxg-1,2:ny-1, 2:nzg-1) )
    Call compute_rhs_w( U, V, W, rhs_wf(2:nxg-1,2:nyg-1,2:nz-1 ) )

    ! compute rhs at the boundaries
    Call apply_periodic_bc_x   ( rhs_uf, 1 )
    Call compute_rhs_v_boundary( rhs_vf    )
    Call apply_periodic_bc_z   ( rhs_wf, 3 )

    ! compute divergence of rhs terms
    rhs_p = 0d0
    Do k = 2, nzg-1
       Do j = 2, nyg-1
          Do i = 2, nxg-1
             rhs_p(i,j,k)  = ( rhs_uf(i,j,k) - rhs_uf(i-1,j,k) )/( x(i)-x(i-1) ) + & ! d rhs_u/dx
                             ( rhs_vf(i,j,k) - rhs_vf(i,j-1,k) )/( y(j)-y(j-1) ) + & ! d rhs_v/dy
                             ( rhs_wf(i,j,k) - rhs_wf(i,j,k-1) )/( z(k)-z(k-1) )     ! d rhs_w/dz
          End Do
       End Do
    End Do

  End Subroutine compute_pressure_rhs

  !----------------------------------------------------!
  !              Solve pressure equation               !
  !                                                    !
  !   Equation:                                        !
  !     Laplacian(p) = rhs_p = div(rhs terms)          !
  !                                                    !
  !   Boundary conditions: (valid if V=0 at the wall)  !
  !     dP/dy = d/dy(nu dV/dy) at bottom wall          !
  !     dP/dy = d/dy(nu dV/dy) at top    wall          !
  !                                                    !
  !         This is particular for channels:           !
  !               Fast Poisson solver                  !
  !                                                    !
  ! Input:  rhs_p (right-hand side for the pressure)   !
  ! Output: P                                          !
  !                                                    !
  ! Note:                  doesn't affect prms         !
  !                                |                   !
  !       p_total = p - dPdx*x - <V^2>                 !
  !                 |     |        |                   !
  !              (kx,kz) Grad    (0,0) mode            !
  !----------------------------------------------------!
  Subroutine solve_poisson_equation_for_pressure

    Integer(Int32) :: i, j, k, i_global, k_global, info
    Real   (Int64) :: maxerr

    External :: zgesv

    ! 2D Fourier transform interior points rhs_p
    Do j = 2, nyg-1
       plane = dcmplx( rhs_p(2:nxp+1,j,2:nzp+1) ) ! nxp+1 = nxg-2, nzp+1 = nzg-2
       Call fftw_mpi_execute_dft(plan_d,plane,plane_hat)
       rhs_p_hat(:,j,:) = plane_hat
    End Do

    ! 2D Fourier transform boundary conditions
    plane = dcmplx( bc_1(2:nxp+1,2:nzp+1) )
    Call fftw_mpi_execute_dft(plan_d,plane,plane_hat)
    bc_1_hat = plane_hat

    plane = dcmplx( bc_2(2:nxp+1,2:nzp+1) )
    Call fftw_mpi_execute_dft(plan_d,plane,plane_hat)
    bc_2_hat = plane_hat

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
          ! remove singularity 00 mode (set a reference pressure)
          if ( i_global==0 .And. k_global==0 ) D(2) = 3d0/2d0*D(2)
          ! rhs with boundary conditions for pressure
          rhs_aux        = rhs_p_hat(i,:,k) 
          rhs_aux(    2) = rhs_aux(    2) + coef_bc_1*(yg(  2)-yg(    1))*bc_1_hat(i,k)
          rhs_aux(nyg-1) = rhs_aux(nyg-1) - coef_bc_2*(yg(nyg)-yg(nyg-1))*bc_2_hat(i,k)
          ! solve M*u = rhs (solution stored in rhs_p_hat)
          Call Zgtsv( nr, nrhs, DL, D, DU, rhs_aux, nr, info)
          rhs_p_hat(i,:,k) = rhs_aux
       End Do
    End Do

    ! 2D inverse Fourier transform
    Do j = 2, nyg-1
       plane_hat = rhs_p_hat(:,j,:)
       Call fftw_mpi_execute_dft(plan_i,plane_hat,plane)
       P(2:nxg-2,j,2:nzg-2) = plane/Real(nxp_global*nzp_global,8)
    End Do

    ! extend values to ghost cells in y
    ! p(1)   = p(2)     - ( yg(  2)-yg(    1) )*d/dy(nu dV/dy)
    ! p(nyg) = p(nyg-1) + ( yg(nyg)-yg(nyg-1) )*d/dy(nu dV/dy) 
    ! bc_1   = d/dy(nu dV/dy) bottom wall
    ! bc_2   = d/dy(nu dV/dy) top    wall
    P(:,  1,:) = P(:,    2,:) - ( yg(  2)-yg(    1) )*bc_1
    P(:,nyg,:) = P(:,nyg-1,:) + ( yg(nyg)-yg(nyg-1) )*bc_2

    ! apply periodicity in x and z
    Call apply_periodic_xz_pressure

    ! update ghost interior planes
    Call update_ghost_interior_planes_pressure   

    ! Add mean pressure gradient (not for Prms)
    !Do i = 2, nxg-1
    !   P(i,:,:) = P(i,:,:) + dPdx*xg(i)
    !End Do

  End Subroutine solve_poisson_equation_for_pressure

  !--------------------------------------------------!
  ! Periodic boundary conditions for pseudo-pressure !
  !--------------------------------------------------!
  Subroutine apply_periodic_xz_pressure

    ! apply periodicity in x (All processors, no MPI needed)
    P ( nxg-1, :, : ) = P ( 2, :, : )

    ! apply periodicity in z (Only first and last processor, MPI needed) 
    If     ( myid==0 ) Then
       buffer_p = P ( 2:nxg-1, :, 2 ) 
       ! send data to nprocs-1
       Call Mpi_send(buffer_p, (nxg-2)*(nyg-2), MPI_real8, nprocs-1, 0, &
             MPI_COMM_WORLD,ierr)
    Elseif ( myid==nprocs-1 ) Then
       ! receive data from 0
       Call Mpi_recv(buffer_p, (nxg-2)*(nyg-2), MPI_real8, 0, 0, &
            MPI_COMM_WORLD,istat,ierr)
       P ( 2:nxg-1, :, nzp+1+1 ) = buffer_p
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
    buffer_ps = P(:,:,2)  ! send buffer
    Call Mpi_sendrecv(buffer_ps, (nxg-2)*(nyg-2), Mpi_real8, sendto, tagto,        &
         buffer_pr, (nxg-2)*(nyg-2), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
         istat, ierr)   
    If ( myid/=nprocs-1 ) P(:,:,nzg) = buffer_pr ! received buffer 

  End Subroutine update_ghost_interior_planes_pressure

  !----------------------------------------------------------!
  !          Compute rhs for V at the y-boundaries           !
  !                                                          !
  ! Equation: (valid if V=0 at the wall)                     !
  !    rhs_v(boundary) = d/dy(nu dV/dy)                      !
  !                                                          !
  ! Input:  rhs_v (interior points)                          !
  ! Output: rhs_v (interior+boundary points)                 !
  !                                                          !
  !----------------------------------------------------------!
  Subroutine compute_rhs_v_boundary(rhs_v)

    Real(Int64), Dimension(1:nxg,1:ny,1:nzg), Intent(InOut) :: rhs_v
    
    ! bottom wall
    rhs_v(2:nxg-1, 1,2:nzg-1) = bc_1

    ! top wall
    rhs_v(2:nxg-1,ny,2:nzg-1) = bc_2

  End Subroutine compute_rhs_v_boundary

  !----------------------------------------------------------!
  !     Compute boundary conditions for actual pressure      !
  !                                                          !
  ! Boundary conditions: (valid if V=0 at the wall)          !
  !     dP/dy = d/dy(nu dV/dy) at bottom wall                !
  !     dP/dy = d/dy(nu dV/dy) at top    wall                !
  !                                                          !
  ! Boundary conditions general case: <- Not used here       !     
  !     dP/dy = d/dy(nu dV/dy) - dV/dt + NL_v                !
  !                                                          !           
  ! d/dy(nu dV/dy) = dnu/dy*dV/dy  +  nu*d^2V/dy^2           !
  !                                                          !
  ! Numerical approx: 2nd order                              !
  !                                                          !
  ! Input:  V, nu, nu_t                                      !
  ! Output: bc_1 ( d/dy(nu dV/dy) at bottom wall )           !
  !         bc_2 ( d/dy(nu dV/dy) at top    wall )           !
  !                                                          !
  ! NOTE: using avg_nu_t instead of nu_t                     !
  !                                                          !
  !----------------------------------------------------------!
  Subroutine compute_pressure_BC

    ! local variables
    Integer(Int32) :: i, j, k
    Real   (Int64) :: nui, dnu, dV, ddV, dzeta, g, g2, maxerr

    dzeta  = 1d0 ! arbitrary
    maxerr = 0d0

    !-------------------------------------------------!
    ! Part 1: bottom wall

    ! metric factor for first derivative dy/dzeta 
    g  = ( -3d0*y(1) + 4d0*y(2) - y(3) )/(2d0*dzeta)

    ! metric factor for second derivative d^2y/dzeta^2
    g2 = ( 2d0*y(1) - 5d0*y(2) + 4d0*y(3) - 1d0*y(4) )/dzeta**2d0

    Do k = 2, nzg-1
       Do i = 2, nxg-1

          ! interpolated nu at faces 
          nui = nu 

          ! derivative of nu at the wall 
          dnu = 0d0 

          ! derivative of V at the wall 
          dV = ( -3d0*V(i,1,k) + 4d0*V(i,2,k) - V(i,3,k) )/(2d0*dzeta)*1d0/g 

          ! second derivative of V at the wall 
          ddV = ( 2d0*V(i,1,k) - 5d0*V(i,2,k) + 4d0*V(i,3,k) - 1d0*V(i,4,k) )/dzeta**2d0
          ddV = ( ddV - dV*g2 )/g**2d0 

          ! boundary conditions
          bc_1(i,k) = dnu*dV + nui*ddV

       End Do
    End Do

    !-------------------------------------------------!
    ! Part 2: top wall

    ! metric factor for first derivative dy/dzeta 
    g  = ( 3d0*y(ny) - 4d0*y(ny-1) + y(ny-2) )/(2d0*dzeta)

    ! metric factor for second derivative d^2y/dzeta^2 
    g2 = ( 2d0*y(ny) - 5d0*y(ny-1) + 4d0*y(ny-2) - 1d0*y(ny-3) )/dzeta**2d0

    Do k = 2, nzg-1
       Do i = 2, nxg-1

          ! interpolated nu at faces 
          nui = nu 

          ! derivative of nu at the wall 
          dnu = 0d0 

          ! derivative of V at the wall 
          dV = ( 3d0*V(i,ny,k) - 4d0*V(i,ny-1,k) + V(i,ny-2,k) )/(2d0*dzeta)*1d0/g 

          ! second derivative of V at the wall 
          ddV = ( 2d0*V(i,ny,k) - 5d0*V(i,ny-1,k) + 4d0*V(i,ny-2,k) - 1d0*V(i,ny-3,k) )/dzeta**2d0
          ddV = ( ddV - dV*g2 )/g**2d0 

          ! boundary conditions
          bc_2(i,k) = dnu*dV + nui*ddV

       End Do
    End Do

  End Subroutine compute_pressure_BC

End Module pressure

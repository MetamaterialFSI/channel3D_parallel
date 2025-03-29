!----------------------------------------------------!
!     Module for immersed boundary operators         !
!----------------------------------------------------!
Module immersed_boundary_operators

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use boundary_conditions

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------------------------------------------------------------------------!
  !      Compute the DDF weights and flow grid indices used in the IB regularization and interpolation       !
  !----------------------------------------------------------------------------------------------------------!
  Subroutine setup_IB_operators
    Integer(Int32) :: i, j, k, l, ii, jj, kk, count

    !-----------------Find the indices of the flow grid points closest to the body points---------------------!
    Do l = 1, nb
      ! index of largest x value that is smaller than xb(l)
      x_pivot_index(l) = Int(xb(l) / dx) + 1
      ! index of largest xm value that is smaller than xb(l), accounting for periodicity, but not for ghost cells
      xm_pivot_index(l) = Modulo(Int((xb(l) + 0.5 * dx) / dx) - 1, nxm_global - 1) + 1
      ! index of largest y value that is smaller than yb(l)
      y_pivot_index(l) = Int((yb(l) - (1.d0 - nd * dymin)) / dymin) + (ny_global + 1) / 2 - nd
      ! index of largest ym value that is smaller than yb(l), not accounting for ghost cells
      ym_pivot_index(l) = Int((yb(l) - (1.d0 - nd * dymin + 0.5 * dymin)) / dymin) + (ny_global + 1) / 2 - nd 
      ! index of largest z value that is smaller than zb(l)
      z_pivot_index(l) = Int(zb(l) / dz) + 1
      ! index of largest zm value that is smaller than zb(l), accounting for periodicity, but not for ghost cells
      zm_pivot_index(l) = Modulo(Int((zb(l) + 0.5 * dz) / dz) - 1, nzm_global - 1) + 1 
    End Do

    !----Find the flow grid indices and corresponding DDF weights within the support of the DDF at each body point----!
    ! u velocity locations
    Do l = 1, nb
      count = 0
      Do k = -suppz, suppz
        Do j = -suppy, suppy
          Do i = -suppx, suppx
            ii = Modulo(x_pivot_index(l) + i - 2, nx_global - 2 ) + 2
            jj = ym_pivot_index(l) + j
            kk = Modulo(zm_pivot_index(l) + k - 1, nzm_global - 1) + 1

            count = count + 1
            u_x_indices(count, l) = ii
            u_y_indices(count, l) = jj + 1 ! plus one for ghost cell
            u_z_indices(count, l) = kk + 1 ! plus one for ghost cell

            u_weights(count, l) =  dx * dymin * dz &
              * deltafnc( x_global(ii), xb(l), dx,    Lxp) &
              * deltafnc(ym_global(jj), yb(l), dymin, 0.d0) &
              * deltafnc(zm_global(kk), zb(l), dz,    Lzp)
          End Do
        End Do
      End Do
    End Do
  
    ! v velocity locations
    Do l = 1, nb
      count = 0
      Do k = -suppz, suppz
        Do j = -suppy, suppy
          Do i = -suppx, suppx
            ii = Modulo(xm_pivot_index(l) + i - 1, nxm_global - 1) + 1
            jj = y_pivot_index(l) + j
            kk = Modulo(zm_pivot_index(l) + k - 1, nzm_global - 1) + 1

            count = count+1
            v_x_indices(count, l) = ii + 1 ! plus one for ghost cell
            v_y_indices(count, l) = jj
            v_z_indices(count, l) = kk + 1 ! plus one for ghost cell

            v_weights(count, l) = dx * dymin * dz &
              * deltafnc(xm_global(ii), xb(l), dx,    Lxp) &
              * deltafnc( y_global(jj), yb(l), dymin, 0.d0) &
              * deltafnc(zm_global(kk), zb(l), dz,    Lzp)
          End Do
        End Do
      End Do
    End Do
  
    ! w velocity locations
    Do l = 1, nb
      count = 0
      Do k = -suppz, suppz
        Do j = -suppy, suppy
          Do i = -suppx, suppx
            ii = Modulo(xm_pivot_index(l) + i - 1, nxm_global - 1) + 1
            jj = ym_pivot_index(l) + j
            kk = Modulo(z_pivot_index(l) + k - 2, nz_global - 2 ) + 2

            count = count + 1
            w_x_indices(count, l) = ii + 1 ! plus one for ghost cell
            w_y_indices(count, l) = jj + 1 ! plus one for ghost cell
            w_z_indices(count, l) = kk

            w_weights(count, l) = dx * dymin * dz &
              * deltafnc(xm_global(ii), xb(l), dx,    Lxp) &
              * deltafnc(ym_global(jj), yb(l), dymin, 0.d0) &
              * deltafnc( z_global(kk), zb(l), dz,    Lzp)
          End Do
        End Do
      End Do
    End Do

  End Subroutine setup_IB_operators

  !-----------------------------------------------!
  !     Value of the Yang3 1D Discrete Delta      !
  !     function at a distance (x1-x2) from       !
  !     the center, accounting for periodicity    !
  !                                               !
  ! Input:                                        !
  !  - x1: point 1                                !
  !  - x2: point 2                                !
  !  - dr: grid spacing                           !
  !  - L: periodic length of the domain           !
  ! Output: deltafnc                              !
  !                                               !
  !-----------------------------------------------!
  Function deltafnc(x1, x2, dr, L)

    Real(Int64) :: x1, x2, r, dr, deltafnc, L, r1, r2

    r = periodic_distance(x1, x2, L)

    r1 = r / dr
    r2 = r1 * r1

    If (r1 .le. 1.d0) Then
      deltafnc = 17d0 / 48d0 + sqrt(3d0) * pi / 108d0 + &
              r1 / 4d0 - r2/4d0 + (1d0 - 2d0 * r1) / 16d0 * sqrt(-12d0 * r2 + 12d0 * &
              r1 + 1d0) - sqrt(3d0) / 12d0 * asin(sqrt(3d0) / 2d0 * (2d0 * r1 - 1d0))
    ElseIf (r1 .le. 2.d0) Then
      deltafnc = 55d0 / 48d0 - sqrt(3d0) * pi / 108 - &
              13d0 * r1 / 12d0 + r2 / 4d0 + (2d0 * r1 - 3d0) / 48d0 * &
              sqrt(-12d0 * r2 + 36d0 * r1 - 23d0) + sqrt(3d0) / 36d0 * &
              asin(sqrt(3d0) / 2d0 * (2d0 * r1 - 3d0))
    Else
      deltafnc = 0.d0
    End If

    deltafnc = deltafnc / dr

  End Function deltafnc

  ! This function is not ideal, since it might be the root of the problem why small grid sizes don't give weights that sum to one
  Function periodic_distance(x1, x2, L)
    Real   (Int64), INTENT(IN) :: x1, x2, L
    Real   (Int64) periodic_distance
    periodic_distance = min( abs( x1 - x2 ), abs( x1 - x2 - L ), abs ( x1 - x2 + L) )
  End Function periodic_distance

  Subroutine global_regT

    Call apply_boundary_conditions

    ! Gather U, V, W fields into global fields
    Call MPI_Gatherv(U(:, :, 2:nzm+1), local_size_U, MPI_REAL8, &
                 U_global(:, :, 2:nzm_global+1), send_counts_U, displs_U, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Gatherv(V(:, :, 2:nzm+1), local_size_V, MPI_REAL8, &
                 V_global(:, :, 2:nzm_global+1), send_counts_V, displs_V, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Gatherv(W(:, :, 2:nz-1), local_size_W, MPI_REAL8, &
                 W_global(:, :, 2:nz_global-1), send_counts_W, displs_W, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

    fb = 0d0
    If (myid==0) Then
      fb = regT(U_global, V_global, W_global)
    End If

  End Subroutine global_regT

  Subroutine global_reg

    fb = 1d0
    If (myid==0) Then
      U_global = regu(fb)
      V_global = regv(fb)
      W_global = regw(fb)
    End If

    ! Scatter U, V, W fields to local fields
    Call MPI_Scatterv(U_global(:, :, 2:nzm_global+1), send_counts_U, displs_U, MPI_REAL8, &
                  U_reg(:, :, 2:nzm+1), local_size_U, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Scatterv(V_global(:, :, 2:nzm_global+1), send_counts_V, displs_V, MPI_REAL8, &
                  V_reg(:, :, 2:nzm+1), local_size_V, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    Call MPI_Scatterv(W_global(:, :, 2:nz_global-1), send_counts_W, displs_W, MPI_REAL8, &
                  W_reg(:, :, 2:nz-1), local_size_W, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)


    ! This call only affects U, V, W. Should change this to apply to U_reg, etc
    Call apply_boundary_conditions

  End Subroutine global_reg

  !-----------------------------------------------!
  !          Interpolation to body points         !
  !                                               !
  ! Input: U_, V_, W_                             !
  ! Output: regT                                  !
  !                                               !
  !-----------------------------------------------!
  Function regT(U_, V_, W_)
    Implicit None
    Real(Int64), Dimension(nx_global, nyg_global, nzg_global), Intent(In) :: U_
    Real(Int64), Dimension(nxg_global, ny_global, nzg_global), Intent(In) :: V_
    Real(Int64), Dimension(nxg_global, nyg_global, nz_global), Intent(In) :: W_
    Real(Int64), Dimension(3 * nb):: regT

    Integer(Int32) :: i, j
    regT = 0.D0

    Do j = 1, nb
      Do i = 1, nweights
        regT(j         ) = regT(j)          + u_weights(i, j) * U_(u_x_indices(i, j), u_y_indices(i, j), u_z_indices(i, j))
        regT(j + nb    ) = regT(j + nb)     + v_weights(i, j) * V_(v_x_indices(i, j), v_y_indices(i, j), v_z_indices(i, j))
        regT(j + 2 * nb) = regT(j + 2 * nb) + w_weights(i, j) * W_(w_x_indices(i, j), w_y_indices(i, j), w_z_indices(i, j))
      End Do
    End Do

  End Function regT

  !-----------------------------------------------!
  !      Regularization of the x-component of     ! 
  !      f_ to U grid points                      !
  !                                               !
  ! Input: f_                                     !
  ! Output: regu                                  !
  !                                               !
  !-----------------------------------------------!
  Function regu(f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nx_global, nyg_global, nzg_global) :: regu
    Integer(Int32) :: i, j

    regu = 0.D0

    Do j = 1, nb
      Do i = 1, nweights
        regu(u_x_indices(i, j), u_y_indices(i, j), u_z_indices(i, j)) = &
          regu(u_x_indices(i, j), u_y_indices(i, j), u_z_indices(i, j)) &
           + u_weights(i, j) * 1.d0 / (dx * dymin * dz) * sb(j) * f_(j)
      End Do
    End Do

  End Function regu

  !-----------------------------------------------!
  !      Regularization of the y-component of     ! 
  !      f_ to V grid points                      !
  !                                               !
  ! Input: f_                                     !
  ! Output: regv                                  !
  !                                               !
  !-----------------------------------------------!
  Function regv(f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nxg_global, ny_global, nzg_global) :: regv
    Integer(Int32) :: i, j

    regv = 0.D0

    Do j = 1, nb
      Do i = 1, nweights
        regv(v_x_indices(i, j), v_y_indices(i, j), v_z_indices(i, j)) = &
          regv(v_x_indices(i, j), v_y_indices(i, j), v_z_indices(i, j)) &
           + v_weights(i, j) * 1.d0 / (dx * dymin * dz) * sb(j) * f_(j + nb)
      End Do
    End Do

  End Function regv
  
  !-----------------------------------------------!
  !      Regularization of the z-component of     ! 
  !      f_ to W grid points                      !
  !                                               !
  ! Input: f_                                     !
  ! Output: regw                                  !
  !                                               !
  !-----------------------------------------------!
  Function regw(f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nxg_global, nyg_global, nz_global) :: regw
    Integer(Int32) :: i, j

    regw = 0.D0

    Do j = 1, nb
      Do i = 1, nweights
        regw(w_x_indices(i, j), w_y_indices(i, j), w_z_indices(i, j)) = &
          regw(w_x_indices(i, j), w_y_indices(i, j), w_z_indices(i, j)) &
           + w_weights(i, j) * 1.d0 / (dx * dymin * dz) * sb(j) * f_(j + 2 * nb)
      End Do
    End Do
  
  End Function regw

  Subroutine divergence(diV_, U_, V_, W_)

    Integer(Int32) :: i, j, k
    Real   (Int64), DIMENSION( 2:nxg-1, 2:nyg-1, 2:nzg-1 ), INTENT(OUT) :: div_
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
    ! write(*,*) 'diV_ = ', maxval(abs(diV_))

  End Subroutine divergence

End Module immersed_boundary_operators

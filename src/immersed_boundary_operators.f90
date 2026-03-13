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
    Integer(Int32) :: i, j, k, l, ii, ii_periodic, jj, kk, kk_periodic, count, x_periodic_shifts, z_periodic_shifts
    character(len=50) :: filename
    integer :: k_local,k_supp         ! local index
    integer :: proc_id      ! MPI rank that owns this index
    Real(Int64) :: x_grid, y_grid, z_grid
    integer :: prev, next, k_global

    !-----------------Find the indices of the flow grid points closest to the body points---------------------!
    Do l = nb_start, nb_end
      ! index of largest x value that is smaller than xb(l)
      x_pivot_index(l) = Floor(xb(l) / dx) + 1
      ! index of largest xm value that is smaller than xb(l), not accounting for periodicity. Can be zero
      xm_pivot_index(l) = Floor((xb(l) + 0.5 * dx) / dx)
      ! index of largest y value that is smaller than yb(l). The use of reference points is more efficient here than searching
      ! the y array for the closest point, in case there is grid stretching.
      y_pivot_index(l) = Floor((yb(l) - y(y_ref_index(l))) / dymin) + y_ref_index(l)
      ! index of largest ym value that is smaller than yb(l), not accounting for ghost cells
      ym_pivot_index(l) = Floor((yb(l) - (y(y_ref_index(l)) + 0.5 * dymin)) / dymin) + y_ref_index(l)
      ! index of largest z value that is smaller than zb(l)
      z_pivot_index(l) = Floor(zb(l) / dz) + 1
      ! index of largest zm value that is smaller than zb(l), not accounting for periodicity. Can be zero
      zm_pivot_index(l) = Floor((zb(l) + 0.5 * dz) / dz)+1
    End Do

    !----Find the flow grid indices and corresponding DDF weights within the support of the DDF at each body point----!
    ! u velocity locations
    Do l = nb_start, nb_end
      count = 0
      Do k = -suppz, suppz
        Do j = -suppy, suppy
          Do i = -suppx, suppx
            ii = x_pivot_index(l) + i
            jj = ym_pivot_index(l) + j
            kk = zm_pivot_index(l) + k 

            x_periodic_shifts = Floor(Real(ii - 1, Int64) / (nx_global - 2))
            ii_periodic = ii - x_periodic_shifts * (nx_global - 2)
            z_periodic_shifts = Floor(Real(kk - 2, Int64) / (nzm_global))
            kk_periodic = kk - z_periodic_shifts * (nzm_global-1) !shift z index according to periodicity

            count = count + 1
            If (ii_periodic .eq. 1) Then
              u_x_indices(count, l) = nx_global - 1 ! because of how the periodic boundary conditions are set
            Else
              u_x_indices(count, l) = ii_periodic
            End If
            u_y_indices(count, l) = jj + 1 ! plus one for ghost cell
            if ( (kk_periodic .eq. 2) .and. (myid .eq. nprocs-1) ) then
              ! kk_periodic =2 at the last processor is stored in the local data, not support cell
              u_z_indices(count, l) = nzg_global - 1 ! due to periodicity
            else
              u_z_indices(count, l) = kk_periodic
            end if
            
            ! calculate the support cell index
            If ( moving_z_flag .Or. istep <= 1 ) Then
              call global_to_local_center(u_z_indices(count, l), k_supp, k_local, proc_id)
              u_z_local_indices(count, l) = k_local
              u_z_supp_idx(count, l) = k_supp
              u_proc(count, l) = proc_id 
            End If 

            x_grid = x_global(ii_periodic) + x_periodic_shifts * Lxp
            y_grid = ym_global(jj)
            z_grid = zm_global(kk_periodic - 1) + z_periodic_shifts * Lzp

            u_weights(count, l) =  dx * dymin * dz &
              * deltafnc(x_grid, xb(l), dx) &
              * deltafnc(y_grid, yb(l), dymin) &
              * deltafnc(z_grid, zb(l), dz)

            dxnu(count, l) = &
              normals(l)          * (x_grid - xb(l)) + &
              normals(nb + l)     * (y_grid - yb(l)) + &
              normals(2 * nb + l) * (z_grid - zb(l))
          End Do
        End Do
      End Do
    End Do
  
    ! v velocity locations
    Do l = nb_start, nb_end
      count = 0
      Do k = -suppz, suppz
        Do j = -suppy, suppy
          Do i = -suppx, suppx
            ii = xm_pivot_index(l) + i
            jj = y_pivot_index(l) + j
            kk = zm_pivot_index(l) + k

            x_periodic_shifts = Floor(Real(ii - 1, Int64) / (nxm_global - 1))
            ii_periodic = ii - x_periodic_shifts * (nxm_global - 1)
            z_periodic_shifts = Floor(Real(kk - 2, Int64) / (nzm_global))
            kk_periodic = kk - z_periodic_shifts * (nzm_global-1)

            count = count + 1
            v_x_indices(count, l) = ii_periodic + 1 ! plus one for ghost cell
            v_y_indices(count, l) = jj
            if ( (kk_periodic .eq. 2) .and. (myid .eq. nprocs-1) ) then
              v_z_indices(count, l) = nzg_global - 1 ! due to periodicity
            else
              v_z_indices(count, l) = kk_periodic
            end if

            ! calculate the support cell index
            If ( moving_z_flag .Or. istep <= 1 ) Then
              call global_to_local_center(v_z_indices(count, l), k_supp, k_local, proc_id)
              v_z_local_indices(count, l) = k_local
              v_z_supp_idx(count, l) = k_supp
              v_proc(count, l) = proc_id 
            End If 

            x_grid = xm_global(ii_periodic) + x_periodic_shifts * Lxp
            y_grid = y_global(jj)
            z_grid = zm_global(kk_periodic - 1) + z_periodic_shifts * Lzp

            v_weights(count, l) = dx * dymin * dz &
              * deltafnc(x_grid, xb(l),    dx) &
              * deltafnc(y_grid, yb(l), dymin) &
              * deltafnc(z_grid, zb(l),    dz)

            dxnv(count, l) = &
              normals(l)          * (x_grid - xb(l)) + &
              normals(nb + l)     * (y_grid - yb(l)) + &
              normals(2 * nb + l) * (z_grid - zb(l))
          End Do
        End Do
      End Do
    End Do

    ! w velocity locations
    Do l = nb_start, nb_end
      count = 0
      Do k = -suppz, suppz
        Do j = -suppy, suppy
          Do i = -suppx, suppx
            ii = xm_pivot_index(l) + i - 1
            jj = ym_pivot_index(l) + j
            kk = z_pivot_index(l) + k

            x_periodic_shifts = Floor(Real(ii - 1, Int64) / (nxm_global - 1))
            ii_periodic = ii - x_periodic_shifts * (nxm_global - 1)
            z_periodic_shifts = Floor(Real(kk - 2, Int64) / (nz_global - 2))
            kk_periodic = kk - z_periodic_shifts * (nz_global - 2)

            count = count + 1 
            w_x_indices(count, l) = ii_periodic + 1 ! plus one for ghost cell
            w_y_indices(count, l) = jj + 1 ! plus one for ghost cell
            w_z_indices(count, l) = kk_periodic

            ! calculate the support cell index
            If ( moving_z_flag .Or. istep <= 1 ) Then
              call global_to_local_face(w_z_indices(count, l), k_supp, k_local, proc_id)
              w_z_local_indices(count, l) = k_local
              w_z_supp_idx(count, l) = k_supp
              w_proc(count, l) = proc_id 
            End If 

            x_grid = xm_global(ii_periodic) + x_periodic_shifts * Lxp
            y_grid = ym_global(jj)
            z_grid = z_global(kk_periodic) + z_periodic_shifts * Lzp

            w_weights(count, l) = dx * dymin * dz &
              * deltafnc(x_grid, xb(l),    dx) &
              * deltafnc(y_grid, yb(l), dymin) &
              * deltafnc(z_grid, zb(l),    dz)

            dxnw(count, l) = &
              normals(l)          * (x_grid - xb(l)) + &
              normals(nb + l)     * (y_grid - yb(l)) + &
              normals(2 * nb + l) * (z_grid - zb(l))
          End Do
        End Do
      End Do
    End Do

    ! compute weightc (pressure, Gux, Gvy, Gwz locations)
    Do l = nb_start, nb_end
      count = 0
      Do k = -suppz, suppz
        Do j = -suppy, suppy
          Do i = -suppx, suppx
            ii = xm_pivot_index(l) + i
            jj = ym_pivot_index(l) + j
            kk = zm_pivot_index(l) + k

            x_periodic_shifts = Floor(Real(ii - 1, Int64) / (nxm_global - 1))
            ii_periodic = ii - x_periodic_shifts * (nxm_global - 1)
            z_periodic_shifts = Floor(Real(kk - 2, Int64) / (nzm_global))
            kk_periodic = kk - z_periodic_shifts * (nzm_global-1)

            count = count + 1
            c_x_indices(count, l) = ii_periodic + 1 ! plus one for ghost cell
            c_y_indices(count, l) = jj + 1 ! plus one for ghost cell
            if ( (kk_periodic .eq. 2) .and. (myid .eq. nprocs-1) ) then
              c_z_indices(count, l) = nzg_global - 1 ! due to periodicity
            else
              c_z_indices(count, l) = kk_periodic ! plus one for ghost cell
            end if

            ! calculate the support cell index
            If ( moving_z_flag .Or. istep <= 1 ) Then
              call global_to_local_center(c_z_indices(count, l), k_supp, k_local, proc_id)
              c_z_local_indices(count, l) = k_local
              c_z_supp_idx(count, l) = k_supp
              c_proc(count, l) = proc_id 
            End If 

            x_grid = xm_global(ii_periodic) + x_periodic_shifts * Lxp
            y_grid = ym_global(jj)
            z_grid = zm_global(kk_periodic-1) + z_periodic_shifts * Lzp

            c_weights(count, l) = dx * dymin * dz &
              * deltafnc(x_grid, xb(l),    dx) &
              * deltafnc(y_grid, yb(l), dymin) &
              * deltafnc(z_grid, zb(l),    dz)

            dxnc(count, l) = &
              normals(l)          * (x_grid - xb(l)) + &
              normals(nb + l)     * (y_grid - yb(l)) + &
              normals(2 * nb + l) * (z_grid - zb(l))
          End Do
        End Do
      End Do
    End Do

    local_size_nb = nb_end - nb_start + 1

    ! Gather send_counts
    Call MPI_Allgather(local_size_nb, 1, MPI_INT, send_counts_nb, 1, MPI_INT, MPI_COMM_WORLD, ierr)

    displs_nb(1) = 0
    Do i = 2, nprocs
      displs_nb(i) = displs_nb(i-1) + send_counts_nb(i-1)
    End Do

  End Subroutine setup_IB_operators

  !-----------------------------------------------!
  !     Value of the Yang3 1D Discrete Delta      !
  !     function at a distance (x1-x2) from       !
  !     the center                                !
  !                                               !
  ! Input:                                        !
  !  - x1: point 1                                !
  !  - x2: point 2                                !
  !  - dr: grid spacing                           !
  ! Output: deltafnc                              !
  !                                               !
  !-----------------------------------------------!
  Function deltafnc(x1, x2, dr)

    Real(Int64) :: x1, x2, r, dr, deltafnc, r1, r2

    r = abs(x1 - x2)

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

  Function regT(U_, V_, W_)
    Implicit None
    Real(Int64), CONTIGUOUS, INTENT(INOUT)  :: U_(:, :, :)
    Real(Int64), CONTIGUOUS, INTENT(INOUT)  :: V_(:, :, :)
    Real(Int64), CONTIGUOUS, INTENT(INOUT)  :: W_(:, :, :)
    Real(Int64), Dimension(3 * nb):: regT

    ! update support cells from interior points
    Call interior_planes_update_support(U_, U_supp, 1)
    Call interior_planes_update_support(V_, V_supp, 2)
    Call interior_planes_update_support(W_, W_supp, 3)

    ! compute local_regT
    regT_buffer_vector = local_regT(U_, V_, W_, U_supp, V_supp, W_supp)
    ! Gather regT values from all partitions
    Call MPI_Allgatherv(regT_buffer_vector(nb_start : nb_end), local_size_nb, MPI_REAL8, &
                 regT(1 : nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    Call MPI_Allgatherv(regT_buffer_vector(nb + nb_start : nb + nb_end), local_size_nb, MPI_REAL8, &
                 regT(nb + 1 : 2 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    Call MPI_Allgatherv(regT_buffer_vector(2 * nb + nb_start : 2 * nb + nb_end), local_size_nb, MPI_REAL8, &
                 regT(2 * nb + 1 : 3 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)

  End Function regT

  Function regT_1n(U_, V_, W_)
    Implicit None
    Real(Int64), CONTIGUOUS, INTENT(INOUT)  :: U_(:, :, :)
    Real(Int64), CONTIGUOUS, INTENT(INOUT)  :: V_(:, :, :)
    Real(Int64), CONTIGUOUS, INTENT(INOUT)  :: W_(:, :, :)
    Real(Int64), Dimension(3 * nb):: regT_1n

    ! update support cells from interior points
    Call interior_planes_update_support(U_, U_supp, 1)
    Call interior_planes_update_support(V_, V_supp, 2)
    Call interior_planes_update_support(W_, W_supp, 3)

    ! compute local_regT
    regT_buffer_vector = local_regT_1n(U_, V_, W_, U_supp, V_supp, W_supp)
    ! Gather regT values from all partitions
    Call MPI_Allgatherv(regT_buffer_vector(nb_start : nb_end), local_size_nb, MPI_REAL8, &
                 regT_1n(1 : nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    Call MPI_Allgatherv(regT_buffer_vector(nb + nb_start : nb + nb_end), local_size_nb, MPI_REAL8, &
                 regT_1n(nb + 1 : 2 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    Call MPI_Allgatherv(regT_buffer_vector(2 * nb + nb_start : 2 * nb + nb_end), local_size_nb, MPI_REAL8, &
                 regT_1n(2 * nb + 1 : 3 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)

  End Function regT_1n

  Function regTc_1n(P_)
    Implicit None
    !Real(Int64), CONTIGUOUS, INTENT(INOUT)  :: P_(:, :, :)
    Real(Int64), Dimension(2:nxg-1, 2:nyg-1, 2:nzg ), Intent(INOUT) :: P_
    Real(Int64), Dimension(nb):: regTc_1n

    ! update support cells from interior points
    !Call interior_planes_update_support(P_, P_supp, 4)
    Call interior_planes_update_support_pressure(P_, P_supp)

    ! compute local_regT
    regT_buffer_scalar = local_regTc_1n(P_, P_supp)
    ! Gather regT values from all partitions
    Call MPI_Allgatherv(regT_buffer_scalar(nb_start : nb_end), local_size_nb, MPI_REAL8, &
                 regTc_1n(1 : nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)

  End Function regTc_1n

  Subroutine regu(U_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nx, nyg, nzg), Intent(Out) :: U_

    U_ = 0d0

    ! compute local_reg for each partition
    Call local_reg(U_, U_supp,f_, 1)
    ! update the interior points from other partitions
    Call support_update_interior_planes(U_, U_supp, 1)

  End Subroutine regu

  Subroutine regv(V_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nxg, ny, nzg), Intent(Out) :: V_

    V_ = 0d0
    
    ! compute local_reg for each partition
    Call local_reg(V_, V_supp, f_, 2)
    ! update the interior points from other partitions
    Call support_update_interior_planes(V_, V_supp, 2)

  End Subroutine regv

  Subroutine regw(W_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nxg, nyg, nz), Intent(Out) :: W_
    Integer(Int32) :: i, j

    W_ = 0d0

    ! compute local_reg for each partition
    Call local_reg(W_,W_supp,f_,3)
    ! update the interior points from other partitions
    Call support_update_interior_planes(W_,W_supp,3)

  End Subroutine regw

  !-----------------------------------------------!
  !          Interpolation to body points         !
  !                                               !
  ! Input: U_, V_, W_                             !
  ! Output: local_regT                                  !
  !                                               !
  !-----------------------------------------------!
  Function local_regT(U_, V_, W_, U_supp_, V_supp_, W_supp_)
    Implicit None
    Real(Int64), Dimension(nx, nyg, nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny, nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg, nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nx, nyg, suppz*2+1), Intent(In) :: U_supp_
    Real(Int64), Dimension(nxg, ny, suppz*2+1), Intent(In) :: V_supp_
    Real(Int64), Dimension(nxg, nyg, suppz*2+1), Intent(In) :: W_supp_
    Real(Int64), Dimension(3 * nb):: local_regT
    integer :: proc_idx

    Integer(Int32) :: i, j
    local_regT = 0.D0
     
    Do j = nb_start, nb_end
      Do i = 1, nweights
        ! get data of U
        proc_idx = u_proc(i, j)
        if ( proc_idx .eq. myid ) then
          local_regT(j) = local_regT(j)+ &
          u_weights(i, j) * U_(u_x_indices(i, j),u_y_indices(i, j),u_z_local_indices(i, j))
        else
          local_regT(j) = local_regT(j)+ &
          u_weights(i, j) * U_supp_(u_x_indices(i, j),u_y_indices(i, j),u_z_supp_idx(i, j))
        end if
        ! get data of V
        proc_idx = v_proc(i, j)
        if ( proc_idx .eq. myid ) then
          local_regT(j + nb) = local_regT(j + nb)+ &
          v_weights(i, j) * V_(v_x_indices(i, j),v_y_indices(i, j),v_z_local_indices(i, j))
        else
          local_regT(j + nb) = local_regT(j + nb)+ &
          v_weights(i, j) * V_supp_(v_x_indices(i, j),v_y_indices(i, j),v_z_supp_idx(i, j))
        end if
        ! get data of W
        proc_idx = w_proc(i, j)
        if ( proc_idx .eq. myid ) then
          local_regT(j + 2 * nb) = local_regT(j + 2 * nb)+ &
          w_weights(i, j) * W_(w_x_indices(i, j),w_y_indices(i, j),w_z_local_indices(i, j))
        else
          local_regT(j + 2 * nb) = local_regT(j + 2 * nb)+ &
          w_weights(i, j) * W_supp_(w_x_indices(i, j),w_y_indices(i, j),w_z_supp_idx(i, j))
        end if
      End Do
    End Do
  End Function local_regT

  Function local_regT_1n(U_, V_, W_, U_supp_, V_supp_, W_supp_)
    Implicit None
    Real(Int64), Dimension(nx, nyg, nzg), Intent(In) :: U_
    Real(Int64), Dimension(nxg, ny, nzg), Intent(In) :: V_
    Real(Int64), Dimension(nxg, nyg, nz), Intent(In) :: W_
    Real(Int64), Dimension(nx, nyg, suppz*2+1), Intent(In) :: U_supp_
    Real(Int64), Dimension(nxg, ny, suppz*2+1), Intent(In) :: V_supp_
    Real(Int64), Dimension(nxg, nyg, suppz*2+1), Intent(In) :: W_supp_
    Real(Int64), Dimension(3 * nb):: local_regT_1n
    integer :: proc_idx

    Integer(Int32) :: i, j
    local_regT_1n = 0.D0
     
    Do j = nb_start, nb_end
      Do i = 1, nweights
        ! get data of U
        proc_idx = u_proc(i, j)
        if ( proc_idx .eq. myid ) then
          local_regT_1n(j) = local_regT_1n(j)+ &
          u_weights(i, j) * U_(u_x_indices(i, j),u_y_indices(i, j),u_z_local_indices(i, j)) * dxnu(i, j)
        else
          local_regT_1n(j) = local_regT_1n(j)+ &
          u_weights(i, j) * U_supp_(u_x_indices(i, j),u_y_indices(i, j),u_z_supp_idx(i, j)) * dxnu(i, j)
        end if
        ! get data of V
        proc_idx = v_proc(i, j)
        if ( proc_idx .eq. myid ) then
          local_regT_1n(j + nb) = local_regT_1n(j + nb)+ &
          v_weights(i, j) * V_(v_x_indices(i, j),v_y_indices(i, j),v_z_local_indices(i, j)) * dxnv(i, j)
        else
          local_regT_1n(j + nb) = local_regT_1n(j + nb)+ &
          v_weights(i, j) * V_supp_(v_x_indices(i, j),v_y_indices(i, j),v_z_supp_idx(i, j)) * dxnv(i, j)
        end if
        ! get data of W
        proc_idx = w_proc(i, j)
        if ( proc_idx .eq. myid ) then
          local_regT_1n(j + 2 * nb) = local_regT_1n(j + 2 * nb)+ &
          w_weights(i, j) * W_(w_x_indices(i, j),w_y_indices(i, j),w_z_local_indices(i, j)) * dxnw(i, j)
        else
          local_regT_1n(j + 2 * nb) = local_regT_1n(j + 2 * nb)+ &
          w_weights(i, j) * W_supp_(w_x_indices(i, j),w_y_indices(i, j),w_z_supp_idx(i, j)) * dxnw(i, j)
        end if
      End Do
    End Do
  End Function local_regT_1n

  Function local_regTc_1n(P_, P_supp_)
    Implicit None
    Real(Int64), Dimension(2:nxg-1, 2:nyg-1, 2:nzg ), Intent(In) :: P_
    Real(Int64), Dimension(2:nxg-1, 2:nyg-1, suppz * 2 + 1), Intent(In) :: P_supp_
    Real(Int64), Dimension(nb):: local_regTc_1n
    integer :: proc_idx

    Integer(Int32) :: i, j
    local_regTc_1n = 0.D0

    Do j = nb_start, nb_end
      Do i = 1, nweights
        ! get data of P
        proc_idx = c_proc(i, j)
        if ( proc_idx .eq. myid ) then
          local_regTc_1n(j) = local_regTc_1n(j)+ &
          c_weights(i, j) * P_(c_x_indices(i, j),c_y_indices(i, j),c_z_local_indices(i, j)) * dxnc(i, j)
        else
          local_regTc_1n(j) = local_regTc_1n(j)+ &
          c_weights(i, j) * P_supp_(c_x_indices(i, j),c_y_indices(i, j),c_z_supp_idx(i, j)) * dxnc(i, j)
        end if
      End Do
    End Do
  End Function local_regTc_1n

  !-----------------------------------------------------------------------
  !  Generalized “subset_regu” for U, V, or W:
  !    - id = 1 -> U
  !    - id = 2 -> V
  !    - id = 3 -> W
  !
  !  Inputs:
  !    id        : 1->U, 2->V, 3->W
  !    f_        : real(8) array of length nb (f_ for each body point)
  !
  !  Output:
  !    subset    : real(8)(nweights, nb)  (subset_regu for chosen field)
  !
  !  Uses global weights u_weights, v_weights, w_weights; and
  !  dx, dy, dymin, dz, sb(:), nb_start, nb_end, nweights, nb.
  !-----------------------------------------------------------------------
  subroutine local_reg(F, F_supp, f_, id)
    implicit none
    integer, intent(in)  :: id
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Intent(InOut) :: F(:,:,:)    
    Real(Int64), Intent(InOut) :: F_supp(:,:,:)   
    integer :: i, j
    integer :: proc_idx
    integer :: xi, yi, zi,zi_supp,zi_loc
    Real(Int64) :: factor, weight

    ! Initialize output to zero
    F = 0.0
    F_supp =0.0

    ! Determine scaling factor based on id
    factor = 1.0 / (dx * dymin * dz)

    ! Loop over body points and accumulate using the appropriate weight array
    do j = nb_start, nb_end
      do i = 1, nweights
        select case (id)
        case (1)
          proc_idx = u_proc(i, j)
          xi       = u_x_indices(i, j)
          yi       = u_y_indices(i, j)
          zi_loc   = u_z_local_indices(i, j)
          zi_supp  = u_z_supp_idx(i, j)
          weight = u_weights(i,j)
        case (2)
          proc_idx = v_proc(i, j)
          xi       = v_x_indices(i, j)
          yi       = v_y_indices(i, j)
          zi_loc   = v_z_local_indices(i, j)
          zi_supp  = v_z_supp_idx(i, j)
          weight = v_weights(i,j)
        case (3)
          proc_idx = w_proc(i, j)
          xi       = w_x_indices(i, j)
          yi       = w_y_indices(i, j)
          zi_loc   = w_z_local_indices(i, j)
          zi_supp  = w_z_supp_idx(i, j)
          weight = w_weights(i,j)
        case default
          print *, 'Error in scatter_IB_subset: invalid id =', id
          stop
        end select
        if (proc_idx == myid) then
          ! Write into local F
          F(xi, yi, zi_loc) = F(xi, yi, zi_loc)+weight * factor* sb(j) *f_(j+(id-1)*nb)
        else
          if (zi_supp < 1 .or. zi_supp > 2*suppz+1) Then
            WRITE(*,*) 'myid',myid,'proc_idx',proc_idx,'zi_supp',zi_supp
            stop 'Error: zi_supp out of [1..suppz]'
          END IF
          F_supp(xi, yi, zi_supp) = F_supp(xi, yi, zi_supp)+weight * factor* sb(j) *f_(j+(id-1)*nb)
          !print *, 'Error: unexpected proc_idx =', proc_idx, ' for myid =', myid
          !stop
        end if
      end do
    end do
  end subroutine local_reg

  !-------------------------------------------------------------------------------
  !  global_to_local_face
  !
  !  Purpose:
  !    Determines the owner rank, local index, and support index (if needed)
  !    of a given global face index in the z-direction.
  !
  !  Inputs:
  !    k_global : Global face index in z-direction
  !
  !  Outputs:
  !    k_sup    : Index in F_supp if the face is in a neighboring rank
  !               (1..suppz+1 from left, suppz+2..2*suppz+1 from right), -1 if local
  !    k_loc    : Local index in F if the face is local
  !    rank     : Rank (myid, prev, or next) that owns k_global
  !
  !  Notes:
  !    - Special handling when nprocs = 2 (prev == next)
  !-------------------------------------------------------------------------------
  subroutine global_to_local_face(k_global, k_sup, k_loc, rank)
    implicit none
    integer, intent(in)  :: k_global
    integer, intent(out) :: k_sup, k_loc, rank
    integer :: prev, next

    ! Initialize outputs
    rank  = -1
    k_sup = -1
    k_loc = -1

    ! Determine periodic neighbors
    prev = myid - 1
    if (prev < 0) prev = nprocs - 1
    next = myid + 1
    if (next == nprocs) next = 0

    ! Check ownership among {prev, myid, next}
    if (k_global >= k1_global(myid)+1 .and. k_global <= k2_global(myid)-1) then
      rank = myid
    elseif (k_global >= k1_global(prev)+1 .and. k_global <= k2_global(prev)-1) then
      rank = prev
    elseif (k_global >= k1_global(next)+1 .and. k_global <= k2_global(next)-1) then
      rank = next
    elseif (k_global .eq. 1) Then
        rank=0
    else
      print *, 'Error: Face index ', k_global, ' not in {', prev, ',', myid, ',', next, '}.'
      stop
    end if
    k_loc = k_global - k1_global(rank) + 1
    if (rank == myid) then
      ! Locally owned -> compute k_loc
      k_sup = -1
    else  ! rank == next or rank == prev
      if ( nprocs .eq. 2 ) then
        ! special case for only 2 processors: prev == next
        if ( k_loc>nz/2 ) then
          k_sup = k_global - (k2_global(prev) - suppz)+2
        else
          k_sup = (k_global - k1_global(next)) + (suppz+1)
        end if
      else 
        if ( rank .eq. prev ) then
          k_sup = k_global - (k2_global(prev) - suppz)+2
        elseif (rank .eq. next) then
          k_sup = (k_global - k1_global(next)) + (suppz+1)
        end if
      end if
      if (k_sup < 1 .or. k_sup > 2*suppz+1) Then
        stop 'Error: zi_supp out of [1..suppz*2+1] for face'
      END IF
    end if
  end subroutine global_to_local_face

  !-------------------------------------------------------------------------------
  !  global_to_local_center
  !
  !  Purpose:
  !    Determines the owner rank, local index, and support index (if needed)
  !    of a given global center index in the z-direction.
  !
  !  Inputs:
  !    k_global : Global center index in z-direction
  !
  !  Outputs:
  !    k_sup    : Index in F_supp if the center is in a neighboring rank
  !               (1..suppz from left, suppz+1..2*suppz from right), -1 if local
  !    k_loc    : Local index in F if the center is local
  !    rank     : Rank (myid, prev, or next) that owns k_global
  !
  !  Notes:
  !    - Special handling when nprocs = 2 (prev == next)
  !-------------------------------------------------------------------------------
  subroutine global_to_local_center(k_global, k_sup, k_loc, rank)
    implicit none
    integer, intent(in)  :: k_global
    integer, intent(out) :: k_sup, k_loc, rank
    integer :: prev, next, neighbor_top_start

    ! Initialize outputs
    rank  = -1
    k_sup = -1
    k_loc = -1

    ! Determine periodic neighbors
    prev = myid - 1
    if (prev < 0) prev = nprocs - 1
    next = myid + 1
    if (next == nprocs) next = 0

    ! Check ownership among {prev, myid, next}
    if (k_global >= (kg1_global(prev)+1) .and. &
        k_global <= (kg2_global(prev)-1)) then
      rank = prev
    elseif (k_global >= (kg1_global(myid)+1) .and. &
            k_global <= (kg2_global(myid)-1)) then
      rank = myid
    elseif (k_global >= (kg1_global(next)+1) .and. &
            k_global <= (kg2_global(next)-1)) then
      rank = next
    else
      write(*,*) 'Error: Center index ', k_global, ' not in {', prev, ',', myid, ',', next, '}.'
      stop
    end if
    k_loc = k_global - kg1_global(rank) + 1
    if (rank == myid) then
      ! Locally owned -> compute k_loc
      k_sup = -1
    else  ! rank == next
      if ( nprocs .eq. 2 ) then
        ! special case for only 2 processors: prev=next
        if ( k_loc>nzg/2 ) then
          neighbor_top_start = (kg2_global(prev) - suppz - 2)
          k_sup = k_global - neighbor_top_start
          if ( rank .eq. nprocs-1 ) then
            k_sup=k_sup+1
          end if
        else
          k_sup = (k_global - (kg1_global(next) + 1)) + (suppz + 2)
          if (rank .eq. 0) then
            k_sup=k_sup-1
          end if
        end if
      else 
        if ( rank .eq. prev ) then
          k_sup = k_global - kg2_global(prev) + suppz + 2
          if ( rank .eq. nprocs-1 ) then
            k_sup=k_sup+1
          end if
        elseif (rank .eq. next) then
          k_sup = k_global - kg1_global(next) + suppz + 1
          if (rank .eq. 0) then
            k_sup=k_sup-1
          end if
        end if
      end if
      if (k_sup < 1 .or. k_sup > 2*suppz+1) Then
        stop 'Error: zi_supp out of [1..suppz*2+1] for center'
      END IF
    end if
  end subroutine global_to_local_center

End Module immersed_boundary_operators

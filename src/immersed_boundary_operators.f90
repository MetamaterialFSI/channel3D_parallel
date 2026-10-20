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
    Real(Int64) :: x_grid, y_grid, z_grid, x_offset, y_offset, z_offset
    integer :: prev, next, k_global

    !-----------------Find the indices of the flow grid points closest to the body points---------------------!
    E2n1c = 0d0
    M2dir = 0d0

    Do l = 1, nb
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
    Do l = 1, nb
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
            If ( l >= nb_start .And. l <= nb_end .And. (moving_z_flag .Or. istep <= 1) ) Then
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
    Do l = 1, nb
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
            If ( l >= nb_start .And. l <= nb_end .And. (moving_z_flag .Or. istep <= 1) ) Then
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
    Do l = 1, nb
      count = 0
      Do k = -suppz, suppz
        Do j = -suppy, suppy
          Do i = -suppx, suppx
            ii = xm_pivot_index(l) + i
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
            If ( l >= nb_start .And. l <= nb_end .And. (moving_z_flag .Or. istep <= 1) ) Then
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
    Do l = 1, nb
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
            If ( l >= nb_start .And. l <= nb_end .And. (moving_z_flag .Or. istep <= 1) ) Then
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

            x_offset = x_grid - xb(l)
            y_offset = y_grid - yb(l)
            z_offset = z_grid - zb(l)
            dxnc(count, l) = &
              normals(l)          * x_offset + &
              normals(nb + l)     * y_offset + &
              normals(2 * nb + l) * z_offset
            E2n1c(l) = E2n1c(l) + c_weights(count, l) * dxnc(count, l) ** 2
            M2dir(l) = M2dir(l) + c_weights(count, l) * x_offset ** 2
            M2dir(nb + l) = M2dir(nb + l) + c_weights(count, l) * y_offset ** 2
            M2dir(2 * nb + l) = M2dir(2 * nb + l) + c_weights(count, l) * z_offset ** 2
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

    Call compute_ecross

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

  Real(Int64) Function dxw(a, b)
    Integer(Int32), Intent(In) :: a, b

    dxw = xb(a) - xb(b)
    If (dxw > 0.5d0 * Lxp) dxw = dxw - Lxp
    If (dxw < -0.5d0 * Lxp) dxw = dxw + Lxp
  End Function dxw

  Real(Int64) Function dzw(a, b)
    Integer(Int32), Intent(In) :: a, b

    dzw = zb(a) - zb(b)
    If (dzw > 0.5d0 * Lzp) dzw = dzw - Lzp
    If (dzw < -0.5d0 * Lzp) dzw = dzw + Lzp
  End Function dzw

  Subroutine compute_ecross
    Integer(Int32) :: k, kp, km, m, f
    Real(Int64) :: ds, tx, ty, tz

    Ecross = 0d0
    Do k = 1, nb
      f = ring_start(k)
      m = ring_len(k)
      If (m < 3) Cycle
      kp = f + Mod(k - f + 1, m)
      km = f + Mod(k - f - 1 + m, m)
      ds = Sqrt(dxw(kp, km) ** 2 + (yb(kp) - yb(km)) ** 2 + dzw(kp, km) ** 2)
      If (ds <= 0d0) Cycle
      tx = dxw(kp, km) / ds
      ty = (yb(kp) - yb(km)) / ds
      tz = dzw(kp, km) / ds
      Ecross(k) = normals(k) * tx * M2dir(k) &
                + normals(nb + k) * ty * M2dir(nb + k) &
                + normals(2 * nb + k) * tz * M2dir(2 * nb + k)
    End Do
  End Subroutine compute_ecross

  Subroutine dpmean_ds(p_, cfield, pj_, dpds, tang)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: p_(2:, 2:, 2:)
    Real(Int64), Intent(In) :: cfield
    Real(Int64), Intent(In) :: pj_(:)
    Real(Int64), Intent(Out) :: dpds(:)
    Real(Int64), Intent(Out) :: tang(:)
    Real(Int64) :: pmean(nb), ds
    Integer(Int32) :: k, kp, km, m, f

    dpds = 0d0
    tang = 0d0
    Call regTc(pmean, p_)
    pmean = cfield * pmean - (EHc_exterior - 0.5d0) * pj_

    Do k = 1, nb
      f = ring_start(k)
      m = ring_len(k)
      If (m < 3) Cycle
      kp = f + Mod(k - f + 1, m)
      km = f + Mod(k - f - 1 + m, m)
      ds = Sqrt(dxw(kp, km) ** 2 + (yb(kp) - yb(km)) ** 2 + dzw(kp, km) ** 2)
      If (ds <= 0d0) Cycle
      tang(k) = dxw(kp, km) / ds
      tang(nb + k) = (yb(kp) - yb(km)) / ds
      tang(2 * nb + k) = dzw(kp, km) / ds
      dpds(k) = (pmean(kp) - pmean(km)) / ds
    End Do

  End Subroutine dpmean_ds

  Subroutine frame_to_cart_t(f_, dudn)
    Implicit None
    Real(Int64), Intent(In) :: f_(:)
    Real(Int64), Intent(Out) :: dudn(:)
    Integer(Int32) :: ic, j0

    Do ic = 1, 3
      j0 = (ic - 1) * nb
      dudn(j0 + 1 : j0 + nb) = f_(1 : nb) * tangents_1(j0 + 1 : j0 + nb) &
                              + f_(nb + 1 : 2 * nb) * tangents_2(j0 + 1 : j0 + nb)
    End Do

  End Subroutine frame_to_cart_t

  Subroutine cart_to_frame_t(v3, out)
    Implicit None
    Real(Int64), Intent(In) :: v3(:)
    Real(Int64), Intent(Out) :: out(:)

    out(1 : nb) = v3(1 : nb) * tangents_1(1 : nb) &
                + v3(nb + 1 : 2 * nb) * tangents_1(nb + 1 : 2 * nb) &
                + v3(2 * nb + 1 : 3 * nb) * tangents_1(2 * nb + 1 : 3 * nb)
    out(nb + 1 : 2 * nb) = v3(1 : nb) * tangents_2(1 : nb) &
                         + v3(nb + 1 : 2 * nb) * tangents_2(nb + 1 : 2 * nb) &
                         + v3(2 * nb + 1 : 3 * nb) * tangents_2(2 * nb + 1 : 3 * nb)

  End Subroutine cart_to_frame_t

  Subroutine regT(Eu_, U_, V_, W_)
    Implicit None
    Real(Int64), Contiguous, Intent(In)  :: U_(:, :, :)
    Real(Int64), Contiguous, Intent(In)  :: V_(:, :, :)
    Real(Int64), Contiguous, Intent(In)  :: W_(:, :, :)
    Real(Int64), Contiguous, Intent(Out) :: Eu_(:)

    ! update support cells from interior points
    Call interior_planes_update_support(U_, U_supp, 1)
    Call interior_planes_update_support(V_, V_supp, 2)
    Call interior_planes_update_support(W_, W_supp, 3)

    ! compute local_regT
    Call local_regT(regT_buffer_vector, U_, V_, W_, U_supp, V_supp, W_supp)
    ! Gather regT values from all partitions
    Call MPI_Allgatherv(regT_buffer_vector(nb_start : nb_end), local_size_nb, MPI_REAL8, &
                 Eu_(1 : nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    Call MPI_Allgatherv(regT_buffer_vector(nb + nb_start : nb + nb_end), local_size_nb, MPI_REAL8, &
                 Eu_(nb + 1 : 2 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    Call MPI_Allgatherv(regT_buffer_vector(2 * nb + nb_start : 2 * nb + nb_end), local_size_nb, MPI_REAL8, &
                 Eu_(2 * nb + 1 : 3 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)

  End Subroutine regT

  Subroutine regT_1n(E1nu_, U_, V_, W_)
    Implicit None
    Real(Int64), Contiguous, Intent(In)  :: U_(:, :, :)
    Real(Int64), Contiguous, Intent(In)  :: V_(:, :, :)
    Real(Int64), Contiguous, Intent(In)  :: W_(:, :, :)
    Real(Int64), Contiguous, Intent(Out) :: E1nu_(:)

    ! update support cells from interior points
    Call interior_planes_update_support(U_, U_supp, 1)
    Call interior_planes_update_support(V_, V_supp, 2)
    Call interior_planes_update_support(W_, W_supp, 3)

    ! compute local_regT
    Call local_regT_1n(regT_buffer_vector, U_, V_, W_, U_supp, V_supp, W_supp)
    ! Gather regT values from all partitions
    Call MPI_Allgatherv(regT_buffer_vector(nb_start : nb_end), local_size_nb, MPI_REAL8, &
                 E1nu_(1 : nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    Call MPI_Allgatherv(regT_buffer_vector(nb + nb_start : nb + nb_end), local_size_nb, MPI_REAL8, &
                 E1nu_(nb + 1 : 2 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    Call MPI_Allgatherv(regT_buffer_vector(2 * nb + nb_start : 2 * nb + nb_end), local_size_nb, MPI_REAL8, &
                 E1nu_(2 * nb + 1 : 3 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)

  End Subroutine regT_1n

  Subroutine regTc_1n(E1np_, P_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: P_(2:, 2:, 2:)
    Real(Int64), Contiguous, Intent(Out) :: E1np_(:)

    ! update support cells from interior points
    Call interior_planes_update_support_pressure(P_, P_supp)
    ! compute local_regT
    Call local_regTc_1n(regT_buffer_scalar, P_, P_supp)
    ! Gather regT values from all partitions
    Call MPI_Allgatherv(regT_buffer_scalar(nb_start : nb_end), local_size_nb, MPI_REAL8, &
                 E1np_, send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)

  End Subroutine regTc_1n

  Subroutine regTc(Ec_, P_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: P_(2:, 2:, 2:)
    Real(Int64), Contiguous, Intent(Out) :: Ec_(:)

    Call interior_planes_update_support_pressure(P_, P_supp)
    Call local_regTc(regT_buffer_scalar, P_, P_supp)
    Call MPI_Allgatherv(regT_buffer_scalar(nb_start : nb_end), local_size_nb, MPI_REAL8, &
                 Ec_, send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)

  End Subroutine regTc

  Subroutine regu(U_, f_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: U_(:, :, :)

    ! compute local_reg for each partition
    Call local_reg(U_, U_supp, f_, 1)
    ! update the interior points from other partitions
    Call support_update_interior_planes(U_, U_supp, 1)

  End Subroutine regu

  Subroutine regv(V_, f_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: V_(:, :, :)

    ! compute local_reg for each partition
    Call local_reg(V_, V_supp, f_, 2)
    ! update the interior points from other partitions
    Call support_update_interior_planes(V_, V_supp, 2)

  End Subroutine regv

  Subroutine regw(W_, f_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: W_(:, :, :)
    ! compute local_reg for each partition
    Call local_reg(W_, W_supp, f_, 3)
    ! update the interior points from other partitions
    Call support_update_interior_planes(W_, W_supp, 3)

  End Subroutine regw

  Subroutine regu_1n(U_, f_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: U_(:, :, :)

    ! compute local_reg for each partition
    Call local_reg_1n(U_, U_supp, f_, 1)
    ! update the interior points from other partitions
    Call support_update_interior_planes(U_, U_supp, 1)

  End Subroutine regu_1n

  Subroutine regv_1n(V_, f_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: V_(:, :, :)

    ! compute local_reg for each partition
    Call local_reg_1n(V_, V_supp, f_, 2)
    ! update the interior points from other partitions
    Call support_update_interior_planes(V_, V_supp, 2)

  End Subroutine regv_1n

  Subroutine regw_1n(W_, f_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: W_(:, :, :)
    ! compute local_reg for each partition
    Call local_reg_1n(W_, W_supp, f_, 3)
    ! update the interior points from other partitions
    Call support_update_interior_planes(W_, W_supp, 3)

  End Subroutine regw_1n

  !-----------------------------------------------!
  !          Interpolation to body points         !
  !                                               !
  ! Input: U_, V_, W_                             !
  ! Output: local_regT                                  !
  !                                               !
  !-----------------------------------------------!
  Subroutine local_regT(Eu_, U_, V_, W_, U_supp_, V_supp_, W_supp_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: U_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: V_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: W_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: U_supp_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: V_supp_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: W_supp_(:, :, :)
    Real(Int64), Contiguous, Intent(Out) :: Eu_(:)
    integer :: proc_idx

    Integer(Int32) :: i, j
    Eu_(nb_start : nb_end) = 0d0
    Eu_(nb + nb_start : nb + nb_end) = 0d0
    Eu_(2 * nb + nb_start : 2 * nb + nb_end) = 0d0
     
    Do j = nb_start, nb_end
      Do i = 1, nweights
        ! get data of U
        proc_idx = u_proc(i, j)
        if ( proc_idx .eq. myid ) then
          Eu_(j) = Eu_(j)+ &
          u_weights(i, j) * U_(u_x_indices(i, j),u_y_indices(i, j),u_z_local_indices(i, j))
        else
          Eu_(j) = Eu_(j)+ &
          u_weights(i, j) * U_supp_(u_x_indices(i, j),u_y_indices(i, j),u_z_supp_idx(i, j))
        end if
        ! get data of V
        proc_idx = v_proc(i, j)
        if ( proc_idx .eq. myid ) then
          Eu_(j + nb) = Eu_(j + nb)+ &
          v_weights(i, j) * V_(v_x_indices(i, j),v_y_indices(i, j),v_z_local_indices(i, j))
        else
          Eu_(j + nb) = Eu_(j + nb)+ &
          v_weights(i, j) * V_supp_(v_x_indices(i, j),v_y_indices(i, j),v_z_supp_idx(i, j))
        end if
        ! get data of W
        proc_idx = w_proc(i, j)
        if ( proc_idx .eq. myid ) then
          Eu_(j + 2 * nb) = Eu_(j + 2 * nb)+ &
          w_weights(i, j) * W_(w_x_indices(i, j),w_y_indices(i, j),w_z_local_indices(i, j))
        else
          Eu_(j + 2 * nb) = Eu_(j + 2 * nb)+ &
          w_weights(i, j) * W_supp_(w_x_indices(i, j),w_y_indices(i, j),w_z_supp_idx(i, j))
        end if
      End Do
    End Do
  End Subroutine local_regT

  Subroutine local_regT_1n(E1nu_, U_, V_, W_, U_supp_, V_supp_, W_supp_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: U_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: V_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: W_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: U_supp_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: V_supp_(:, :, :)
    Real(Int64), Contiguous, Intent(In) :: W_supp_(:, :, :)
    Real(Int64), Contiguous, Intent(Out) :: E1nu_(:)
    integer :: proc_idx

    Integer(Int32) :: i, j
    E1nu_(nb_start : nb_end) = 0d0
    E1nu_(nb + nb_start : nb + nb_end) = 0d0
    E1nu_(2 * nb + nb_start : 2 * nb + nb_end) = 0d0
     
    Do j = nb_start, nb_end
      Do i = 1, nweights
        ! get data of U
        proc_idx = u_proc(i, j)
        if ( proc_idx .eq. myid ) then
          E1nu_(j) = E1nu_(j)+ &
          u_weights(i, j) * U_(u_x_indices(i, j),u_y_indices(i, j),u_z_local_indices(i, j)) * dxnu(i, j)
        else
          E1nu_(j) = E1nu_(j)+ &
          u_weights(i, j) * U_supp_(u_x_indices(i, j),u_y_indices(i, j),u_z_supp_idx(i, j)) * dxnu(i, j)
        end if
        ! get data of V
        proc_idx = v_proc(i, j)
        if ( proc_idx .eq. myid ) then
          E1nu_(j + nb) = E1nu_(j + nb)+ &
          v_weights(i, j) * V_(v_x_indices(i, j),v_y_indices(i, j),v_z_local_indices(i, j)) * dxnv(i, j)
        else
          E1nu_(j + nb) = E1nu_(j + nb)+ &
          v_weights(i, j) * V_supp_(v_x_indices(i, j),v_y_indices(i, j),v_z_supp_idx(i, j)) * dxnv(i, j)
        end if
        ! get data of W
        proc_idx = w_proc(i, j)
        if ( proc_idx .eq. myid ) then
          E1nu_(j + 2 * nb) = E1nu_(j + 2 * nb)+ &
          w_weights(i, j) * W_(w_x_indices(i, j),w_y_indices(i, j),w_z_local_indices(i, j)) * dxnw(i, j)
        else
          E1nu_(j + 2 * nb) = E1nu_(j + 2 * nb)+ &
          w_weights(i, j) * W_supp_(w_x_indices(i, j),w_y_indices(i, j),w_z_supp_idx(i, j)) * dxnw(i, j)
        end if
      End Do
    End Do
  End Subroutine local_regT_1n

  Subroutine local_regTc_1n(E1np_, P_, P_supp_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: P_(2:, 2:, 2:)
    Real(Int64), Contiguous, Intent(In) :: P_supp_(2:, 2:, :)
    Real(Int64), Contiguous, Intent(Out) :: E1np_(:)
    integer :: proc_idx

    Integer(Int32) :: i, j
    E1np_(nb_start : nb_end) = 0d0

    Do j = nb_start, nb_end
      Do i = 1, nweights
        ! get data of P
        proc_idx = c_proc(i, j)
        if ( proc_idx .eq. myid ) then
          E1np_(j) = E1np_(j)+ &
          c_weights(i, j) * P_(c_x_indices(i, j),c_y_indices(i, j),c_z_local_indices(i, j)) * dxnc(i, j)
        else
          E1np_(j) = E1np_(j)+ &
          c_weights(i, j) * P_supp_(c_x_indices(i, j),c_y_indices(i, j),c_z_supp_idx(i, j)) * dxnc(i, j)
        end if
      End Do
    End Do
  End Subroutine local_regTc_1n

  Subroutine local_regTc(Ec_, P_, P_supp_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: P_(2:, 2:, 2:)
    Real(Int64), Contiguous, Intent(In) :: P_supp_(2:, 2:, :)
    Real(Int64), Contiguous, Intent(Out) :: Ec_(:)
    Integer(Int32) :: i, j, proc_idx

    Ec_(nb_start : nb_end) = 0d0
    Do j = nb_start, nb_end
      Do i = 1, nweights
        proc_idx = c_proc(i, j)
        If (proc_idx == myid) Then
          Ec_(j) = Ec_(j) + c_weights(i, j) * &
            P_(c_x_indices(i, j), c_y_indices(i, j), c_z_local_indices(i, j))
        Else
          Ec_(j) = Ec_(j) + c_weights(i, j) * &
            P_supp_(c_x_indices(i, j), c_y_indices(i, j), c_z_supp_idx(i, j))
        End If
      End Do
    End Do

  End Subroutine local_regTc

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
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: F(:,:,:)    
    Real(Int64), Contiguous, Intent(Out) :: F_supp(:,:,:)   
    integer :: i, j
    integer :: proc_idx
    integer :: xi, yi, zi,zi_supp,zi_loc
    Real(Int64) :: factor, weight

    ! Initialize output to zero
    F = 0.0
    F_supp = 0.0

    ! Determine scaling factor based on id
    factor = 1.0 / (dx * dymin * dz)

    ! Loop over body points and accumulate using the appropriate weight array
    select case (id)
    case (1)
      do j = nb_start, nb_end
        do i = 1, nweights
          proc_idx = u_proc(i, j)
          xi       = u_x_indices(i, j)
          yi       = u_y_indices(i, j)
          zi_loc   = u_z_local_indices(i, j)
          zi_supp  = u_z_supp_idx(i, j)
          weight = u_weights(i,j)
          if (proc_idx == myid) then
            F(xi, yi, zi_loc) = F(xi, yi, zi_loc)+weight * factor* sb(j) *f_(j)
          else
            if (zi_supp < 1 .or. zi_supp > 2*suppz+1) Then
              WRITE(*,*) 'myid',myid,'proc_idx',proc_idx,'zi_supp',zi_supp
              stop 'Error: zi_supp out of [1..suppz]'
            END IF
            F_supp(xi, yi, zi_supp) = F_supp(xi, yi, zi_supp)+weight * factor* sb(j) *f_(j)
          end if
        end do
      end do
    case (2)
      do j = nb_start, nb_end
        do i = 1, nweights
          proc_idx = v_proc(i, j)
          xi       = v_x_indices(i, j)
          yi       = v_y_indices(i, j)
          zi_loc   = v_z_local_indices(i, j)
          zi_supp  = v_z_supp_idx(i, j)
          weight = v_weights(i,j)
          if (proc_idx == myid) then
            F(xi, yi, zi_loc) = F(xi, yi, zi_loc)+weight * factor* sb(j) *f_(j)
          else
            if (zi_supp < 1 .or. zi_supp > 2*suppz+1) Then
              WRITE(*,*) 'myid',myid,'proc_idx',proc_idx,'zi_supp',zi_supp
              stop 'Error: zi_supp out of [1..suppz]'
            END IF
            F_supp(xi, yi, zi_supp) = F_supp(xi, yi, zi_supp)+weight * factor* sb(j) *f_(j)
          end if
        end do
      end do
    case (3)
      do j = nb_start, nb_end
        do i = 1, nweights
          proc_idx = w_proc(i, j)
          xi       = w_x_indices(i, j)
          yi       = w_y_indices(i, j)
          zi_loc   = w_z_local_indices(i, j)
          zi_supp  = w_z_supp_idx(i, j)
          weight = w_weights(i,j)
          if (proc_idx == myid) then
            F(xi, yi, zi_loc) = F(xi, yi, zi_loc)+weight * factor* sb(j) *f_(j)
          else
            if (zi_supp < 1 .or. zi_supp > 2*suppz+1) Then
              WRITE(*,*) 'myid',myid,'proc_idx',proc_idx,'zi_supp',zi_supp
              stop 'Error: zi_supp out of [1..suppz]'
            END IF
            F_supp(xi, yi, zi_supp) = F_supp(xi, yi, zi_supp)+weight * factor* sb(j) *f_(j)
          end if
        end do
      end do
    case default
      print *, 'Error in scatter_IB_subset: invalid id =', id
      stop
    end select
  end subroutine local_reg

  subroutine local_reg_1n(F, F_supp, f_, id)
    implicit none
    integer, intent(in)  :: id
    Real(Int64), Contiguous, Intent(In) :: f_(:)
    Real(Int64), Contiguous, Intent(Out) :: F(:,:,:)    
    Real(Int64), Contiguous, Intent(Out) :: F_supp(:,:,:)   
    integer :: i, j
    integer :: proc_idx
    integer :: xi, yi, zi,zi_supp,zi_loc
    Real(Int64) :: factor, weight

    ! Initialize output to zero
    F = 0.0
    F_supp = 0.0

    ! Determine scaling factor based on id
    factor = 1.0 / (dx * dymin * dz)

    ! Loop over body points and accumulate using the appropriate weight array
    select case (id)
    case (1)
      do j = nb_start, nb_end
        do i = 1, nweights
          proc_idx = u_proc(i, j)
          xi       = u_x_indices(i, j)
          yi       = u_y_indices(i, j)
          zi_loc   = u_z_local_indices(i, j)
          zi_supp  = u_z_supp_idx(i, j)
          weight = u_weights(i,j) * dxnu(i,j)
          if (proc_idx == myid) then
            F(xi, yi, zi_loc) = F(xi, yi, zi_loc)+weight * factor* sb(j) *f_(j)
          else
            if (zi_supp < 1 .or. zi_supp > 2*suppz+1) Then
              WRITE(*,*) 'myid',myid,'proc_idx',proc_idx,'zi_supp',zi_supp
              stop 'Error: zi_supp out of [1..suppz]'
            END IF
            F_supp(xi, yi, zi_supp) = F_supp(xi, yi, zi_supp)+weight * factor* sb(j) *f_(j)
          end if
        end do
      end do
    case (2)
      do j = nb_start, nb_end
        do i = 1, nweights
          proc_idx = v_proc(i, j)
          xi       = v_x_indices(i, j)
          yi       = v_y_indices(i, j)
          zi_loc   = v_z_local_indices(i, j)
          zi_supp  = v_z_supp_idx(i, j)
          weight = v_weights(i,j) * dxnv(i,j)
          if (proc_idx == myid) then
            F(xi, yi, zi_loc) = F(xi, yi, zi_loc)+weight * factor* sb(j) *f_(j)
          else
            if (zi_supp < 1 .or. zi_supp > 2*suppz+1) Then
              WRITE(*,*) 'myid',myid,'proc_idx',proc_idx,'zi_supp',zi_supp
              stop 'Error: zi_supp out of [1..suppz]'
            END IF
            F_supp(xi, yi, zi_supp) = F_supp(xi, yi, zi_supp)+weight * factor* sb(j) *f_(j)
          end if
        end do
      end do
    case (3)
      do j = nb_start, nb_end
        do i = 1, nweights
          proc_idx = w_proc(i, j)
          xi       = w_x_indices(i, j)
          yi       = w_y_indices(i, j)
          zi_loc   = w_z_local_indices(i, j)
          zi_supp  = w_z_supp_idx(i, j)
          weight = w_weights(i,j) * dxnw(i,j)
          if (proc_idx == myid) then
            F(xi, yi, zi_loc) = F(xi, yi, zi_loc)+weight * factor* sb(j) *f_(j)
          else
            if (zi_supp < 1 .or. zi_supp > 2*suppz+1) Then
              WRITE(*,*) 'myid',myid,'proc_idx',proc_idx,'zi_supp',zi_supp
              stop 'Error: zi_supp out of [1..suppz]'
            END IF
            F_supp(xi, yi, zi_supp) = F_supp(xi, yi, zi_supp)+weight * factor* sb(j) *f_(j)
          end if
        end do
      end do
    case default
      print *, 'Error in scatter_IB_subset: invalid id =', id
      stop
    end select
  end subroutine local_reg_1n

  Subroutine apply_momentum_forcing(Fx_, Fy_, Fz_, dudn_jump_, p_jump_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: dudn_jump_(:)
    Real(Int64), Contiguous, Intent(In) :: p_jump_(:)
    Real(Int64), Contiguous, Intent(Out) :: Fx_(:, :, :)
    Real(Int64), Contiguous, Intent(Out) :: Fy_(:, :, :)
    Real(Int64), Contiguous, Intent(Out) :: Fz_(:, :, :)
    Real(Int64) :: fwork(nb)
    Integer(Int32) :: i, j, k

    Fx_ = 0d0
    Fy_ = 0d0
    Fz_ = 0d0

    fwork = nu * dudn_jump_(1 : nb) * normals(1 : nb) ** 2
    Call regu(U_ib_tmp, fwork)
    fwork = nu * dudn_jump_(1 : nb) * normals(nb + 1 : 2 * nb) ** 2
    Call regv(V_ib_tmp, fwork)
    fwork = nu * dudn_jump_(1 : nb) * normals(2 * nb + 1 : 3 * nb) ** 2
    Call regw(W_ib_tmp, fwork)
    Call apply_boundary_conditions(U_ib_tmp, V_ib_tmp, W_ib_tmp)
    Do k = 2, nzg - 1
      Do j = 2, nyg - 1
        Do i = 2, nx - 1
          Fx_(i, j, k) = 0.25d0 * (U_ib_tmp(i - 1, j, k) + 2d0 * U_ib_tmp(i, j, k) + U_ib_tmp(i + 1, j, k)) &
                       + 0.25d0 * (V_ib_tmp(i, j - 1, k) + V_ib_tmp(i + 1, j - 1, k) &
                                 + V_ib_tmp(i, j, k) + V_ib_tmp(i + 1, j, k)) &
                       + 0.25d0 * (W_ib_tmp(i, j, k - 1) + W_ib_tmp(i + 1, j, k - 1) &
                                 + W_ib_tmp(i, j, k) + W_ib_tmp(i + 1, j, k))
        End Do
      End Do
    End Do

    fwork = nu * dudn_jump_(nb + 1 : 2 * nb) * normals(1 : nb) ** 2
    Call regu(U_ib_tmp, fwork)
    fwork = nu * dudn_jump_(nb + 1 : 2 * nb) * normals(nb + 1 : 2 * nb) ** 2
    Call regv(V_ib_tmp, fwork)
    fwork = nu * dudn_jump_(nb + 1 : 2 * nb) * normals(2 * nb + 1 : 3 * nb) ** 2
    Call regw(W_ib_tmp, fwork)
    Call apply_boundary_conditions(U_ib_tmp, V_ib_tmp, W_ib_tmp)
    Do k = 2, nzg - 1
      Do j = 2, ny - 1
        Do i = 2, nxg - 1
          Fy_(i, j, k) = 0.25d0 * (U_ib_tmp(i - 1, j, k) + U_ib_tmp(i - 1, j + 1, k) &
                                 + U_ib_tmp(i, j, k) + U_ib_tmp(i, j + 1, k)) &
                       + 0.25d0 * (V_ib_tmp(i, j - 1, k) + 2d0 * V_ib_tmp(i, j, k) + V_ib_tmp(i, j + 1, k)) &
                       + 0.25d0 * (W_ib_tmp(i, j, k - 1) + W_ib_tmp(i, j + 1, k - 1) &
                                 + W_ib_tmp(i, j, k) + W_ib_tmp(i, j + 1, k))
        End Do
      End Do
    End Do

    fwork = nu * dudn_jump_(2 * nb + 1 : 3 * nb) * normals(1 : nb) ** 2
    Call regu(U_ib_tmp, fwork)
    fwork = nu * dudn_jump_(2 * nb + 1 : 3 * nb) * normals(nb + 1 : 2 * nb) ** 2
    Call regv(V_ib_tmp, fwork)
    fwork = nu * dudn_jump_(2 * nb + 1 : 3 * nb) * normals(2 * nb + 1 : 3 * nb) ** 2
    Call regw(W_ib_tmp, fwork)
    Call apply_boundary_conditions(U_ib_tmp, V_ib_tmp, W_ib_tmp)
    Do k = 2, nz - 1
      Do j = 2, nyg - 1
        Do i = 2, nxg - 1
          Fz_(i, j, k) = 0.25d0 * (U_ib_tmp(i - 1, j, k) + U_ib_tmp(i - 1, j, k + 1) &
                                 + U_ib_tmp(i, j, k) + U_ib_tmp(i, j, k + 1)) &
                       + 0.25d0 * (V_ib_tmp(i, j - 1, k) + V_ib_tmp(i, j - 1, k + 1) &
                                 + V_ib_tmp(i, j, k) + V_ib_tmp(i, j, k + 1)) &
                       + 0.25d0 * (W_ib_tmp(i, j, k - 1) + 2d0 * W_ib_tmp(i, j, k) + W_ib_tmp(i, j, k + 1))
        End Do
      End Do
    End Do

    Call apply_boundary_conditions(Fx_, Fy_, Fz_)
    Call accumulate_DR1n_forcing(Fx_, Fy_, Fz_, dudn_jump_)
    Call apply_boundary_conditions(Fx_, Fy_, Fz_)

    fwork = -p_jump_ * normals(1 : nb)
    Call regu(U_ib_tmp, fwork)
    fwork = -p_jump_ * normals(nb + 1 : 2 * nb)
    Call regv(V_ib_tmp, fwork)
    fwork = -p_jump_ * normals(2 * nb + 1 : 3 * nb)
    Call regw(W_ib_tmp, fwork)
    Fx_ = Fx_ + U_ib_tmp
    Fy_ = Fy_ + V_ib_tmp
    Fz_ = Fz_ + W_ib_tmp
    Call apply_boundary_conditions(Fx_, Fy_, Fz_)

  End Subroutine apply_momentum_forcing

  Subroutine apply_continuity_forcing(F_, dudn_jump_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: dudn_jump_(:)
    Real(Int64), Contiguous, Intent(Out) :: F_(2:, 2:, 2:)
    Real(Int64) :: fwork(nb)
    Integer(Int32) :: i, j, k

    fwork = dudn_jump_(1 : nb) * normals(1 : nb)
    Call regu_1n(U_ib_tmp, fwork)
    fwork = dudn_jump_(nb + 1 : 2 * nb) * normals(nb + 1 : 2 * nb)
    Call regv_1n(V_ib_tmp, fwork)
    fwork = dudn_jump_(2 * nb + 1 : 3 * nb) * normals(2 * nb + 1 : 3 * nb)
    Call regw_1n(W_ib_tmp, fwork)
    Call apply_boundary_conditions(U_ib_tmp, V_ib_tmp, W_ib_tmp)

    F_ = 0d0
    Do k = 2, nzg - 1
      Do j = 2, nyg - 1
        Do i = 2, nxg - 1
          F_(i, j, k) = 0.5d0 * (U_ib_tmp(i - 1, j, k) + U_ib_tmp(i, j, k)) &
                       + 0.5d0 * (V_ib_tmp(i, j - 1, k) + V_ib_tmp(i, j, k)) &
                       + 0.5d0 * (W_ib_tmp(i, j, k - 1) + W_ib_tmp(i, j, k))
        End Do
      End Do
    End Do

  End Subroutine apply_continuity_forcing

  Integer(Int32) Function wrap_x_ib(i)
    Implicit None
    Integer(Int32), Intent(In) :: i

    wrap_x_ib = Modulo(i - 2, nx_global - 2) + 2

  End Function wrap_x_ib

  Integer(Int32) Function wrap_z_ib(k)
    Implicit None
    Integer(Int32), Intent(In) :: k

    wrap_z_ib = Modulo(k - 2, nz_global - 2) + 2

  End Function wrap_z_ib

  Subroutine accumulate_DR1n_forcing(Fx_, Fy_, Fz_, dudn_jump_)
    Implicit None
    Real(Int64), Contiguous, Intent(In) :: dudn_jump_(:)
    Real(Int64), Contiguous, Intent(InOut) :: Fx_(:, :, :)
    Real(Int64), Contiguous, Intent(InOut) :: Fy_(:, :, :)
    Real(Int64), Contiguous, Intent(InOut) :: Fz_(:, :, :)
    Integer(Int32), Parameter :: nlx = 2 * suppx + 4
    Integer(Int32), Parameter :: nly = 2 * suppy + 4
    Integer(Int32), Parameter :: nlz = 2 * suppz + 4
    Real(Int64) :: Ru(nlx, nly, nlz), Rv(nlx, nly, nlz), Rw(nlx, nly, nlz)
    Real(Int64) :: Du(nlx, nly, nlz), Dv(nlx, nly, nlz), Dw(nlx, nly, nlz)
    Real(Int64) :: T1u(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Real(Int64) :: T2u(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Real(Int64) :: T3u(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Real(Int64) :: T1v(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Real(Int64) :: T2v(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Real(Int64) :: T3v(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Real(Int64) :: T1w(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Real(Int64) :: T2w(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Real(Int64) :: T3w(0:nlx + 1, 0:nly + 1, 0:nlz + 1)
    Integer(Int32) :: ig(nlx), kg(nlz)
    Integer(Int32) :: marker, l, m, n, count, p, q, r, pn, qn, rn, pm, rm
    Integer(Int32) :: i0, j0, k0, dxoff, dyoff, dzoff, Lx, Lz, pmax, qmax, rmax
    Integer(Int32) :: pu, qu, su, pv, qv, sv, pw, qw, sw, ii, jj, kk, kl
    Logical(Int32) :: marker_is_local
    Real(Int64) :: dvol, sbk, nrmx, nrmy, nrmz, bu, bv, bw
    Real(Int64) :: IxRu, IyRu, IzRu, IxRv, IyRv, IzRv, IxRw, IyRw, IzRw
    Real(Int64) :: IxDu, IyDu, IzDu, IxDv, IyDv, IzDv, IxDw, IyDw, IzDw

    dvol = dx * dymin * dz
    Lx = Min(nlx, nx_global - 2)
    Lz = Min(nlz, nz_global - 2)

    Do marker = 1, nb
      dxoff = xm_pivot_index(marker) + 1 - x_pivot_index(marker)
      dyoff = ym_pivot_index(marker) + 1 - y_pivot_index(marker)
      dzoff = zm_pivot_index(marker) - z_pivot_index(marker)
      i0 = x_pivot_index(marker) - suppx - 2
      j0 = y_pivot_index(marker) - suppy - 2
      k0 = z_pivot_index(marker) - suppz - 2

      Do p = 1, Lx
        ig(p) = wrap_x_ib(i0 + p)
      End Do
      Do r = 1, Lz
        kg(r) = wrap_z_ib(k0 + r)
      End Do

      rmax = Min(2 * suppz + 2 + dzoff, Lz)
      marker_is_local = .False.
      Do r = 1, Min(rmax + 1, Lz)
        kk = kg(r)
        If ((kk >= kg1_global(myid) + 1 .And. kk <= kg2_global(myid) - 1) .Or. &
            (kk >= k1_global(myid) + 1 .And. kk <= k2_global(myid) - 1)) Then
          marker_is_local = .True.
          Exit
        End If
      End Do
      If (.Not. marker_is_local) Cycle

      sbk = sb(marker)
      nrmx = normals(marker)
      nrmy = normals(nb + marker)
      nrmz = normals(2 * nb + marker)
      bu = nu * dudn_jump_(marker)
      bv = nu * dudn_jump_(nb + marker)
      bw = nu * dudn_jump_(2 * nb + marker)

      Ru = 0d0
      Rv = 0d0
      Rw = 0d0
      Du = 0d0
      Dv = 0d0
      Dw = 0d0
      count = 0
      Do l = -suppz, suppz
        Do m = -suppy, suppy
          Do n = -suppx, suppx
            count = count + 1
            pu = Modulo(n + suppx + 1, Lx) + 1
            pv = Modulo(n + suppx + 1 + dxoff, Lx) + 1
            pw = pv
            qu = m + suppy + 2 + dyoff
            qv = m + suppy + 2
            qw = qu
            su = Modulo(l + suppz + 1 + dzoff, Lz) + 1
            sv = su
            sw = Modulo(l + suppz + 1, Lz) + 1
            Ru(pu, qu, su) = Ru(pu, qu, su) + u_weights(count, marker) * sbk * nrmx / dvol
            Rv(pv, qv, sv) = Rv(pv, qv, sv) + v_weights(count, marker) * sbk * nrmy / dvol
            Rw(pw, qw, sw) = Rw(pw, qw, sw) + w_weights(count, marker) * sbk * nrmz / dvol
            Du(pu, qu, su) = bu * dxnu(count, marker)
            Dv(pv, qv, sv) = bv * dxnv(count, marker)
            Dw(pw, qw, sw) = bw * dxnw(count, marker)
          End Do
        End Do
      End Do

      pmax = Min(2 * suppx + 2 + dxoff, Lx)
      qmax = 2 * suppy + 2 + dyoff
      T1u = 0d0
      T2u = 0d0
      T3u = 0d0
      T1v = 0d0
      T2v = 0d0
      T3v = 0d0
      T1w = 0d0
      T2w = 0d0
      T3w = 0d0

      Do r = 1, rmax
        rn = Modulo(r, Lz) + 1
        Do q = 1, qmax
          qn = q + 1
          Do p = 1, pmax
            pn = Modulo(p, Lx) + 1
            IxRu = 0.5d0 * (Ru(p, q, r) + Ru(pn, q, r))
            IyRu = 0.5d0 * (Ru(p, q, r) + Ru(p, qn, r))
            IzRu = 0.5d0 * (Ru(p, q, r) + Ru(p, q, rn))
            IxRv = 0.5d0 * (Rv(p, q, r) + Rv(pn, q, r))
            IyRv = 0.5d0 * (Rv(p, q, r) + Rv(p, qn, r))
            IzRv = 0.5d0 * (Rv(p, q, r) + Rv(p, q, rn))
            IxRw = 0.5d0 * (Rw(p, q, r) + Rw(pn, q, r))
            IyRw = 0.5d0 * (Rw(p, q, r) + Rw(p, qn, r))
            IzRw = 0.5d0 * (Rw(p, q, r) + Rw(p, q, rn))
            IxDu = 0.5d0 * (Du(p, q, r) + Du(pn, q, r))
            IyDu = 0.5d0 * (Du(p, q, r) + Du(p, qn, r))
            IzDu = 0.5d0 * (Du(p, q, r) + Du(p, q, rn))
            IxDv = 0.5d0 * (Dv(p, q, r) + Dv(pn, q, r))
            IyDv = 0.5d0 * (Dv(p, q, r) + Dv(p, qn, r))
            IzDv = 0.5d0 * (Dv(p, q, r) + Dv(p, q, rn))
            IxDw = 0.5d0 * (Dw(p, q, r) + Dw(pn, q, r))
            IyDw = 0.5d0 * (Dw(p, q, r) + Dw(p, qn, r))
            IzDw = 0.5d0 * (Dw(p, q, r) + Dw(p, q, rn))
            T1u(p, q, r) = IxRu * IxDu
            T2u(p, q, r) = IxRv * IyDu
            T3u(p, q, r) = IxRw * IzDu
            T1v(p, q, r) = IyRu * IxDv
            T2v(p, q, r) = IyRv * IyDv
            T3v(p, q, r) = IyRw * IzDv
            T1w(p, q, r) = IzRu * IxDw
            T2w(p, q, r) = IzRv * IyDw
            T3w(p, q, r) = IzRw * IzDw
          End Do
        End Do
      End Do

      Do r = 1, Min(rmax + 1, Lz)
        rm = Modulo(r - 2, Lz) + 1
        kk = kg(r)
        Do q = 1, qmax + 1
          jj = j0 + q
          Do p = 1, Min(pmax + 1, Lx)
            pm = Modulo(p - 2, Lx) + 1
            ii = ig(p)
            If (kk >= kg1_global(myid) + 1 .And. kk <= kg2_global(myid) - 1) Then
              kl = kk - kg1_global(myid) + 1
              If (jj >= 2 .And. jj <= nyg - 1) Then
                Fx_(ii, jj, kl) = Fx_(ii, jj, kl) &
                  + (T1u(p, q, r) - T1u(pm, q, r)) / (xg_global(ii + 1) - xg_global(ii)) &
                  + (T2u(p, q, r) - T2u(p, q - 1, r)) / (y_global(jj) - y_global(jj - 1)) &
                  + (T3u(p, q, r) - T3u(p, q, rm)) / (z_global(kk) - z_global(kk - 1))
              End If
              If (jj >= 2 .And. jj <= ny - 1) Then
                Fy_(ii, jj, kl) = Fy_(ii, jj, kl) &
                  + (T1v(p, q, r) - T1v(pm, q, r)) / (x_global(ii) - x_global(ii - 1)) &
                  + (T2v(p, q, r) - T2v(p, q - 1, r)) / (yg_global(jj + 1) - yg_global(jj)) &
                  + (T3v(p, q, r) - T3v(p, q, rm)) / (z_global(kk) - z_global(kk - 1))
              End If
            End If
            If (kk >= k1_global(myid) + 1 .And. kk <= k2_global(myid) - 1) Then
              kl = kk - k1_global(myid) + 1
              If (jj >= 2 .And. jj <= nyg - 1) Then
                Fz_(ii, jj, kl) = Fz_(ii, jj, kl) &
                  + (T1w(p, q, r) - T1w(pm, q, r)) / (x_global(ii) - x_global(ii - 1)) &
                  + (T2w(p, q, r) - T2w(p, q - 1, r)) / (y_global(jj) - y_global(jj - 1)) &
                  + (T3w(p, q, r) - T3w(p, q, rm)) / (zg_global(kk + 1) - zg_global(kk))
              End If
            End If
          End Do
        End Do
      End Do
    End Do

  End Subroutine accumulate_DR1n_forcing

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

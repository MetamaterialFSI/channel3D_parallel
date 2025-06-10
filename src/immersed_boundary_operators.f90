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
    character(len=50) :: filename
    integer :: k_local,k_supp         ! local index
    integer :: proc_id      ! MPI rank that owns this index

    !-----------------Find the indices of the flow grid points closest to the body points---------------------!
    Do l = nb_start, nb_end
      ! index of largest x value that is smaller than xb(l)
      x_pivot_index(l) = Int(xb(l) / dx) + 1
      ! index of largest xm value that is smaller than xb(l), accounting for periodicity, but not for ghost cells
      xm_pivot_index(l) = Modulo(Int((xb(l) + 0.5 * dx) / dx) - 1, nxm_global - 1) + 1
      ! index of largest y value that is smaller than yb(l). The use of reference points is more efficient here than searching
      ! the y array for the closest point, in case there is grid stretching.
      y_pivot_index(l) = Int((yb(l) - y(y_ref_index(l))) / dymin) + y_ref_index(l)
      ! index of largest ym value that is smaller than yb(l), not accounting for ghost cells
      ym_pivot_index(l) = Int((yb(l) - (y(y_ref_index(l)) + 0.5 * dymin)) / dymin) + y_ref_index(l)
      ! index of largest z value that is smaller than zb(l)
      z_pivot_index(l) = Int(zb(l) / dz) + 1
      ! index of largest zm value that is smaller than zb(l), accounting for periodicity, but not for ghost cells
      zm_pivot_index(l) = Modulo(Int((zb(l) + 0.5 * dz) / dz) - 1, nzm_global - 1) + 1 
    End Do
    ! debug line
    ! Do l=2,nzg_global-1
    !   call global_to_local_center(l, k_supp, k_local, proc_id)
    !   WRITE(*,*) 'Center: myid',myid,'proc_id',proc_id,'k_global',l,'k_local',k_local,'supp_idx',k_supp        
    ! end do
    ! Do l=2,nz_global-1
    !   call global_to_local_face(l, k_supp, k_local, proc_id)
    !   WRITE(*,*) 'Face: myid',myid,'proc_id',proc_id,'k_global',l,'k_local',k_local,'supp_idx',k_supp        
    ! end do

    !----Find the flow grid indices and corresponding DDF weights within the support of the DDF at each body point----!
    ! u velocity locations
    Do l = nb_start, nb_end
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
            if ( moving_z_flag ) then
              call global_to_local_center(kk + 1, k_supp, k_local, proc_id)
              u_z_local_indices(count, l) = k_local
              u_z_supp_idx(count, l) = k_supp
              u_proc(count, l) = proc_id 
            Else
              if ( istep .le. 1 ) then
                call global_to_local_center(kk + 1, k_supp, k_local, proc_id)
                u_z_local_indices(count, l) = k_local
                u_z_supp_idx(count, l) = k_supp
                u_proc(count, l) = proc_id
              end if
            end if 
            ! if ( myid .ne. proc_id ) then
            !   if ( k_supp.le.0 .or. k_supp .gt. 2*(suppz)+1  ) then
            !     WRITE(*,*) 'U: myid',myid,'proc_id',proc_id,'k_global',kk+1,'k_local',k_local,'supp_idx',k_supp
            !   end if 
            ! end if
            

            u_weights(count, l) =  dx * dymin * dz &
              * deltafnc( x_global(ii), xb(l), dx,    Lxp) &
              * deltafnc(ym_global(jj), yb(l), dymin, 0.d0) &
              * deltafnc(zm_global(kk), zb(l), dz,    Lzp)
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
            ii = Modulo(xm_pivot_index(l) + i - 1, nxm_global - 1) + 1
            jj = y_pivot_index(l) + j
            kk = Modulo(zm_pivot_index(l) + k - 1, nzm_global - 1) + 1

            count = count+1          
            v_x_indices(count, l) = ii + 1 ! plus one for ghost cell
            v_y_indices(count, l) = jj
            v_z_indices(count, l) = kk + 1 ! plus one for ghost cell
            if ( moving_z_flag ) then
              call global_to_local_center(kk + 1, k_supp, k_local, proc_id)
              v_z_local_indices(count, l) = k_local
              v_z_supp_idx(count, l) = k_supp
              v_proc(count, l) = proc_id 
            Else
              if ( istep .le. 1 ) then
                call global_to_local_center(kk + 1, k_supp, k_local, proc_id)
                v_z_local_indices(count, l) = k_local
                v_z_supp_idx(count, l) = k_supp
                v_proc(count, l) = proc_id
              end if
            end if 

            v_weights(count, l) = dx * dymin * dz &
              * deltafnc(xm_global(ii), xb(l), dx,    Lxp) &
              * deltafnc( y_global(jj), yb(l), dymin, 0.d0) &
              * deltafnc(zm_global(kk), zb(l), dz,    Lzp)
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
            ii = Modulo(xm_pivot_index(l) + i - 1, nxm_global - 1) + 1
            jj = ym_pivot_index(l) + j
            kk = Modulo(z_pivot_index(l) + k - 2, nz_global - 2 ) + 2

            count = count + 1           
            w_x_indices(count, l) = ii + 1 ! plus one for ghost cell
            w_y_indices(count, l) = jj + 1 ! plus one for ghost cell
            w_z_indices(count, l) = kk
            if ( moving_z_flag ) then
              call global_to_local_face(kk, k_supp, k_local, proc_id)
              w_z_local_indices(count, l) = k_local
              w_z_supp_idx(count, l) = k_supp
              w_proc(count, l) = proc_id 
            Else
              if ( istep .le. 1 ) then
                call global_to_local_face(kk, k_supp, k_local, proc_id)
                w_z_local_indices(count, l) = k_local
                w_z_supp_idx(count, l) = k_supp
                w_proc(count, l) = proc_id
              end if
            end if 
            ! if ( myid .ne. proc_id ) then
            !   if ( k_supp.le.0 .or. k_supp .gt. 2*(suppz)+1 ) then
            !     WRITE(*,*) 'W: myid',myid,'proc_id',proc_id,'k_global',kk+1,'k_local',k_local,'supp_idx',k_supp
            !   end if 
            ! end if

            w_weights(count, l) = dx * dymin * dz &
              * deltafnc(xm_global(ii), xb(l), dx,    Lxp) &
              * deltafnc(ym_global(jj), yb(l), dymin, 0.d0) &
              * deltafnc( z_global(kk), zb(l), dz,    Lzp)
          End Do
        End Do
      End Do
    End Do

    !local_size_nb = nb_end - nb_start + 1

    ! ! Gather send_counts
    ! Call MPI_Allgather(local_size_nb, 1, MPI_INT, send_counts_nb, 1, MPI_INT, MPI_COMM_WORLD, ierr)

    ! displs_nb(1) = 0
    ! Do i = 2, nprocs
    !   displs_nb(i) = displs_nb(i-1) + send_counts_nb(i-1)
    ! End Do

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

  Function regT(U_, V_, W_)
    Implicit None
    ! Real(Int64), Dimension(nx, nyg, nzg), Intent(In) :: U_
    ! Real(Int64), Dimension(nxg, ny, nzg), Intent(In) :: V_
    ! Real(Int64), Dimension(nxg, nyg, nz), Intent(In) :: W_
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: U_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: V_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: W_(:, :, :)
    Real(Int64), Dimension(3 * nb):: regT
    Real(Int64), Dimension(3 * nb):: regT_global

    ! Gather U, V, W fields into global fields
    ! Call MPI_Allgatherv(U_(:, :, 2:nzm+1), local_size_U, MPI_REAL8, &
    !              U_global(:, :, 2:nzm_global+1), send_counts_U, displs_U, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ! Call MPI_AllGatherv(V_(:, :, 2:nzm+1), local_size_V, MPI_REAL8, &
    !              V_global(:, :, 2:nzm_global+1), send_counts_V, displs_V, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ! Call MPI_AllGatherv(W_(:, :, 2:nz-1), local_size_W, MPI_REAL8, &
    !              W_global(:, :, 2:nz_global-1), send_counts_W, displs_W, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ! !Gather U, V, W fields into the subset field
    ! update the support vector for U, V,W
    !WRITE(*,*) 'myid',myid,'update the support vector'
    prev_internal = MPI_WTIME()
    Call interior_planes_update_support(U_,U_supp,1)
    Call interior_planes_update_support(V_,V_supp,2)
    Call interior_planes_update_support(W_,W_supp,3)
    last_internal = MPI_WTIME()
    IF (E_internal_flag) then
      E_transfer = E_transfer +last_internal-prev_internal
    end if
    ! update the subset matrix
    prev_internal = MPI_WTIME()
    !WRITE(*,*) 'myid',myid,'update the subset matrix U'
    Call update_IB_subset(U_,U_supp,U_subset,1)
    !WRITE(*,*) 'myid',myid,'update the subset matrix V'
    Call update_IB_subset(V_,V_supp,V_subset,2)
    !WRITE(*,*) 'myid',myid,'update the subset matrix W'
    Call update_IB_subset(W_,W_supp,W_subset,3)
    last_internal = MPI_WTIME()
    IF (E_internal_flag) then
      E_update_subset = E_update_subset +last_internal-prev_internal
    end if

    !regT_global = global_regT(U_global, V_global, W_global)
    !WRITE(*,*) 'myid',myid,'compute regT'
    prev_internal = MPI_WTIME()
    regT = subset_regT(U_subset, V_subset, W_subset)
    last_internal = MPI_WTIME()
    IF (E_internal_flag) then
      E_subset = E_subset +last_internal-prev_internal
    end if



    ! Gather regT values from all partitions
    ! Call MPI_Allgatherv(aux_surface_vector(nb_start : nb_end), local_size_nb, MPI_REAL8, &
    !              regT_global(1 : nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ! Call MPI_Allgatherv(aux_surface_vector(nb + nb_start : nb + nb_end), local_size_nb, MPI_REAL8, &
    !              regT_global(nb + 1 : 2 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ! Call MPI_Allgatherv(aux_surface_vector(2 * nb + nb_start : 2 * nb + nb_end), local_size_nb, MPI_REAL8, &
    !              regT_global(2 * nb + 1 : 3 * nb), send_counts_nb, displs_nb, MPI_REAL8, MPI_COMM_WORLD, ierr)
    ! if ( NORM2(regT-regT_global) .gt. 1E-8 ) then
    !   if ( myid .eq. 0 ) then
    !     WRITE(*,*) 'global_regT not equal to subset_regT'
    !     WRITE(*,*) 'global_regT',regT_global
    !     WRITE(*,*) 'subset_regT',regT
    !   end if
    ! end if
  End Function regT

  Subroutine regu(U_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nx, nyg, nzg), Intent(Out) :: U_
    Real(Int64), Dimension(nx, nyg, nzg):: U_local

    prev_internal = MPI_WTIME()
    Call subset_reg(U_,U_supp,f_,1)
    last_internal = MPI_WTIME()
    IF (R_internal_flag) then
      R_subset = R_subset +last_internal-prev_internal
    end if
    
    ! Call global_regu(U_global, f_)

    ! If (myid == 0) Then
    !   Call MPI_Reduce(MPI_IN_PLACE, U_global, nx_global * nyg_global * nzg_global, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! Else
    !   Call MPI_Reduce(U_global, U_global, nx_global * nyg_global * nzg_global, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! End If

    ! Call MPI_Scatterv(U_global(:, :, 2:nzm_global+1), send_counts_U, displs_U, MPI_REAL8, &
    !               U_local(:, :, 2:nzm+1), local_size_U, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
    !Call scatter_IB_subset(U_, U_supp, U_subset, 1)
    prev_internal = MPI_WTIME()
    Call support_update_interior_planes(U_,U_supp,1)
    last_internal = MPI_WTIME()
    IF (R_internal_flag) then
      R_transfer = R_transfer +last_internal-prev_internal
    end if

  End Subroutine regu

  Subroutine regv(V_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nxg, ny, nzg), Intent(Out) :: V_

    Call subset_reg(V_,V_supp,f_,2)
    !Call scatter_IB_subset(V_, V_supp, V_subset, 2)
    Call support_update_interior_planes(V_,V_supp,2)
    !Call global_regv(V_global, f_)

    ! If (myid == 0) Then
    !   Call MPI_Reduce(MPI_IN_PLACE, V_global, nxg_global * ny_global * nzg_global, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! Else
    !   Call MPI_Reduce(V_global, V_global, nxg_global * ny_global * nzg_global, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! End If

    ! Call MPI_Scatterv(V_global(:, :, 2:nzm_global+1), send_counts_V, displs_V, MPI_REAL8, &
    !               V_(:, :, 2:nzm+1), local_size_V, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

  End Subroutine regv

  Subroutine regw(W_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nxg, nyg, nz), Intent(Out) :: W_
    Integer(Int32) :: i, j

    Call subset_reg(W_,W_supp,f_,3)
    !Call scatter_IB_subset(W_, W_supp, W_subset, 3)
    Call support_update_interior_planes(W_,W_supp,3)
    !Call global_regw(W_global, f_)

    ! If (myid == 0) Then
    !   Call MPI_Reduce(MPI_IN_PLACE, W_global, nxg_global * nyg_global * nz_global, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! Else
    !   Call MPI_Reduce(W_global, W_global, nxg_global * nyg_global * nz_global, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    ! End If

    ! Call MPI_Scatterv(W_global(:, :, 2:nz_global-1), send_counts_W, displs_W, MPI_REAL8, &
    !               W_(:, :, 2:nz-1), local_size_W, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  
  End Subroutine regw

  !-----------------------------------------------!
  !          Interpolation to body points         !
  !                                               !
  ! Input: U_, V_, W_                             !
  ! Output: regT                                  !
  !                                               !
  !-----------------------------------------------!
  Function global_regT(U_, V_, W_)
    Implicit None
    Real(Int64), Dimension(nx_global, nyg_global, nzg_global), Intent(In) :: U_
    Real(Int64), Dimension(nxg_global, ny_global, nzg_global), Intent(In) :: V_
    Real(Int64), Dimension(nxg_global, nyg_global, nz_global), Intent(In) :: W_
    Real(Int64), Dimension(3 * nb):: global_regT

    Integer(Int32) :: i, j
    global_regT = 0.D0
     
    Do j = nb_start, nb_end
      Do i = 1, nweights
        global_regT(j         ) = global_regT(j)          &
          + u_weights(i, j) * U_(u_x_indices(i, j), u_y_indices(i, j), u_z_indices(i, j))
        global_regT(j + nb    ) = global_regT(j + nb)     &
          + v_weights(i, j) * V_(v_x_indices(i, j), v_y_indices(i, j), v_z_indices(i, j))
        global_regT(j + 2 * nb) = global_regT(j + 2 * nb) &
          + w_weights(i, j) * W_(w_x_indices(i, j), w_y_indices(i, j), w_z_indices(i, j))
      End Do
    End Do

  End Function global_regT

  !-----------------------------------------------!
  !          Interpolation to body points         !
  !                                               !
  ! Input: U_, V_, W_                             !
  ! Output: regT                                  !
  !                                               !
  !-----------------------------------------------!
  Function subset_regT(U_, V_, W_)
    Implicit None
    Real(Int64), Dimension(nweights, nb), Intent(In) :: U_
    Real(Int64), Dimension(nweights, nb), Intent(In) :: V_
    Real(Int64), Dimension(nweights, nb), Intent(In) :: W_
    Real(Int64), Dimension(3 * nb):: subset_regT

    Integer(Int32) :: i, j
    subset_regT = 0.D0
     
    Do j = nb_start, nb_end
      Do i = 1, nweights
        subset_regT(j         ) = subset_regT(j)          &
          + u_weights(i, j) * U_(i,j)
          subset_regT(j + nb    ) = subset_regT(j + nb)     &
          + v_weights(i, j) * V_(i,j)
          subset_regT(j + 2 * nb) = subset_regT(j + 2 * nb) &
          + w_weights(i, j) * W_(i,j)
      End Do
    End Do

  End Function subset_regT
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
  subroutine subset_reg(F,F_supp, f_,id)
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
  
  end subroutine subset_reg

  !-----------------------------------------------!
  !      Regularization of the x-component of     ! 
  !      f_ to U grid points                      !
  !                                               !
  ! Input: f_                                     !
  ! Output: global_regu                           !
  !                                               !
  !-----------------------------------------------!
  Subroutine global_regu(U_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nx_global, nyg_global, nzg_global), Intent(Out) :: U_
    Integer(Int32) :: i, j

    U_ = 0.D0

    Do j = nb_start, nb_end
      Do i = 1, nweights
        U_(u_x_indices(i, j), u_y_indices(i, j), u_z_indices(i, j)) = &
          U_(u_x_indices(i, j), u_y_indices(i, j), u_z_indices(i, j)) &
           + u_weights(i, j) * 1.d0 / (dx * dymin * dz) * sb(j) * f_(j)
      End Do
    End Do

  End Subroutine global_regu

  !-----------------------------------------------!
  !      Regularization of the y-component of     ! 
  !      f_ to V grid points                      !
  !                                               !
  ! Input: f_                                     !
  ! Output: global_regv                           !
  !                                               !
  !-----------------------------------------------!
  Subroutine global_regv(V_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nxg_global, ny_global, nzg_global), Intent(Out) :: V_
    Integer(Int32) :: i, j

    V_ = 0.D0

    Do j = nb_start, nb_end
      Do i = 1, nweights
        V_(v_x_indices(i, j), v_y_indices(i, j), v_z_indices(i, j)) = &
          V_(v_x_indices(i, j), v_y_indices(i, j), v_z_indices(i, j)) &
           + v_weights(i, j) * 1.d0 / (dx * dymin * dz) * sb(j) * f_(j + nb)
      End Do
    End Do

  End Subroutine global_regv
  
  !-----------------------------------------------!
  !      Regularization of the z-component of     ! 
  !      f_ to W grid points                      !
  !                                               !
  ! Input: f_                                     !
  ! Output: global_regw                           !
  !                                               !
  !-----------------------------------------------!
  Subroutine global_regw(W_, f_)
    Implicit None
    Real(Int64), Dimension(3 * nb), Intent(In) :: f_
    Real(Int64), Dimension(nxg_global, nyg_global, nz_global), Intent(Out) :: W_
    Integer(Int32) :: i, j

    W_ = 0.D0

    Do j = nb_start, nb_end
      Do i = 1, nweights
        W_(w_x_indices(i, j), w_y_indices(i, j), w_z_indices(i, j)) = &
          W_(w_x_indices(i, j), w_y_indices(i, j), w_z_indices(i, j)) &
           + w_weights(i, j) * 1.d0 / (dx * dymin * dz) * sb(j) * f_(j + 2 * nb)
      End Do
    End Do
  
  End Subroutine global_regw


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
    if (k_global >= k1_global(prev)+1 .and. k_global <= k2_global(prev)-1) then
      rank = prev
    elseif (k_global >= k1_global(myid)+1 .and. k_global <= k2_global(myid)-1) then
      rank = myid
    elseif (k_global >= k1_global(next)+1 .and. k_global <= k2_global(next)-1) then
      rank = next
    else
      print *, 'Error: Face index ', k_global, ' not in {', prev, ',', myid, ',', next, '}.'
      stop
    end if
  
    if (rank == myid) then
      ! Locally owned -> compute k_loc
      k_loc = k_global - k1_global(myid) + 1
      k_sup = -1
    else  ! rank == next
      k_loc = k_global - k1_global(rank) + 1
      if ( k_loc>nz/2 ) then
        k_sup = k_global - (k2_global(prev) - suppz)+2
      else
        k_sup = (k_global - k1_global(next)) + (suppz+1)
      end if
    
    end if
  end subroutine global_to_local_face
  
  
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
    elseif (k_global >= (kg1_global(next+1)+1) .and. &
            k_global <= (kg2_global(next+1)-1)) then
      rank = next
    else
      print *, 'Error: Center index ', k_global, ' not in {', prev, ',', myid, ',', next, '}.'
      stop
    end if
  
    if (rank == myid) then
      ! Locally owned -> compute k_loc
      k_loc = k_global - kg1_global(myid) + 1
      k_sup = -1
    else  ! rank == next
      k_loc = k_global - kg1_global(rank) + 1
      if ( k_loc>nzg/2 ) then
        neighbor_top_start = (kg2_global(prev) - suppz - 2)
        k_sup = k_global - neighbor_top_start
        if ( rank .eq. nprocs-1 ) then
          k_sup=k_sup+1
        end if
      else
        k_sup = (k_global - (kg1_global(next) + 1)) + (suppz + 2)
      end if
    
    end if
  end subroutine global_to_local_center
  
  
!-------------------------------------------------------------------------------
!  update_IB_subset
!
!  Inputs:
!    F(:, :, :)     : full 3D field (real(8)) for U, V, or W
!    F_supp(:, :, :) : 3D support array of size (dimx, dimy, 2*suppz)
!                      where the first suppz slices come from myid−1, and
!                      the next suppz slices come from myid+1.
!    id              : 1->U, 2->V, 3->W
!
!  Outputs:
!    F_subset(i, j−nb_start+1) : extracted scalar values for each weight i at body point j
!
!  Notes:
!    Uses u_proc, u_x_indices, etc., from module “globals”.
!    Assumes any nonlocal index zi is in 1..suppz, so layer = zi (prev rank)
!    or layer = suppz+zi (next rank).
!-------------------------------------------------------------------------------
  subroutine update_IB_subset(F, F_supp, F_subset, id)
    implicit none
  
    Real   (Int64), intent(in)  :: F(:, :, :)
    Real   (Int64), intent(in)  :: F_supp(:, :, :)   ! dims depend on id:
                                            !  id=1: dims = (nx, nyg, 2*suppz)
                                            !  id=2: dims = (nxg, ny, 2*suppz)
                                            !  id=3: dims = (nxg, nyg, 2*suppz)
    Real   (Int64), intent(out) :: F_subset(:, :)   ! dims = (nweights, nb_end−nb_start+1)
    integer,  intent(in) :: id
  
    integer :: i, j
    integer :: proc_idx,prev,next
    integer :: xi, yi, zi,zi_supp, zi_global
    integer :: layer
    Real (Int64) :: f_global

    ! prev = myid - 1
    ! if (prev < 0) prev = nprocs - 1
    ! next = myid + 1
    ! if (next == nprocs) next = 0
  
    do j = nb_start, nb_end
      do i = 1, nweights
  
        select case (id)
        case (1)
          proc_idx = u_proc(i, j)
          xi = u_x_indices(i, j)
          yi = u_y_indices(i, j)
          zi = u_z_local_indices(i, j)
          zi_supp = u_z_supp_idx(i, j)
          !zi_global = u_z_indices(i, j)
          !f_global=U_global(xi,yi,zi_global)
        case (2)
          proc_idx = v_proc(i, j)
          xi = v_x_indices(i, j)
          yi = v_y_indices(i, j)
          zi = v_z_local_indices(i, j)
          !zi_global = v_z_indices(i, j)
          zi_supp = v_z_supp_idx(i, j)
          !f_global=V_global(xi,yi,zi_global)
        case (3)
          proc_idx = w_proc(i, j)
          xi = w_x_indices(i, j)
          yi = w_y_indices(i, j)
          zi = w_z_local_indices(i, j)
          !zi_global = w_z_indices(i, j)
          zi_supp = w_z_supp_idx(i, j)
          !f_global=W_global(xi,yi,zi_global)
        case default
          print *, 'Error in update_IB_subset: invalid id =', id
          stop
        end select
  
        if (proc_idx == myid) then
          ! Local: read directly from F
          F_subset(i, j) = F(xi, yi, zi)
        else
          IF (zi_supp.le.0 .or. zi_supp .gt. 2*suppz+1) Then
            WRITE(*,*) 'myid',myid,'proc_idx',proc_idx,'zi_supp',zi_supp
            Stop 'Error: support index is wrong for left support cell'
          END IF
          F_subset(i, j) = F_supp(xi, yi, zi_supp)
          !print *, 'Error: unexpected proc_idx =', proc_idx, ' for myid =', myid
          !stop
        END IF

        ! debug line
        ! if ( ABS(F_subset(i,j)-f_global) .gt. 1E-8 ) then
        !   WRITE(*,*) 'F_subset does not equal to f_global,F_subset=',F_subset(i,j),'f_global=',f_global
        !   WRITE(*,*) 'ID',id,'myid',myid,'proc_id',proc_idx,'k_global',zi_global,'k_local',zi,'supp_idx',zi_supp
        !   !stop  
        ! end if
  
      end do
    end do
  
  end subroutine update_IB_subset

  
  


End Module immersed_boundary_operators

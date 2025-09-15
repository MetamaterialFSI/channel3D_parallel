!------------------------------------!
!   Module for boundary conditions   !
!------------------------------------!
Module boundary_conditions

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi

  ! prevent implicit typing
  Implicit None

Contains

  !--------------------------------------------------!
  !     Apply boundary conditions to velocities      !
  !              in the 3 directions                 !
  !--------------------------------------------------!
  Subroutine apply_boundary_conditions(U_, V_, W_)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: U_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: V_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: W_(:, :, :)

    ! interior region
    Call update_ghost_interior_planes(U_,1)
    Call update_ghost_interior_planes(V_,2)
    Call update_ghost_interior_planes(W_,3)
   
    ! apply periodicity in x
    Call apply_periodic_bc_x(U_,1)
    Call apply_periodic_bc_x(V_,2)
    Call apply_periodic_bc_x(W_,2)

    ! apply periodicity in z 
    Call apply_periodic_bc_z(U_,1)
    Call apply_periodic_bc_z(V_,2)
    Call apply_periodic_bc_z(W_,3)

    ! U boundary condition at the wall
    Call apply_Dirichlet_bc_y(U_,2)

    ! V boundary condition at the wall
    !Call apply_OppControl_bc_y(V,1)
    Call apply_Dirichlet_bc_y(V_,1)

    ! W boundary condition at the wall
    Call apply_Dirichlet_bc_y(W_,2)

    ! compute boundary conditions for pseudo-pressure
    !If ( i_BC_V==2 ) Then
    !   Call compute_pseudo_pressure_bc_for_robin_bc
    !End If

  End Subroutine apply_boundary_conditions
  
  !-------------------------------------------------!
  !                Periodicity in x                 !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at x faces          !
  !            id=2-> F defined at x centers        !
  ! Output: F                                       !
  !                                                 !
  !-------------------------------------------------!
  Subroutine apply_periodic_bc_x(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at x faces
       F( 1,:,:) = F(nx-1,:,:)
       F(nx,:,:) = F(   2,:,:)
    Else
       ! F defined at x centers
       F(    1,:,:) = F(nxg-2,:,:)
       F(nxg-1,:,:) = F(    2,:,:) ! see note*
       F(nxg  ,:,:) = F(    3,:,:)
    End If

    ! *Note: this is done in case the initial
    ! condition is not periodic. After the first
    ! step is no longer required
  End Subroutine apply_periodic_bc_x

  !-------------------------------------------------!
  ! Periodicity in z, MPI communication required    !
  !-------------------------------------------------!
  Subroutine apply_periodic_bc_z(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id
    Integer(Int64) :: n(3)
    ! save planes
    If ( myid==0 ) Then
      ! begin planes
      If (id == 3) Then
        ! F defined at z faces
        buffer_wi(:,:)   = F(:,:,2)
      Elseif (id == 1) Then 
        ! F defined at z centers
         buffer_ui(:,:,2) = F(:,:,2) 
         buffer_ui(:,:,3) = F(:,:,3)
      Elseif (id == 2) Then
        ! F defined at z centers
         buffer_vi(:,:,2) = F(:,:,2) 
         buffer_vi(:,:,3) = F(:,:,3)
      Elseif (id == 4) Then
        ! F defined at z centers
         buffer_ci(:,:,2) = F(:,:,2) 
         buffer_ci(:,:,3) = F(:,:,3)
      End If      

    End If

    If ( myid==nprocs-1 ) Then
       ! end planes
      If (id == 3) Then
        ! F defined at z faces
        buffer_we(:,:) = F(:,:, nz-1)
      Elseif (id == 1) Then
        ! F defined at z centers
        buffer_ue(:,:) = F(:,:,nzg-2)
      Elseif (id == 2) Then
        ! F defined at z centers
        buffer_ve(:,:) = F(:,:,nzg-2)
      Elseif (id == 4) Then
        ! F defined at z centers
        buffer_ce(:,:) = F(:,:,nzg-2)
      End If
    End If
    
    ! communicate planes
    If ( myid==0 ) Then
      If (id == 3) Then
        ! Send/receive W
        Call Mpi_sendrecv(buffer_wi, nxg*nyg, Mpi_real8, nprocs-1, 3,  &
        buffer_we, nxg*nyg, Mpi_real8, nprocs-1, 3, MPI_COMM_WORLD,    &
        istat, ierr)
      Elseif (id == 1) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_ui, nx*nyg*2, Mpi_real8, nprocs-1, 1, &
        buffer_ue, nx*nyg, Mpi_real8, nprocs-1, 1, MPI_COMM_WORLD,     &
        istat, ierr)
      Elseif (id == 2) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_vi, nxg*ny*2, Mpi_real8, nprocs-1, 2, &
        buffer_ve, nxg*ny, Mpi_real8, nprocs-1, 2, MPI_COMM_WORLD,     &
        istat, ierr)
      Elseif (id == 4) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_ci, nxg*nyg*2, Mpi_real8, nprocs-1, nprocs-1, &
        buffer_ce, nxg*nyg, Mpi_real8, nprocs-1, nprocs-1, MPI_COMM_WORLD,     &
        istat, ierr)
      End If
    End If

    If ( myid==nprocs-1 ) Then
      If (id == 3) Then
        ! Send/receive W
        Call Mpi_sendrecv(buffer_we, nxg*nyg, Mpi_real8, 0, 3, &
        buffer_wi, nxg*nyg, Mpi_real8, 0, 3, MPI_COMM_WORLD,   &
        istat, ierr)
      Elseif (id == 1) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_ue, nx*nyg, Mpi_real8, 0, 1,  &
        buffer_ui, nx*nyg*2, Mpi_real8, 0, 1, MPI_COMM_WORLD,  &
        istat, ierr)
      Elseif (id == 2) Then
        ! Send/receive V 
        Call Mpi_sendrecv(buffer_ve, nxg*ny, Mpi_real8, 0, 2,  &
        buffer_vi, nxg*ny*2, Mpi_real8, 0, 2, MPI_COMM_WORLD,  &
        istat, ierr)
      Elseif (id == 4) Then
        ! Send/receive U 
        Call Mpi_sendrecv(buffer_ce, nxg*nyg, Mpi_real8, 0, nprocs-1,  &
        buffer_ci, nxg*nyg*2, Mpi_real8, 0, nprocs-1, MPI_COMM_WORLD,  &
        istat, ierr)
      End If
    End If

    If (id == 4) Then
      If (myid == 0) Then
        Call Mpi_sendrecv(buffer_ci, nxg*nyg*2, Mpi_real8, nprocs-1, 4, &
        buffer_ce, nxg*nyg, Mpi_real8, nprocs-1, 4, MPI_COMM_WORLD,     &
        istat, ierr)
      End If
      If (myid == nprocs-1) Then
        Call Mpi_sendrecv(buffer_ce, nxg*nyg, Mpi_real8, 0, 4,  &
        buffer_ci, nxg*nyg*2, Mpi_real8, 0, 4, MPI_COMM_WORLD,  &
        istat, ierr)
      End If
   End If

    ! apply conditions
    If ( myid==0 ) Then
      If (id == 3) Then 
        F(:,:,1) = buffer_we(:,:)       ! W_global(:,:,nz_global-1)       
      Elseif (id == 1) Then
        F(:,:,1) = buffer_ue(:,:)       ! U_global(:,:,nzg_global-2)
      Elseif (id == 2) Then
        F(:,:,1) = buffer_ve(:,:)       ! U_global(:,:,nzg_global-2)
      Elseif (id == 4) Then
        F(:,:,1) = buffer_ce(:,:)       ! U_global(:,:,nzg_global-2)
      End If
    End If
    If ( myid==nprocs-1 ) Then
      If (id == 3) Then
        F(:,:,nz   ) = buffer_wi(:,:)   ! W_global(:,:,2) 
      Elseif (id == 1) Then
        F(:,:,nzg-1) = buffer_ui(:,:,2) ! U_global(:,:,2)
        F(:,:,nzg  ) = buffer_ui(:,:,3) ! U_global(:,:,3)
      Elseif (id == 2) Then
        F(:,:,nzg-1) = buffer_vi(:,:,2) ! U_global(:,:,2)
        F(:,:,nzg  ) = buffer_vi(:,:,3) ! U_global(:,:,3)
      Elseif (id == 4) Then
        F(:,:,nzg-1) = buffer_ci(:,:,2) ! U_global(:,:,2)
        F(:,:,nzg  ) = buffer_ci(:,:,3) ! U_global(:,:,3)
      End If    
    End If   
    Call Mpi_barrier(MPI_COMM_WORLD, ierr)
  End Subroutine apply_periodic_bc_z

  !-------------------------------------------------!
  !        Dirichlet boundary condition in y        !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at y faces          !
  !            id=2-> F defined at y centers        !
  ! Output: F                                       !
  !                                                 !
  ! For now only wall                               !
  !-------------------------------------------------!
  Subroutine apply_Dirichlet_bc_y(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces
       F(:, 1,:) = 0d0
       F(:,ny,:) = 0d0
    Else
       ! F defined at y centers
       F(:,  1,:) = -F(:,    2,:)
       F(:,nyg,:) = -F(:,nyg-1,:)
    End If

  End Subroutine apply_Dirichlet_bc_y

  !-------------------------------------------------!
  !   Opposition control boundary condition in y    !
  !          No MPI communication required          !
  !                                                 !
  ! Input:  F  (array to apply boundary conditions) !
  !         id id=1-> F defined at y faces          !
  !            id=2-> SHOULD NOT BE USED FOR THIS   !
  ! Output: F                                       !
  !                                                 !
  ! For now only wall-normal velocity               !
  !-------------------------------------------------!
  Subroutine apply_OppControl_bc_y(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)
    Integer(Int32), Intent(In)    :: id

    If ( id==1 ) Then
       ! F defined at y faces
       F(:, 1,:) = Vw(:,1,:) 
       F(:,ny,:) = Vw(:,2,:)
    Else 
       Stop 'Error: Opposition contol on wrong comp'
    End If

  End Subroutine apply_OppControl_bc_y


  !--------------------------------------------------!
  !          Update ghost interior planes            !
  !--------------------------------------------------!
  Subroutine update_ghost_interior_planes(F,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)    
    Integer(Int32), Intent(In)    :: id

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom
    
    If (id == 1) Then
      !----------------------update U-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_us = F(:,:,nzg-1) ! send buffer
      Call Mpi_sendrecv(buffer_us, nx*nyg, Mpi_real8, sendto, tagto,        &
           buffer_ur, nx*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_ur ! received buffer
      
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
      buffer_us = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_us, nx*nyg, Mpi_real8, sendto, tagto,        &
           buffer_ur, nx*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_ur ! received buffer

    Elseif (id == 2) Then
      !----------------------update V-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_vs = F(:,:,nzg-1) ! send buffer
      Call Mpi_sendrecv(buffer_vs, nxg*ny, Mpi_real8, sendto, tagto,        &
           buffer_vr, nxg*ny, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_vr ! received buffer
      
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
      buffer_vs = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_vs, nxg*ny, Mpi_real8, sendto, tagto,        &
           buffer_vr, nxg*ny, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_vr ! received buffer
      
    Elseif (id == 3) Then
      !----------------------update W-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_ws = F(:,:,nz-1)    ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_wr ! received buffer
      
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
      buffer_ws = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nz) = buffer_wr ! received buffer     

    Elseif (id == 4) Then
      !----------------------update term-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         recvfrom = MPI_PROC_NULL
         tagfrom  = MPI_ANY_TAG
      End If
      If ( myid==nprocs-1 ) Then
         sendto = MPI_PROC_NULL
         tagto  = 0
      End If
      buffer_ws = F(:,:,nzg-1)    ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=0 ) F(:,:,1) = buffer_wr ! received buffer
      
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
      buffer_ws = F(:,:,2)  ! send buffer
      Call Mpi_sendrecv(buffer_ws, nxg*nyg, Mpi_real8, sendto, tagto,        &
           buffer_wr, nxg*nyg, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      If ( myid/=nprocs-1 ) F(:,:,nzg) = buffer_wr ! received buffer     
    End if  
    
  End Subroutine update_ghost_interior_planes

  !-------------------------------------------------------------------------------
  !  interior_planes_update_support
  !
  !  Purpose:
  !    Sends interior velocity planes to update the ghost/support layers of
  !    neighboring processors (F_supp). Handles periodic boundaries.
  !
  !  Inputs:
  !    F       : Field array (U, V, or W)
  !    F_supp  : Support array to be updated
  !    id      : Component ID (1=U, 2=V, 3=W)
  !
  !  Notes:
  !    - Periodic z-direction assumed
  !    - Communication occurs in both directions
  !-------------------------------------------------------------------------------
  Subroutine interior_planes_update_support(F,F_supp,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)    
    Real   (Int64), Intent(InOut) :: F_supp(:,:,:)    
    Integer(Int32), Intent(In)    :: id

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom
    
    If (id == 1) Then
      !----------------------update U-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
        recvfrom = nprocs-1
        tagfrom  = 0
        buffer_usupp_s = F(:,:,nzg-suppz-1:nzg-1) ! send buffer
      ELSEIf ( myid==nprocs-1 ) Then
        sendto = 0
        tagto  = 0
        if (nzg-suppz-1 .lt. 2) Stop 'Error: nzg-suppz-1 less than 2'
        buffer_usupp_s = F(:,:,nzg-suppz-2:nzg-2) ! send buffer
      ELSE
        buffer_usupp_s = F(:,:,nzg-suppz-1:nzg-1) ! send buffer
      End If
      
      Call Mpi_sendrecv(buffer_usupp_s, nx*nyg*(suppz+1), Mpi_real8, sendto, tagto,        &
      buffer_usupp_r, nx*nyg*(suppz+1), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      F_supp(:,:,1:suppz+1) = buffer_usupp_r ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
        sendto = nprocs-1
        tagto  = nprocs-1
        if (2+suppz .gt. nzg-1) Stop 'Error: 2+suppz greater than nzg-1'
        buffer_usupp_s(:,:,1:suppz) = F(:,:,3:2+suppz)  ! send buffer
      ELSEIf ( myid==nprocs-1 ) Then
        recvfrom = 0
        tagfrom  = nprocs-1
        buffer_usupp_s(:,:,1:suppz) = F(:,:,2:1+suppz)  ! send buffer
      ELSE
        buffer_usupp_s(:,:,1:suppz) = F(:,:,2:1+suppz)  ! send buffer
      End If
      Call Mpi_sendrecv(buffer_usupp_s(:,:,1:suppz), nx*nyg*suppz, Mpi_real8, sendto, tagto,        &
      buffer_usupp_r(:,:,1:suppz), nx*nyg*suppz, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      F_supp(:,:,suppz+2:2*suppz+1) = buffer_usupp_r(:,:,1:suppz) ! received buffer

    Elseif (id == 2) Then
      !----------------------update V-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
        recvfrom = nprocs-1
        tagfrom  = 0
        buffer_vsupp_s = F(:,:,nzg-suppz-1:nzg-1) ! send buffer
      ELSEIf ( myid==nprocs-1 ) Then
        sendto = 0
        tagto  = 0
        if (nzg-suppz-2 .lt. 2) Stop 'Error: nzg-suppz-2 less than 2'
        buffer_vsupp_s = F(:,:,nzg-suppz-2:nzg-2) ! send buffer
      ELSE
        buffer_vsupp_s = F(:,:,nzg-suppz-1:nzg-1) ! send buffer
      End If
      Call Mpi_sendrecv(buffer_vsupp_s, nxg*ny*(suppz+1), Mpi_real8, sendto, tagto,        &
      buffer_vsupp_r, nxg*ny*(suppz+1), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      F_supp(:,:,1:suppz+1) = buffer_vsupp_r ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
        sendto = nprocs-1
        tagto  = nprocs-1
        buffer_vsupp_s(:,:,1:suppz) = F(:,:,3:2+suppz)  ! send buffer
      ELSEIf ( myid==nprocs-1 ) Then
        recvfrom = 0
        tagfrom  = nprocs-1
        buffer_vsupp_s(:,:,1:suppz) = F(:,:,2:1+suppz)  ! send buffer
      ELSE
        buffer_vsupp_s(:,:,1:suppz) = F(:,:,2:1+suppz)  ! send buffer
      End If
        Call Mpi_sendrecv(buffer_vsupp_s(:,:,1:suppz), nxg*ny*suppz, Mpi_real8, sendto, tagto,        &
        buffer_vsupp_r(:,:,1:suppz), nxg*ny*suppz, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
             istat, ierr)   
        F_supp(:,:,suppz+2:2*suppz+1) = buffer_vsupp_r(:,:,1:suppz) ! received buffer
      
    Elseif (id == 3) Then
      !----------------------update W-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
        recvfrom = nprocs-1
        tagfrom  = 0
      ELSEIf ( myid==nprocs-1 ) Then
        sendto = 0
        tagto  = 0
      ELSE
      End If
      if (nz-suppz-1 .lt. 2) Stop 'Error: nz-suppz-1 less than 2'
      buffer_wsupp_s = F(:,:,nz-suppz-1:nz-1)    ! send buffer
      Call Mpi_sendrecv(buffer_wsupp_s, nxg*nyg*(suppz+1), Mpi_real8, sendto, tagto,        &
      buffer_wsupp_r, nxg*nyg*(suppz+1), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      F_supp(:,:,1:suppz+1) = buffer_wsupp_r ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
        sendto = nprocs-1
        tagto  = nprocs-1
      ELSEIf ( myid==nprocs-1 ) Then
        recvfrom = 0
        tagfrom  = nprocs-1
      End If
      buffer_wsupp_s(:,:,1:suppz) = F(:,:,2:suppz+1)  ! send buffer
      Call Mpi_sendrecv(buffer_wsupp_s(:,:,1:suppz), nxg*nyg*suppz, Mpi_real8, sendto, tagto,        &
      buffer_wsupp_r(:,:,1:suppz), nxg*nyg*suppz, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      F_supp(:,:,suppz+2:2*suppz+1) = buffer_wsupp_r(:,:,1:suppz) ! received buffer     
    End if  
    Call Mpi_barrier(MPI_COMM_WORLD, ierr) 
    
  End Subroutine interior_planes_update_support

  !-------------------------------------------------------------------------------
  !  support_update_interior_planes
  !
  !  Purpose:
  !    Updates interior planes (F) using ghost/support planes (F_supp)
  !    from neighboring processors in the z-direction. This subroutine
  !    performs a **summation** of received support layers onto the
  !    corresponding interior layers — additive contribution is critical
  !    for correct formulation (e.g., immersed boundary forcing).
  !
  !  Inputs:
  !    F       : Field array (U, V, or W) to be updated
  !    F_supp  : Support/ghost array
  !    id      : Component ID (1=U, 2=V, 3=W)
  !
  !  Notes:
  !    - Assumes periodic domain in z-direction.
  !    - Performs additive update, **not overwrite**.
  !-------------------------------------------------------------------------------
 Subroutine support_update_interior_planes(F,F_supp,id)

    Real   (Int64), Intent(InOut) :: F(:,:,:)    
    Real   (Int64), Intent(InOut) :: F_supp(:,:,:)    
    Integer(Int32), Intent(In)    :: id

    Integer(Int32) :: sendto, recvfrom
    Integer(Int32) :: tagto,  tagfrom
    
    !WRITE(*,*) 'supp2interior myid:',myid,'id:',id
    If (id == 1) Then
      !----------------------update U-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
         recvfrom = nprocs-1
         tagfrom  = 0
         buffer_usupp_s(:,:,1:suppz) = F_supp(:,:,suppz+2:2*suppz+1) ! send buffer
         Call Mpi_sendrecv(buffer_usupp_s(:,:,1:suppz), nx*nyg*suppz, Mpi_real8, sendto, tagto,        &
            buffer_usupp_r, nx*nyg*(suppz+1), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr) 
          F(:,:,2:2+suppz) = F(:,:,2:2+suppz)+ buffer_usupp_r ! received buffer
      ELSEIf ( myid==nprocs-1 ) Then
        !Also transfering the periodic slices to the first partition
         sendto = 0
         tagto  = 0
         buffer_usupp_s(:,:,2:suppz+1) = F_supp(:,:,suppz+2:2*suppz+1) ! send buffer
         buffer_usupp_s(:,:,1) = F(:,:,nxg-1);
         Call Mpi_sendrecv(buffer_usupp_s, nx*nyg*(suppz+1), Mpi_real8, sendto, tagto,        &
            buffer_usupp_r(:,:,1:suppz), nx*nyg*suppz, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr) 
         F(:,:,2:1+suppz) = F(:,:,2:1+suppz)+ buffer_usupp_r(:,:,1:suppz) ! received buffer
      ELSE
        buffer_usupp_s(:,:,1:suppz) = F_supp(:,:,suppz+2:2*suppz+1) ! send buffer
        Call Mpi_sendrecv(buffer_usupp_s(:,:,1:suppz), nx*nyg*suppz, Mpi_real8, sendto, tagto,        &
        buffer_usupp_r(:,:,1:suppz), nx*nyg*suppz, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr) 
        F(:,:,2:1+suppz) = F(:,:,2:1+suppz)+ buffer_usupp_r(:,:,1:suppz) ! received buffer
      End If
      
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
         sendto = nprocs-1
         tagto  = nprocs-1
      End If
      If ( myid==nprocs-1 ) Then
         recvfrom = 0
         tagfrom  = nprocs-1
      End If
      buffer_usupp_s = F_supp(:,:,1:suppz+1)  ! send buffer
      Call Mpi_sendrecv(buffer_usupp_s, nx*nyg*(suppz+1), Mpi_real8, sendto, tagto,        &
      buffer_usupp_r, nx*nyg*(suppz+1), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      IF (myid==nprocs-1) THEN
        if (nzg-suppz-2 .lt. 2) Stop 'Error: nzg-suppz-1 less than 2'
        F(:,:,nzg-suppz-2:nzg-2) = F(:,:,nzg-suppz-2:nzg-2)+buffer_usupp_r ! received buffer
      ELSE
        F(:,:,nzg-suppz-1:nzg-1) = F(:,:,nzg-suppz-1:nzg-1)+buffer_usupp_r ! received buffer
      END IF
      

    Elseif (id == 2) Then
      !----------------------update V-----------------------!
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
        recvfrom = nprocs-1
        tagfrom  = 0
        buffer_vsupp_s(:,:,1:suppz) = F_supp(:,:,suppz+2:2*suppz+1) ! send buffer
        Call Mpi_sendrecv(buffer_vsupp_s(:,:,1:suppz), nxg*ny*suppz, Mpi_real8, sendto, tagto,        &
        buffer_vsupp_r, nxg*ny*(suppz+1), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
          istat, ierr) 
         F(:,:,2:2+suppz) = F(:,:,2:2+suppz)+ buffer_vsupp_r ! received buffer
     ELSEIf ( myid==nprocs-1 ) Then
       !Also transfering the periodic slices to the first partition
        sendto = 0
        tagto  = 0
        buffer_vsupp_s(:,:,2:suppz+1) = F_supp(:,:,suppz+2:2*suppz+1) ! send buffer
        buffer_vsupp_s(:,:,1) = F(:,:,nxg-1);
        Call Mpi_sendrecv(buffer_vsupp_s, nxg*ny*(suppz+1), Mpi_real8, sendto, tagto,        &
        buffer_vsupp_r(:,:,1:suppz), nxg*ny*suppz, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
          istat, ierr) 
        F(:,:,2:1+suppz) = F(:,:,2:1+suppz)+ buffer_vsupp_r(:,:,1:suppz) ! received buffer
     ELSE
      buffer_vsupp_s(:,:,1:suppz) = F_supp(:,:,suppz+2:2*suppz+1) ! send buffer
       Call Mpi_sendrecv(buffer_vsupp_s(:,:,1:suppz), nxg*ny*suppz, Mpi_real8, sendto, tagto,        &
       buffer_vsupp_r(:,:,1:suppz), nxg*ny*suppz, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
          istat, ierr) 
       F(:,:,2:1+suppz) = F(:,:,2:1+suppz)+ buffer_vsupp_r(:,:,1:suppz) ! received buffer
     End If
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
        sendto = nprocs-1
        tagto  = nprocs-1
     End If
     If ( myid==nprocs-1 ) Then
        recvfrom = 0
        tagfrom  = nprocs-1
     End If
      buffer_vsupp_s = F_supp(:,:,1:1+suppz)  ! send buffer
      Call Mpi_sendrecv(buffer_vsupp_s, nxg*ny*(suppz+1), Mpi_real8, sendto, tagto,        &
      buffer_vsupp_r, nxg*ny*(suppz+1), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      IF (myid==nprocs-1) THEN
        if (nzg-suppz-2 .lt. 2) Stop 'Error: nzg-suppz-2 grater than 2'
        F(:,:,nzg-suppz-2:nzg-2) = F(:,:,nzg-suppz-2:nzg-2)+buffer_vsupp_r ! received buffer
      ELSE
        F(:,:,nzg-suppz-1:nzg-1) = F(:,:,nzg-suppz-1:nzg-1)+buffer_vsupp_r ! received buffer
      END IF
      
    Elseif (id == 3) Then
      !----------------------update W-----------------------!
      ! send to top processor, receive from bottom one
      sendto   = myid + 1
      tagto    = myid + 1
      recvfrom = myid - 1
      tagfrom  = myid 
      If ( myid==0 ) Then 
        recvfrom = nprocs-1
        tagfrom  = 0
     End If
     If ( myid==nprocs-1 ) Then
        sendto = 0
        tagto  = 0
     End If
      buffer_wsupp_s(:,:,1:suppz) = F_supp(:,:,suppz+2:2*suppz+1)  ! send buffer
      Call Mpi_sendrecv(buffer_wsupp_s(:,:,1:suppz), nxg*nyg*suppz, Mpi_real8, sendto, tagto,        &
      buffer_wsupp_r(:,:,1:suppz), nxg*nyg*suppz, Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      F(:,:,2:1+suppz) = F(:,:,2:1+suppz)+ buffer_wsupp_r(:,:,1:suppz) ! received buffer
      
      ! send to bottom processor, receive from top one
      sendto   = myid - 1
      tagto    = myid - 1
      recvfrom = myid + 1
      tagfrom  = myid 
      If ( myid==0 ) Then
        sendto = nprocs-1
        tagto  = nprocs-1
     End If
     If ( myid==nprocs-1 ) Then
        recvfrom = 0
        tagfrom  = nprocs-1
     End If
      buffer_wsupp_s = F_supp(:,:,1:suppz+1)  ! send buffer
      Call Mpi_sendrecv(buffer_wsupp_s, nxg*nyg*(suppz+1), Mpi_real8, sendto, tagto,        &
      buffer_wsupp_r, nxg*nyg*(suppz+1), Mpi_real8, recvfrom, tagfrom, MPI_COMM_WORLD, &
           istat, ierr)   
      F(:,:,nz-suppz-1:nz-1) = F(:,:,nz-suppz-1:nz-1)+ buffer_wsupp_r ! received buffer     
    End if 
    Call Mpi_barrier(MPI_COMM_WORLD, ierr) 
    !WRITE(*,*) 'Done:supp2interior myid:',myid,'id:',id
    
  End Subroutine support_update_interior_planes

End Module boundary_conditions

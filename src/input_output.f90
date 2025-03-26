!--------------------------------------!
!          Module for I/O              !
!--------------------------------------!
Module input_output

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use ifport
  Use pressure

  ! prevent implicit typing
  Implicit None

Contains

  !----------------------------------------!
  !  Read input parameters from a txt file !
  !----------------------------------------!
  Subroutine read_input_parameters

    Character(200) :: dummy_line, msg
    Real(Int64)    :: Rossby_plus, utau_
    Integer(Int32) :: ioerr, iounit

    namelist /params/ &
      nx_global, ny_global, nz_global, &
      nxb, nzb, nd, &
      CFL, &
      nu, &
      dPdx, dPdz, x_mass_cte, y_mass_cte, &
      nsteps, nsave, nstats, nmonitor, &
      filein, fileout, &
      nstep_init, &
      init_type, grid_type, body_type

    ! processor 0 reads the data
    If ( myid==0 ) Then
      Write(*,*) 'reading input parameters...'

      Open(newunit=iounit, file=fileparams, status="old", action="read", iostat=ioerr, iomsg=msg)
      If (ierr /= 0) then
          Print *, "Error opening input file. IOSTAT =", ioerr, " IOMSG = ", msg
          call abort
      End If
      Read(iounit, Nml=params, IOSTAT=ioerr, IOMSG=msg)
      If (ierr /= 0) then
          Print *, "Error reading namelist. IOSTAT =", ioerr, " IOMSG = ", msg
          call abort
      End If

      utau_    = dPdx**0.5d0
      dPdx_ref = dPdx

      Print params
       
    End If

    ! broadcast data to all processors
    Call Mpi_bcast ( nx_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ny_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nz_global,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( nxb,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( nzb,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  nd,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (      CFL,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (       nu,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (     dPdx,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (     dPdz,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( dPdx_ref,1,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast (  nstep_init,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (      nsteps,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (       nsave,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (      nstats,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (    nmonitor,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   init_type,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   grid_type,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (   body_type,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  x_mass_cte,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  y_mass_cte,1,MPI_integer,0,MPI_COMM_WORLD,ierr )

  End Subroutine read_input_parameters

  !--------------------------------------!
  !    Generates a grid for channel      !
  !                                      !
  ! Output: x_global, y_global, z_global !
  !                                      !
  !--------------------------------------!
  Subroutine create_grid

    Integer(Int32) :: i
    Real   (Int64) :: alpha

    Do i=1,nx_global
       x_global(i) = Real(i-1,8)*0.2d0
    End Do
    x_global = 2d0*pi*x_global/x_global(nx_global-1)

    Do i=1,nz_global
       z_global(i) = Real(i-1,8)*0.1d0
    End Do
    z_global = 1d0*z_global/z_global(nz_global-1)

    Select Case (grid_type)
      Case (0) ! Uniform grid
        If ( myid==0 ) Write(*,*) 'Generating uniform y grid'
        Do i=1,ny_global
           y_global(i) = Real(i-1,8)*0.3d0
        End Do
        y_global = 2d0*y_global/Maxval(y_global)

      Case (1) ! Stretched grid
        If ( myid==0 ) Write(*,*) 'Generating stretched y grid'
        Do i=1,ny_global
           y_global(i) = Real(i-1,8)*0.3d0
        End Do
        y_global = 2d0*y_global/Maxval(y_global)-1d0

        alpha = 2.6d0
        Do i=1,ny_global
           y_global(i) = dtanh(alpha*y_global(i))/dtanh(alpha)
        End Do

        y_global = y_global - Minval(y_global)
        y_global = y_global*2d0/Maxval(y_global)
    End Select

  End Subroutine create_grid

  !------------------------------------------------!
  !    Generates an initial condition for channel  !
  !                                                !
  ! Output: U,V,W                                  !
  !                                                !
  !------------------------------------------------!
  Subroutine init_flow
  
    Integer(Int32) :: ii, jj, kk
    Real(Int64) :: ym_val

    Select Case (init_type)
      Case (0) ! read input data from file
        If ( myid==0 ) Write(*,*) 'Reading input data'
        Call read_input_data
      
      Case (1) ! create grid and initialize velocity to random values
        If ( myid==0 ) Write(*,*) 'Generating zero initial condition'
        Call create_grid
        U = 0d0; V = 0d0; W = 0d0

      Case (2) ! create grid and initialize velocity to zero
        If ( myid==0 ) Write(*,*) 'Generating random initial condition'
        Call create_grid

        ! U
        Do jj=1,ny_global-1
          ym_val = 0.5d0 * (y_global(jj) + y_global(jj + 1))
          U(:,jj+1,:) = dpdx / (2d0 * nu) * ym_val * (2d0 - ym_val)
        end Do
        Do ii=1,nx_global
          Do jj=1,nyg_global
             Do kk=1,nzg
               U(ii,jj,kk) = U(ii,jj,kk) + 0.5*(rand()-0.5)
             End Do
          End Do
        End Do
        U(:,1,:) = -U(:,2,:)
        U(:,ny_global+1,:) = -U(:,ny_global,:)

        ! V
        V = 0d0
        Do ii=1,nxg_global
          Do jj=1,ny_global
             Do kk=1,nzg
               V(ii,jj,kk) = V(ii,jj,kk) + 0.5*(rand()-0.5)
             End Do
          End Do
        End Do
        V(:,1,:) = 0d0
        V(:,ny_global,:) = 0d0

        ! W
        W = 0d0
        Do ii=1,nxg_global
          Do jj=1,ny_global
             Do kk=1,nz
               W(ii,jj,kk) = W(ii,jj,kk) + 0.5*(rand()-0.5)
             End Do
          End Do
        End Do
        W(:,1,:) = -W(:,2,:)
        W(:,ny_global+1,:) = -W(:,ny_global,:)
    End Select

    If ( myid==0 ) Then
      Write(*,*) 'Max U',MaxVal(U)
      Write(*,*) 'Max V',MaxVal(V)
      Write(*,*) 'Max W',MaxVal(W)

      Write(*,*) 'Mean U',sum(U)/Real(nx_global*nyg_global*nzg_global,8)
      Write(*,*) 'Mean V',sum(V)/Real(nxg_global*ny_global*nzg_global,8)
      Write(*,*) 'Mean W',sum(W)/Real(nxg_global*nyg_global*nz_global,8)
    End If

  End Subroutine init_flow

  !--------------------------------------------!
  !    Read binary snapshot: mesh, U,V and W   !
  !                                            !
  ! Input:  filein                             !
  ! Output: U,V,W,x,y,z                        !
  !                                            !
  !--------------------------------------------!
  Subroutine read_input_data

    Integer(Int32) ::  nx_global_f,  ny_global_f,  nz_global_f, iproc, nze, nzge
    Integer(Int32) :: nxm_global_f, nym_global_f, nzm_global_f, nn(3), ndum
    Integer(Int64) :: pos_header, nsize_U, nsize_V, ii, jj, kk

    ! processor 0 Reads the all the data
    If ( myid==0 ) Then

      Write(*,*) 'reading ',Trim(Adjustl(filein)),'...'
      Open(1,file=filein,access='stream',form='unformatted',action='Read',convert='big_endian')

      ! mesh
      Read(1) nx_global_f
      If ( nx_global_f/=nx_global ) Stop 'nx_f/=nx'
      Read(1) x_global

      Read(1) ny_global_f
      If ( ny_global_f/=ny_global ) Stop 'ny_f/=ny'
      Read(1) y_global

      Read(1) nz_global_f
      If ( nz_global_f/=nz_global ) Stop 'nz_f/=nz'
      Read(1) z_global

      Read(1) nxm_global_f
      If ( nxm_global_f/=nxm_global ) Stop 'nxm_f/=nxm'
      Read(1) xm_global

      Read(1) nym_global_f
      If ( nym_global_f/=nym_global ) Stop 'nym_f/=nym'
      Read(1) ym_global

      Read(1) nzm_global_f
      If ( nzm_global_f/=nzm_global ) Stop 'nzm_f/=nzm'
      Read(1) zm_global

      ! get header position and size
      Inquire(1,pos=pos_header)
      pos_header = pos_header - 1
      nsize_U    = nx_global*nyg_global*nzg_global*8
      nsize_V    = nxg_global*ny_global*nzg_global*8

    End If

    ! U
    If ( myid==0 ) Then
      ! read dummy
      Read(1) nn 
      If ( nn(1)/=nx_global .or. nn(2)/=nyg_global .or. nn(3)/=nzg_global ) Then 
        Write(*,*) 'nn',nn
        Stop 'Error! wrong size in input file (U)'
      End If
      ! read data for processor 0
      nzge = kg2_global(myid) - kg1_global(myid) + 1
      Read(1) U(:,:,1:nzge)
      ! data for processor n>0    
      Do iproc = 1, nprocs-1
        nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
        If ( iproc<nprocs-1 ) Then
          ndum = fseek(1,-2*nx_global*nyg_global*8,seek_cur) ! ghost cell
          Read(1) Uo(:,:,1:nzge)
          Call Mpi_send(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        Else ! especial case: U has different size for last processor
          ndum = fseek(1,-2*nx_global*nyg_global*8,seek_cur) ! ghost cell
          Read(1) Uoo(:,:,1:nzge)
          Call Mpi_send(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        End If
      Enddo       
    Else
      Call Mpi_recv(U,nx*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! V
    If ( myid==0 ) Then
      ! go to correct position. I dont know, if I dont do this it gets lost sometimes
      ndum = fseek(1,pos_header+3*4+nsize_U,seek_set)
      ! read dummy
      Read(1) nn
      If ( nn(1)/=nxg_global .or. nn(2)/=ny_global .or. nn(3)/=nzg_global ) Then 
         Write(*,*) 'nn',nn
         Stop 'Error! wrong size in input file (V)'
      End If
      ! read data for processor 0
      nzge = kg2_global(myid) - kg1_global(myid) + 1
      Read(1) V(:,:,1:nzge)
      ! data for processor n>0    
      Do iproc = 1, nprocs-1
        nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
        If ( iproc<nprocs-1 ) Then
          ndum = fseek(1,-2*nxg_global*ny_global*8,seek_cur) ! ghost cell
          Read(1) Vo(:,:,1:nzge) 
          Call Mpi_send(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        Else ! especial case: V has different size for last processor
          ndum = fseek(1,-2*nxg_global*ny_global*8,seek_cur) ! ghost cell
          Read(1) Voo(:,:,1:nzge) 
          Call Mpi_send(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        End If
      Enddo       
    Else
      Call Mpi_recv(V,nxg*ny*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! W
    If ( myid==0 ) Then
      ! go to correct position. I dont know, if I dont do this it gets lost sometimes
      ndum = fseek(1,pos_header+3*4+nsize_U+3*4+nsize_V,seek_set)
      ! read dummy
      Read(1) nn
      If ( nn(1)/=nxg_global .or. nn(2)/=nyg_global .or. nn(3)/=nz_global ) Then 
        Write(*,*) 'nn',nn
        Stop 'Error! wrong size in input file (W)'
      End If
      ! read data for processor 0
      nzge = k2_global(myid) - k1_global(myid) + 1
      Read(1) W(:,:,1:nzge)
      ! data for processor n>0    
      Do iproc = 1, nprocs-1
        nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
        If ( iproc<nprocs-1 ) Then
          ndum = fseek(1,-2*nxg_global*nyg_global*8,seek_cur) ! ghost cell
          Read(1) Wo(:,:,1:nzge)
          Call Mpi_send(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        Else ! especial case: W has different size for last processor
          ndum = fseek(1,-2*nxg_global*nyg_global*8,seek_cur) ! ghost cell
          Read(1) Woo(:,:,1:nzge)
          Call Mpi_send(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,ierr)
        End If
      Enddo       
    Else
      Call Mpi_recv(W,nxg*nyg*nz,Mpi_real8,0,myid,MPI_COMM_WORLD,istat,ierr)
    Endif

    ! close file
    If (myid==0) Then
      Close(1)
    End If

    ! send data to all other processors
    ! mesh
    Call Mpi_bcast ( x_global,nx_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( y_global,ny_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( z_global,nz_global,MPI_real8,0,MPI_COMM_WORLD,ierr )

    Call Mpi_bcast ( xm_global,nxm_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( ym_global,nym_global,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast ( zm_global,nzm_global,MPI_real8,0,MPI_COMM_WORLD,ierr ) 


    ! set solution for zero step
    Uo = U
    Vo = V
    Wo = W

  End Subroutine read_input_data

  !--------------------------------------------!
  !    write binary snapshot: mesh, U,V and W  !
  !                                            !
  ! Input: U,V,W,x,y,z,xm,ym,zm                !
  ! Output: fileout                            !
  !                                            !
  !--------------------------------------------!
  Subroutine output_data

    Character(200)   :: fname
    Character(8)     :: ext
    Integer  (Int32) :: iproc, nze, nzge
    
    If ( Mod(istep,nsave)==0 ) then

      ! P
      !If (pressure_computed==.False.) Then
      !   Call compute_pressure
      !End If

      ! processor 0 writes the data
      If ( myid==0 ) Then
        
        Write(ext,'(I8)') istep + nstep_init
        
        fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))
        Write(*,*) 'writing ',Trim(Adjustl(fname))
        Open(1,file=fname,access='stream',form='unformatted',action='write',convert='big_endian')
        
        ! mesh
        Write(1) Shape(x_global), x_global
        Write(1) Shape(y_global), y_global
        Write(1) Shape(z_global), z_global
        
        Write(1) Shape(xm_global), xm_global
        Write(1) Shape(ym_global), ym_global
        Write(1) Shape(zm_global), zm_global          
       
      End If

      ! U
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(U,nx*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write U size
        Write(1) nx_global,nyg_global,nzg_global
        ! processor 0 writes its data
        Write(1) U(:,:,1:nzg-1) 
        ! processor 0 receives and writes rest data
        Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Uo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Uo(:,:,2:nzge-1)
          Else
            Call Mpi_recv(Uoo,nx*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Uoo(:,:,2:nzge)
          End If
        End Do
      Endif

      ! V
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(V,nxg*ny*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write V size
        Write(1) nxg_global,ny_global,nzg_global
        ! processor 0 writes its data
        Write(1) V(:,:,1:nzg-1)
        ! processor 0 receives and write rest data
        Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Vo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Vo(:,:,2:nzge-1)
          Else
            Call Mpi_recv(Voo,nxg*ny*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Voo(:,:,2:nzge)
          End If
        End Do
      Endif

      ! W
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(W,nxg*nyg*nz,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write W size
        Write(1) nxg_global,nyg_global,nz_global
        ! processor 0 writes its data
        Write(1) W(:,:,1:nz-1)
        ! processor 0 receives and writes rest data
        Do iproc = 1, nprocs-1
          nze = k2_global(iproc) - k1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Wo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Wo(:,:,2:nze-1)
          Else
            Call Mpi_recv(Woo,nxg*nyg*nze,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Woo(:,:,2:nze)
          End If
        End Do
      Endif

      ! P
      If ( myid/=0 ) Then
        ! data from processor n>0    
        Call Mpi_send(P,nxg*nyg*nzg,Mpi_real8,0,myid,MPI_COMM_WORLD,ierr)
      Else
        ! write P size
        Write(1) nxg_global,nyg_global,nzg_global
        ! processor 0 writes its data
        Write(1) P(:,:,1:nzg-1)
        ! processor 0 receives and write rest data
        Do iproc = 1, nprocs-1
          nzge = kg2_global(iproc) - kg1_global(iproc) + 1 ! local size in z for processor iproc
          If ( iproc<nprocs-1 ) Then
            Call Mpi_recv(Po,nxg*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Po(:,:,2:nzge-1)
          Else
            Call Mpi_recv(Poo,nxg*nyg*nzge,Mpi_real8,iproc,iproc,MPI_COMM_WORLD,istat,ierr)
            Write(1) Poo(:,:,2:nzge)
          End If
        End Do
      Endif
         
      ! close file
      If (myid==0) Then
        Close(1)
      End If
      
    End If
       
  End Subroutine output_data

  !----------------------------------------------!
  !   Write some basic statistics in a txt file  !
  !----------------------------------------------!
  Subroutine output_statistics

    Character(200) :: fname
    Character(8)   :: ext
    Integer(Int32) :: jj

    If ( myid==0 ) Then

       Write(ext,'(I8)') istep + nstep_init
       
       fname = Trim(Adjustl(fileout))//'.'//Trim(Adjustl(ext))//'.stats.txt'
       Write(*,*) 'writing ',Trim(Adjustl(fname))
       Open(3,file=fname,form='formatted',action='write') 
       Write(3,'(A,4F15.8,4I)') '%',t, Retau, utau, nu, nx_global, ny_global, nz_global, istep
       Do jj=1,nyg
          Write(3,'(8F15.8)') yg(jj), Umean(jj), Vmean(jj), Wmean(jj), U2mean(jj), V2mean(jj), W2mean(jj), UVmean(jj)
       End Do
       Close(3)

    End If

  End Subroutine output_statistics

End Module input_output

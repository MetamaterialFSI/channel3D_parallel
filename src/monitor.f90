!--------------------------------------------!
! Module to monitor status of the simulation !
!--------------------------------------------!
Module monitor

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use projection, Only : check_divergence

  ! prevent implicit typing
  Implicit None

Contains

  !------------------------------------------!
  ! Output some key values during simulation !
  !------------------------------------------!
  Subroutine output_monitor

    Real(Int64) ::  maxU,  maxV,  maxW, local_sum
    Real(Int64) :: meanU, meanV, meanW
    Real(Int64) :: max_divergence

    If ( Mod(istep,nmonitor)==0 ) Then

       ! compute mean values
       local_sum = Sum ( U(2:nx-1,2:nyg-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,meanU,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

       local_sum = Sum ( V(2:nxg-1,2:ny-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,meanV,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

       local_sum = Sum ( W(2:nxg-1,2:nyg-1,2:nz-1) )
       Call MPI_Reduce (local_sum,meanW,1,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

       ! compute maximum values
       local_sum = Maxval ( U(2:nx-1,2:nyg-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,maxU,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

       local_sum = Maxval ( V(2:nxg-1,2:ny-1,2:nzg-1) )
       Call MPI_Reduce (local_sum,maxV,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

       local_sum = Maxval ( W(2:nxg-1,2:nyg-1,2:nz-1) )
       Call MPI_Reduce (local_sum,maxW,1,MPI_real8,MPI_max,0,MPI_COMM_WORLD,ierr)

       Call check_divergence(max_divergence)

       ! end measure time per step
       time2 = MPI_WTIME()

       ! processor 0 shows the results
       If ( myid==0 ) Then

          Write(*,*) 'step number :', istep
          Write(*,*) 'time        :', t
          Write(*,*) 'time step   :', dt
          
          Write(*,*) ' '
          Write(*,*) 'Retau:      :', Retau

          Write(*,*) ' '          
          Write(*,*) 'Maximum U   :', maxU
          Write(*,*) 'Maximum V   :', maxV
          Write(*,*) 'Maximum W   :', maxW
          
          Write(*,*) 'Mean U      :', meanU/Real( nxm_global*nym_global*nzm_global, 8 )
          Write(*,*) 'Mean V      :', meanV/Real( nxm_global*nym_global*nzm_global, 8 )
          Write(*,*) 'Mean W      :', meanW/Real( nxm_global*nym_global*nzm_global, 8 )

          Write(*,*) ' '
          Write(*,*) 'Mean mass flow in x         :', Qflow_x
          Write(*,*) 'Mean mass flow in y         :', Qflow_y
          Write(*,*) 'Mean pressure gradient in x :', dPdx
          Write(*,*) 'Mean pressure gradient in y :', dPdy
          
          Write(*,*) ' '
          write(*,*) 'Maximum divergence          :', max_divergence
          write(*,*) 'Elapsed time (s)            :', time2-time1
          
          Write(*,*) '------------------------------------------------------'

       End If

       ! start measure time per step
       time1 = MPI_WTIME()

       Call Mpi_barrier(MPI_COMM_WORLD,ierr)

    End If

  End Subroutine output_monitor

  !------------------------------------------!
  !   Output summary of initial parameters   !
  !------------------------------------------!
  Subroutine summary

    If ( myid==0 ) Then

       Write(*,*) '------------------------------------------------------------'
       Write(*,*) '              Summary of initial parameters                 '
       Write(*,*) ' '

       !Write(*,*) 'Note: rhs interpolated to good points'
       Write(*,*) ' '
       Write(*,*) 'Number of processors:',nprocs

       Write(*,*) ' '
       Write(*,*) 'filein:  ',Trim(filein)
       Write(*,*) 'fileout: ',Trim(fileout)
       
       Write(*,*) ' '
       Write(*,*) 'nu     :', nu
       Write(*,*) 'CFL    :', CFL

       Write(*,*) ' '
       If ( x_mass_cte==1 ) Then
          Write(*,*) 'Constant mass flow in x'
       Else
          Write(*,*) 'dPdx    :', dPdx
       End If
       If ( y_mass_cte==1 ) Then
          Write(*,*) 'Constant mass flow in y'
       End If
       Write(*,*) 'dPdz    :', dPdz
       
       Write(*,*) ' '
       Write(*,*) 'nsteps   :', nsteps
       Write(*,*) 'nsave    :', nsave
       Write(*,*) 'nstats   :', nstats
       Write(*,*) 'nmonitor :', nmonitor

       Write(*,*) ' '       
       Write(*,*) 'Lx,Ly,Lz:'
       Write(*,*) Lx,Ly,Lz
       
       Write(*,*) ' '
       Write(*,*) 'nx,nxg,nxm :',nx_global,nxg_global,nxm_global
       Write(*,*) 'ny,nyg,nym :',ny_global,nyg_global,nym_global
       Write(*,*) 'nz,nzg,nzm :',nz_global,nzg_global,nzm_global
       Write(*,*) 'nz,nzg,nzm (local) :',nz,nzg,nzm
       
       Write(*,*) ' '
       Write(*,*) 'xg(1),xg(end) :',xg_global(1), xg_global(nxg_global)
       Write(*,*) 'x (1),x (end) :', x_global(1), x_global ( nx_global)
       
       Write(*,*) 'yg(1),yg(end) :',yg_global(1), yg_global(nyg_global)
       Write(*,*) 'y (1),y (end) :', y_global(1), y_global ( ny_global)
       
       Write(*,*) 'zg(1),xz(end) :',zg_global(1), zg_global(nzg_global)
       Write(*,*) 'z (1),z (end) :', z_global(1), z_global ( nz_global)
       
       Write(*,*) ' '
       Write(*,*) '------------------------------------------------------------'
    
    End If
    
  End Subroutine summary

End Module monitor

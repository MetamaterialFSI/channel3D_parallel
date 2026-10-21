!--------------------------------------------!
! Module to monitor status of the simulation !
!--------------------------------------------!
Module monitor

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use immersed_boundary_operators
  Use differential_operators

  ! prevent implicit typing
  Implicit None

Contains

  !------------------------------------------!
  ! Output some key values during simulation !
  !------------------------------------------!
  Subroutine output_monitor

    Real(Int64) ::  maxU,  maxV,  maxW, local_sum
    Real(Int64) :: meanU, meanV, meanW
    Real(Int64) :: max_divergence, continuity_max, continuity_rms
    Real(Int64), Allocatable :: wall_error(:,:)
    Integer(Int32) :: body, nsolves
    Logical :: report_physical_constraints

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

      nsolves = cg_nsolves
      report_physical_constraints = trim(body_type) /= 'none' .And. nsolves > 0
      Call check_divergence(max_divergence, continuity_max, continuity_rms, report_physical_constraints)
      
      If ( trim(body_type) /= 'none' ) Then
        Allocate(wall_error(4, 0:nbodies))
        Call check_wall_errors(wall_error)
      End If

      ! compute average number of bicgstab iterations over the solves that
      ! actually ran since the last monitor output, rather than assuming a
      ! fixed number of solves per time step
      If ( nsolves > 0 ) Then
        cg_mean_iter = Real(cg_accum_iter, 8) / Real(nsolves, 8)
      Else
        cg_mean_iter = 0d0
      End If
      cg_accum_iter = 0
      cg_nsolves    = 0

      ! end measure time per step
      time2 = MPI_WTIME()

      ! processor 0 shows the results
      If ( myid==0 ) Then

        Write(*,*) 'step number :', istep
        Write(*,*) 'time        :', t
        Write(*,*) 'time step   :', dt
        
        Write(*,*) ' '
        Write(*,*) 'Retau_u:      :', Retau_u
        Write(*,*) 'Retau_w:      :', Retau_w

        Write(*,*) ' '          
        Write(*,*) 'Maximum U   :', maxU
        Write(*,*) 'Maximum V   :', maxV
        Write(*,*) 'Maximum W   :', maxW
        
        Write(*,*) 'Mean U      :', meanU/Real( nxm_global*nym_global*nzm_global, 8 )
        Write(*,*) 'Mean V      :', meanV/Real( nxm_global*nym_global*nzm_global, 8 )
        Write(*,*) 'Mean W      :', meanW/Real( nxm_global*nym_global*nzm_global, 8 )

        Write(*,*) ' '
        Write(*,*) 'Mean masked mass flow in x  :', Qflow_x
        Write(*,*) 'Mean masked mass flow in y  :', Qflow_y
        Write(*,*) 'Mean masked mass flow in z  :', Qflow_z
        Write(*,*) 'Mean pressure gradient in x :', dPdx
        Write(*,*) 'Mean pressure gradient in y :', dPdy
        Write(*,*) 'Mean pressure gradient in z :', dPdz
        
        Write(*,*) ' '
        write(*,*) 'Maximum divergence               :', max_divergence
        If ( trim(body_type) /= 'none' ) Then
          If (report_physical_constraints) Then
            Do body = 0, nbodies
              Write(*,'(A,2(1X,I0),4(1X,ES16.8))') 'IB_WALL_ERROR', istep, body, wall_error(:, body)
            End Do
            Write(*,'(A,1X,I0,2(1X,ES16.8))') 'IB_CONTINUITY_ERROR', istep, &
              continuity_max, continuity_rms
          End If
          If ( nsolves > 0 ) Then
            write(*,'(A,F10.3)') ' Average BiCGSTAB iteration count :', cg_mean_iter
            Call output_bicgstab_diagnostics
          End If
        End If
        write(*,*) 'Elapsed time (s)                 :', time2-time1
        
        Write(*,*) '------------------------------------------------------'

      End If

      If (Allocated(wall_error)) Deallocate(wall_error)

      cg_report_count = 0

      ! start measure time per step
      time1 = MPI_WTIME()

      Call Mpi_barrier(MPI_COMM_WORLD,ierr)

    End If

  End Subroutine output_monitor

  Subroutine output_bicgstab_diagnostics
    Integer(Int32) :: solve_index, iteration, block
    Real(Int64) :: maximum_scaled_residual, relative_to_initial(2)

    Do solve_index = 1, cg_report_count
      Write(*,'(A,3(1X,I0),5(1X,ES16.8))') 'BICGSTAB_SOLVE', istep, &
        cg_report_stages(solve_index), cg_report_iterations(solve_index), &
        cg_report_rhs_rms(:, solve_index), cg_report_tolerances(:, solve_index)
      Do iteration = 0, cg_report_iterations(solve_index)
        maximum_scaled_residual = Maxval(cg_report_history(4:5, iteration, solve_index))
        relative_to_initial = 0d0
        Do block = 1, 2
          If (cg_report_history(3 + block, 0, solve_index) > 0d0) Then
            relative_to_initial(block) = cg_report_history(3 + block, iteration, solve_index) / &
              cg_report_history(3 + block, 0, solve_index)
          Else If (cg_report_history(3 + block, iteration, solve_index) > 0d0) Then
            relative_to_initial(block) = -1d0
          End If
        End Do
        Write(*,'(A,3(1X,I0),8(1X,ES16.8))') 'BICGSTAB_HISTORY', istep, &
          cg_report_stages(solve_index), iteration, cg_report_history(:, iteration, solve_index), &
          maximum_scaled_residual, relative_to_initial
      End Do
    End Do
  End Subroutine output_bicgstab_diagnostics

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

      Write(*,*) 'Input parameters file:',Trim(fileparams)

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
      If ( z_mass_cte==1 ) Then
         Write(*,*) 'Constant mass flow in z'
      Else
        Write(*,*) 'dPdz    :', dPdz
      End If
      
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
      
      Write(*,*) 'zg(1),zg(end) :',zg_global(1), zg_global(nzg_global)
      Write(*,*) 'z (1),z (end) :', z_global(1), z_global ( nz_global)
      
      Write(*,*) ' '
      Write(*,*) 'dxmin,dxmax :',dxmin, Maxval ( xg_global(2:nxg_global) - xg_global(1:nxg_global-1) ) 
      Write(*,*) 'dymin,dymax :',dymin, Maxval ( yg_global(2:nyg_global) - yg_global(1:nyg_global-1) ) 
      Write(*,*) 'dzmin,dzmax :',dzmin, Maxval ( zg_global(2:nzg_global) - zg_global(1:nzg_global-1) ) 

      Write(*,*) ' '
      Write(*,*) 'nxb :', nxb
      Write(*,*) 'nzb :', nzb
      Write(*,*) 'nb  :', nb

      Write(*,*) ' '
      Write(*,*) 'dxb / dx :', dxb / dx
      Write(*,*) 'dzb / dz :', dzb / dz

      ! xb/yb/zb hold the Lagrangian boundary points; there are none (nb = 0)
      ! when body_type = 'none', so indexing them would be out of bounds.
      If ( trim(body_type) /= 'none' ) Then
        Write(*,*) ' '
        Write(*,*) 'xb(1),xb(nb) :',xb(1), xb(nb)
        Write(*,*) 'yb(1),yb(nb) :',yb(1), yb(nb)
        Write(*,*) 'zb(1),zb(nb) :',zb(1), zb(nb)
      End If

      Write(*,*) ' '
      Write(*,*) '------------------------------------------------------------'
    
    End If
    
  End Subroutine summary

  !----------------------------------------------------!
  !   Compute maximum divergence for interior points   !
  !                                                    !
  ! Input:  U, V, W                                    !
  ! Output: div, max_divergence                        !
  !                                                    !
  !----------------------------------------------------!
  Subroutine check_divergence(max_divergence, continuity_max, continuity_rms, compute_continuity)

    Real(Int64), Intent(Out) :: max_divergence, continuity_max, continuity_rms
    Logical, Intent(In) :: compute_continuity
    Integer(Int32) :: i, j, k
    Real(Int64) :: max_divergence_local
    Real(Int64) :: squared_error_local, volume_local, squared_error, volume, cell_volume

    rhs_p = 0d0
    Call divergence(rhs_p, U, V, W)
    max_divergence_local = Maxval(Abs(rhs_p(2:nxg-1, 2:nyg-1, 2:nzg-1)))
    Call MPI_Reduce(max_divergence_local, max_divergence, 1, MPI_real8, MPI_max, 0, MPI_COMM_WORLD,ierr)
    continuity_max = 0d0
    continuity_rms = 0d0
    If (compute_continuity) Then
      squared_error_local = 0d0
      volume_local = 0d0
      Do k = 2, nzg - 1
        Do j = 2, nyg - 1
          Do i = 2, nxg - 1
            cell_volume = (x(i) - x(i - 1)) * (y(j) - y(j - 1)) * (z(k) - z(k - 1))
            squared_error_local = squared_error_local + cell_volume * rhs_p(i, j, k)**2
            volume_local = volume_local + cell_volume
          End Do
        End Do
      End Do
      Call MPI_Reduce(squared_error_local, squared_error, 1, MPI_real8, MPI_sum, 0, MPI_COMM_WORLD, ierr)
      Call MPI_Reduce(volume_local, volume, 1, MPI_real8, MPI_sum, 0, MPI_COMM_WORLD, ierr)
      If (myid == 0) Then
        If (volume <= 0d0) Error Stop 'Continuity diagnostic volume must be positive'
        continuity_max = max_divergence
        continuity_rms = Sqrt(squared_error / volume)
      End If
    End If

  End Subroutine check_divergence

  Subroutine check_wall_errors(wall_error)

    Real(Int64), Intent(Out) :: wall_error(4, 0:nbodies)
    Integer(Int32) :: body, k
    Real(Int64) :: error_x, error_y, error_z
    Real(Int64) :: normal_error, tangent_1_error, tangent_2_error, tangent_squared
    Real(Int64) :: area(0:nbodies), normal_squared(0:nbodies), tangent_sum(0:nbodies)

    Call regT(Eu, U, V, W)
    aux_surface_vector = Eu - ub

    wall_error = 0d0
    area = 0d0
    normal_squared = 0d0
    tangent_sum = 0d0

    Do k = 1, nb
      body = marker_body_index(k)
      error_x = aux_surface_vector(k)
      error_y = aux_surface_vector(nb + k)
      error_z = aux_surface_vector(2 * nb + k)
      normal_error = error_x * normals(k) + error_y * normals(nb + k) + &
        error_z * normals(2 * nb + k)
      tangent_1_error = error_x * tangents_1(k) + error_y * tangents_1(nb + k) + &
        error_z * tangents_1(2 * nb + k)
      tangent_2_error = error_x * tangents_2(k) + error_y * tangents_2(nb + k) + &
        error_z * tangents_2(2 * nb + k)
      tangent_squared = tangent_1_error**2 + tangent_2_error**2

      Call accumulate_wall_error(0, sb(k), normal_error, tangent_squared, wall_error, area, &
        normal_squared, tangent_sum)
      Call accumulate_wall_error(body, sb(k), normal_error, tangent_squared, wall_error, area, &
        normal_squared, tangent_sum)
    End Do

    Do body = 0, nbodies
      If (area(body) <= 0d0) Error Stop 'Wall-error diagnostic area must be positive'
      wall_error(2, body) = Sqrt(normal_squared(body) / area(body))
      wall_error(4, body) = Sqrt(tangent_sum(body) / area(body))
    End Do

  End Subroutine check_wall_errors

  Subroutine accumulate_wall_error(body, marker_area, normal_error, tangent_squared, wall_error, &
      area, normal_squared, tangent_sum)

    Integer(Int32), Intent(In) :: body
    Real(Int64), Intent(In) :: marker_area, normal_error, tangent_squared
    Real(Int64), Intent(InOut) :: wall_error(4, 0:nbodies)
    Real(Int64), Intent(InOut) :: area(0:nbodies), normal_squared(0:nbodies), tangent_sum(0:nbodies)

    area(body) = area(body) + marker_area
    wall_error(1, body) = Max(wall_error(1, body), Abs(normal_error))
    wall_error(3, body) = Max(wall_error(3, body), Sqrt(tangent_squared))
    normal_squared(body) = normal_squared(body) + marker_area * normal_error**2
    tangent_sum(body) = tangent_sum(body) + marker_area * tangent_squared

  End Subroutine accumulate_wall_error

  Integer(Int32) Function marker_body_index(marker) Result(body)

    Integer(Int32), Intent(In) :: marker
    Integer(Int32) :: nxb1, within_ring

    Select Case (trim(body_type))
      Case ('double_cylinders_z')
        nxb1 = Int(2 * 3.14159 * body_param_1 / dxb)
        within_ring = Modulo(marker - 1, nxb) + 1
        If (within_ring <= nxb1) Then
          body = 1
        Else
          body = 2
        End If
      Case ('standing_wave_x', 'traveling_wave_x', 'traveling_wave_z')
        within_ring = Modulo(marker - 1, 2 * nxb) + 1
        If (within_ring <= nxb) Then
          body = 1
        Else
          body = 2
        End If
      Case Default
        If (nbodies /= 1) Error Stop 'Unknown multi-body marker layout'
        body = 1
    End Select

  End Function marker_body_index

End Module monitor

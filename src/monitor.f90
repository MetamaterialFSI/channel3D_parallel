!--------------------------------------------!
! Module to monitor status of the simulation !
!--------------------------------------------!
Module monitor

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use immersed_boundary_operators
  Use immersed_boundary_geometry, Only : compute_average_pressure_per_body
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
    Real(Int64) :: max_divergence
    Real(Int64) :: corrected_divergence_max, corrected_divergence_rms
    Real(Int64), Allocatable :: raw_wall_error(:,:), corrected_wall_error(:,:), gauge_mean(:)
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
      Call check_divergence(max_divergence, corrected_divergence_max, corrected_divergence_rms, &
        report_physical_constraints)
      
      If ( trim(body_type) /= 'none' ) Then
        Allocate(raw_wall_error(4, 0:nbodies), corrected_wall_error(4, 0:nbodies), gauge_mean(nbodies))
        Call check_wall_errors(raw_wall_error, corrected_wall_error)
        If (report_physical_constraints) Then
          Call compute_average_pressure_per_body(gauge_mean, fb(2 * nb + 1 : 3 * nb))
        End If
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
              Write(*,'(A,2(1X,I0),8(1X,ES16.8))') 'IB_WALL_ERROR', istep, body, &
                raw_wall_error(:, body), corrected_wall_error(:, body)
            End Do
            Write(*,'(A,1X,I0,2(1X,ES16.8))') 'IB_CONTINUITY_ERROR', istep, &
              corrected_divergence_max, corrected_divergence_rms
            Do body = 1, nbodies
              Write(*,'(A,2(1X,I0),2(1X,ES16.8))') 'IB_GAUGE', istep, body, gauge_mean(body), &
                -pjv_gscale * gauge_mean(body)
            End Do
          End If
          If ( nsolves > 0 ) Then
            write(*,'(A,F10.3)') ' Average BiCGSTAB iteration count :', cg_mean_iter
            Call output_bicgstab_diagnostics
          End If
        End If
        write(*,*) 'Elapsed time (s)                 :', time2-time1
        
        Write(*,*) '------------------------------------------------------'

      End If

      If (Allocated(raw_wall_error)) Deallocate(raw_wall_error, corrected_wall_error, gauge_mean)

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
  Subroutine check_divergence(max_divergence, corrected_max, corrected_rms, compute_corrected)

    Real(Int64), Intent(Out) :: max_divergence, corrected_max, corrected_rms
    Logical, Intent(In) :: compute_corrected
    Integer(Int32) :: i, j, k
    Real(Int64) :: max_divergence_local
    Real(Int64) :: corrected_max_local, corrected_squared_local, volume_local
    Real(Int64) :: corrected_squared, volume, cell_volume, error

    rhs_p = 0d0
    Call divergence(rhs_p, U, V, W)
    max_divergence_local = Maxval(Abs(rhs_p(2:nxg-1, 2:nyg-1, 2:nzg-1)))

    Call MPI_Reduce(max_divergence_local, max_divergence, 1, MPI_real8, MPI_max, 0, MPI_COMM_WORLD,ierr)

    corrected_max = 0d0
    corrected_rms = 0d0
    If (compute_corrected) Then
      corrected_max_local = 0d0
      corrected_squared_local = 0d0
      volume_local = 0d0
      Do k = 2, nzg - 1
        Do j = 2, nyg - 1
          Do i = 2, nxg - 1
            error = rhs_p(i, j, k) - con_F(i, j, k)
            cell_volume = (x(i) - x(i - 1)) * (y(j) - y(j - 1)) * (z(k) - z(k - 1))
            corrected_max_local = Max(corrected_max_local, Abs(error))
            corrected_squared_local = corrected_squared_local + cell_volume * error**2
            volume_local = volume_local + cell_volume
          End Do
        End Do
      End Do
      Call MPI_Reduce(corrected_max_local, corrected_max, 1, MPI_real8, MPI_max, 0, MPI_COMM_WORLD, ierr)
      Call MPI_Reduce(corrected_squared_local, corrected_squared, 1, MPI_real8, MPI_sum, 0, MPI_COMM_WORLD, ierr)
      Call MPI_Reduce(volume_local, volume, 1, MPI_real8, MPI_sum, 0, MPI_COMM_WORLD, ierr)
      If (myid == 0) Then
        If (volume <= 0d0) Error Stop 'Continuity diagnostic volume must be positive'
        corrected_rms = Sqrt(corrected_squared / volume)
      End If
    End If

  End Subroutine check_divergence

  Subroutine check_wall_errors(raw_error, corrected_error)

    Real(Int64), Intent(Out) :: raw_error(4, 0:nbodies), corrected_error(4, 0:nbodies)
    Integer(Int32) :: body, k
    Real(Int64) :: raw_x, raw_y, raw_z, corrected_x, corrected_y, corrected_z
    Real(Int64) :: raw_normal, raw_tangent_1, raw_tangent_2, raw_tangent_squared
    Real(Int64) :: corrected_normal, corrected_tangent_1, corrected_tangent_2, corrected_tangent_squared
    Real(Int64) :: area(0:nbodies), raw_normal_squared(0:nbodies), raw_tangent_sum(0:nbodies)
    Real(Int64) :: corrected_normal_squared(0:nbodies), corrected_tangent_sum(0:nbodies)

    Call regT(Eu, U, V, W)
    aux_surface_vector = Eu - E1nH_exterior * dudn_jump - ub

    raw_error = 0d0
    corrected_error = 0d0
    area = 0d0
    raw_normal_squared = 0d0
    raw_tangent_sum = 0d0
    corrected_normal_squared = 0d0
    corrected_tangent_sum = 0d0

    Do k = 1, nb
      body = marker_body_index(k)
      raw_x = Eu(k) - ub(k)
      raw_y = Eu(nb + k) - ub(nb + k)
      raw_z = Eu(2 * nb + k) - ub(2 * nb + k)
      corrected_x = aux_surface_vector(k)
      corrected_y = aux_surface_vector(nb + k)
      corrected_z = aux_surface_vector(2 * nb + k)

      raw_normal = raw_x * normals(k) + raw_y * normals(nb + k) + raw_z * normals(2 * nb + k)
      raw_tangent_1 = raw_x * tangents_1(k) + raw_y * tangents_1(nb + k) + &
        raw_z * tangents_1(2 * nb + k)
      raw_tangent_2 = raw_x * tangents_2(k) + raw_y * tangents_2(nb + k) + &
        raw_z * tangents_2(2 * nb + k)
      raw_tangent_squared = raw_tangent_1**2 + raw_tangent_2**2

      corrected_normal = corrected_x * normals(k) + corrected_y * normals(nb + k) + &
        corrected_z * normals(2 * nb + k)
      corrected_tangent_1 = corrected_x * tangents_1(k) + corrected_y * tangents_1(nb + k) + &
        corrected_z * tangents_1(2 * nb + k)
      corrected_tangent_2 = corrected_x * tangents_2(k) + corrected_y * tangents_2(nb + k) + &
        corrected_z * tangents_2(2 * nb + k)
      corrected_tangent_squared = corrected_tangent_1**2 + corrected_tangent_2**2

      Call accumulate_wall_error(0, sb(k), raw_normal, raw_tangent_squared, corrected_normal, &
        corrected_tangent_squared, raw_error, corrected_error, area, raw_normal_squared, &
        raw_tangent_sum, corrected_normal_squared, corrected_tangent_sum)
      Call accumulate_wall_error(body, sb(k), raw_normal, raw_tangent_squared, corrected_normal, &
        corrected_tangent_squared, raw_error, corrected_error, area, raw_normal_squared, &
        raw_tangent_sum, corrected_normal_squared, corrected_tangent_sum)
    End Do

    Do body = 0, nbodies
      If (area(body) <= 0d0) Error Stop 'Wall-error diagnostic area must be positive'
      raw_error(2, body) = Sqrt(raw_normal_squared(body) / area(body))
      raw_error(4, body) = Sqrt(raw_tangent_sum(body) / area(body))
      corrected_error(2, body) = Sqrt(corrected_normal_squared(body) / area(body))
      corrected_error(4, body) = Sqrt(corrected_tangent_sum(body) / area(body))
    End Do

  End Subroutine check_wall_errors

  Subroutine accumulate_wall_error(body, marker_area, raw_normal, raw_tangent_squared, corrected_normal, &
      corrected_tangent_squared, raw_error, corrected_error, area, raw_normal_squared, raw_tangent_sum, &
      corrected_normal_squared, corrected_tangent_sum)

    Integer(Int32), Intent(In) :: body
    Real(Int64), Intent(In) :: marker_area, raw_normal, raw_tangent_squared
    Real(Int64), Intent(In) :: corrected_normal, corrected_tangent_squared
    Real(Int64), Intent(InOut) :: raw_error(4, 0:nbodies), corrected_error(4, 0:nbodies)
    Real(Int64), Intent(InOut) :: area(0:nbodies), raw_normal_squared(0:nbodies), raw_tangent_sum(0:nbodies)
    Real(Int64), Intent(InOut) :: corrected_normal_squared(0:nbodies), corrected_tangent_sum(0:nbodies)

    area(body) = area(body) + marker_area
    raw_error(1, body) = Max(raw_error(1, body), Abs(raw_normal))
    raw_error(3, body) = Max(raw_error(3, body), Sqrt(raw_tangent_squared))
    corrected_error(1, body) = Max(corrected_error(1, body), Abs(corrected_normal))
    corrected_error(3, body) = Max(corrected_error(3, body), Sqrt(corrected_tangent_squared))
    raw_normal_squared(body) = raw_normal_squared(body) + marker_area * raw_normal**2
    raw_tangent_sum(body) = raw_tangent_sum(body) + marker_area * raw_tangent_squared
    corrected_normal_squared(body) = corrected_normal_squared(body) + marker_area * corrected_normal**2
    corrected_tangent_sum(body) = corrected_tangent_sum(body) + marker_area * corrected_tangent_squared

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

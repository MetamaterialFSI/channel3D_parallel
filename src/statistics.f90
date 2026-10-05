!--------------------------------------------!
! Module for computing some basic statistics !
!--------------------------------------------!
Module statistics

  ! Modules
  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global
  Use mpi
  Use interpolation
  Use input_output
  Use mass_flow

  ! prevent implicit typing
  Implicit None

Contains

  !------------------------------------------!
  ! Compute some basic statistics on the fly !
  !------------------------------------------!
  Subroutine compute_statistics    

    Integer(Int32) :: jj
    Real   (Int64) :: dUdy_wall_b, dUdy_wall_t
    Real   (Int64) :: dWdy_wall_b, dWdy_wall_t
    Real   (Int64) ::   UV_wall_b,   UV_wall_t
    Real   (Int64) ::   VW_wall_b,   VW_wall_t

    ! if pressure not computed 
    pressure_computed = .False.

    ! statistics computed at grid y -> U and W interpolated    
    if ( Mod(istep,nstats)==0 .Or. istep==1 ) Then

       ! Compute actual pressure (should be called first, uses term_1,...)
       !Call compute_pressure       
       ! now computed in projection.f90
       pressure_computed = .True.
       
       ! interpolate U in x -> term_1
       Call interpolate_x(U,term_1(2:nxg-1,1:nyg,1:nzg))

       ! interpolate V in y -> term_2
       Call interpolate_y(V,term_2(1:nxg,2:nyg-1,1:nzg))
       ! boundary condition for V
       term_2(:,  1,:) = -term_2(:,    2,:)
       term_2(:,nyg,:) = -term_2(:,nyg-1,:)

       ! compute local statistics
       If ( myid < nprocs-1 ) Then
          Do jj=1,nyg
             Umean  (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-1) )
             Vmean  (jj) = Sum( term_2(2:nxg-2,jj,2:nzg-1) )
             Wmean  (jj) = Sum(      W(2:nxg-2,jj,2:nz -1) )

             U2mean (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-1)**2d0 )
             V2mean (jj) = Sum( term_2(2:nxg-2,jj,2:nzg-1)**2d0 )
             W2mean (jj) = Sum(      W(2:nxg-2,jj,2:nz -1)**2d0 )

             UVmean (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-1)*term_2(2:nxg-2,jj,2:nzg-1) )
          End Do
        Else
          Do jj=1,nyg
             Umean  (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-2) )
             Vmean  (jj) = Sum( term_2(2:nxg-2,jj,2:nzg-2) )
             Wmean  (jj) = Sum(      W(2:nxg-2,jj,2:nz -1) )

             U2mean (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-2)**2d0 )
             V2mean (jj) = Sum( term_2(2:nxg-2,jj,2:nzg-2)**2d0 )
             W2mean (jj) = Sum(      W(2:nxg-2,jj,2:nz -1)**2d0 )

             UVmean (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-2)*term_2(2:nxg-2,jj,2:nzg-2) )
          End Do
        End If

       ! reduce statatistics between processors      
       IF ( myid==0 ) Then

          Call MPI_Reduce(MPI_IN_PLACE,Umean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,Vmean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,Wmean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(MPI_IN_PLACE,U2mean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,V2mean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(MPI_IN_PLACE,W2mean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(MPI_IN_PLACE,UVmean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       Else

          Call MPI_Reduce(Umean,0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(Vmean,0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(Wmean,0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(U2mean,0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(V2mean,0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          Call MPI_Reduce(W2mean,0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
          
          Call MPI_Reduce(UVmean,0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
       End If

       ! These statistics are only good for processor 0
       Umean  = Umean/Real( (nxg_global-3)*(nzg_global-3), 8)
       Vmean  = Vmean/Real( (nxg_global-3)*(nzg_global-3), 8)
       Wmean  = Wmean/Real( (nxg_global-3)*( nz_global-2), 8)

       U2mean = U2mean/Real( (nxg_global-3)*(nzg_global-3), 8)
       V2mean = V2mean/Real( (nxg_global-3)*(nzg_global-3), 8)
       W2mean = W2mean/Real( (nxg_global-3)*( nz_global-2), 8)

       UVmean = UVmean/Real( (nxg_global-3)*(nzg_global-3) ,8)

       ! Mean derivative at the walls (CHECK THIS PLEASE)
       dUdy_wall_b = ( Umean(2) -  Umean(1) )/( yg(  2) - yg(    1)) 
       dUdy_wall_t = ( Umean(nyg) - Umean(nyg-1) )/( yg(nyg) - yg(nyg-1))
       dWdy_wall_b = ( Wmean(2) -  Wmean(1) )/( yg(  2) - yg(    1)) 
       dWdy_wall_t = ( Wmean(nyg) - Wmean(nyg-1) )/( yg(nyg) - yg(nyg-1))

       ! Mean Reynolds stress at the walls
       UV_wall_b = 0d0 
       UV_wall_t = 0d0 
       VW_wall_b = 0d0 
       VW_wall_t = 0d0 

       ! friction velocity
       utau = ( ( UV_wall_t - UV_wall_b - nu*dUdy_wall_t + nu*dUdy_wall_b )/Ly )**0.5d0
       !wtau = ( ( VW_wall_t - VW_wall_b - nu*dWdy_wall_t + nu*dWdy_wall_b )/Ly )**0.5d0
       wtau = Sqrt(Max(0.0d0, &
       (VW_wall_t - VW_wall_b - nu*dWdy_wall_t + nu*dWdy_wall_b)/Ly))
       ! friction Reynolds number
       Retau_u = utau*(y(ny)-y(1))/2d0/nu
       Retau_w = wtau*(y(ny)-y(1))/2d0/nu

       ! mean mass flow in x
       Call compute_mean_mass_flow_U(U, Qflow_x)
       Call compute_mean_mass_flow_V(V, Qflow_y)
       Call compute_mean_mass_flow_W(W, Qflow_z)

       ! write statistics
       Call output_statistics

       ! Sanity check 
       If ( Any( Isnan(U) ) ) Stop 'Error: NaNs!'
       If ( Any( Isnan(V) ) ) Stop 'Error: NaNs!'
       If ( Any( Isnan(W) ) ) Stop 'Error: NaNs!'
       
    End If
    
  End Subroutine compute_statistics
Subroutine compute_statistics_topchannel    

  Integer(Int32) :: jj
  Integer(Int32) :: jIB
  Real(Int64)    :: dUdy_wall_b, dUdy_wall_t
  Real(Int64)    :: dWdy_wall_b, dWdy_wall_t
  Real(Int64)    :: UV_wall_b, UV_wall_t
  Real(Int64)    :: VW_wall_b, VW_wall_t
  Real(Int64)    :: Ly_top, yIB

  yIB = 2.d0

  pressure_computed = .False.

  If ( Mod(istep,nstats)==0 .Or. istep==1 ) Then

     pressure_computed = .True.

     Call interpolate_x(U,term_1(2:nxg-1,1:nyg,1:nzg))

     Call interpolate_y(V,term_2(1:nxg,2:nyg-1,1:nzg))

     term_2(:,  1,:) = -term_2(:,    2,:)
     term_2(:,nyg,:) = -term_2(:,nyg-1,:)

     ! local statistics
     If ( myid < nprocs-1 ) Then
        Do jj=1,nyg
           Umean  (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-1) )
           Vmean  (jj) = Sum( term_2(2:nxg-2,jj,2:nzg-1) )
           Wmean  (jj) = Sum(      W(2:nxg-2,jj,2:nz -1) )

           U2mean (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-1)**2d0 )
           V2mean (jj) = Sum( term_2(2:nxg-2,jj,2:nzg-1)**2d0 )
           W2mean (jj) = Sum(      W(2:nxg-2,jj,2:nz -1)**2d0 )

           UVmean (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-1) * &
                              term_2(2:nxg-2,jj,2:nzg-1) )
        End Do
     Else
        Do jj=1,nyg
           Umean  (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-2) )
           Vmean  (jj) = Sum( term_2(2:nxg-2,jj,2:nzg-2) )
           Wmean  (jj) = Sum(      W(2:nxg-2,jj,2:nz -1) )

           U2mean (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-2)**2d0 )
           V2mean (jj) = Sum( term_2(2:nxg-2,jj,2:nzg-2)**2d0 )
           W2mean (jj) = Sum(      W(2:nxg-2,jj,2:nz -1)**2d0 )

           UVmean (jj) = Sum( term_1(2:nxg-2,jj,2:nzg-2) * &
                              term_2(2:nxg-2,jj,2:nzg-2) )
        End Do
     End If

     ! reduce statistics
     If ( myid==0 ) Then
        Call MPI_Reduce(MPI_IN_PLACE,Umean, nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
        Call MPI_Reduce(MPI_IN_PLACE,Vmean, nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
        Call MPI_Reduce(MPI_IN_PLACE,Wmean, nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

        Call MPI_Reduce(MPI_IN_PLACE,U2mean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
        Call MPI_Reduce(MPI_IN_PLACE,V2mean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
        Call MPI_Reduce(MPI_IN_PLACE,W2mean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

        Call MPI_Reduce(MPI_IN_PLACE,UVmean,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
     Else
        Call MPI_Reduce(Umean,  0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
        Call MPI_Reduce(Vmean,  0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
        Call MPI_Reduce(Wmean,  0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

        Call MPI_Reduce(U2mean, 0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
        Call MPI_Reduce(V2mean, 0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
        Call MPI_Reduce(W2mean, 0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)

        Call MPI_Reduce(UVmean, 0,nyg,MPI_real8,MPI_sum,0,MPI_COMM_WORLD,ierr)
     End If

     ! normalize
     Umean  = Umean /Real( (nxg_global-3)*(nzg_global-3), 8)
     Vmean  = Vmean /Real( (nxg_global-3)*(nzg_global-3), 8)
     Wmean  = Wmean /Real( (nxg_global-3)*( nz_global-2), 8)

     U2mean = U2mean/Real( (nxg_global-3)*(nzg_global-3), 8)
     V2mean = V2mean/Real( (nxg_global-3)*(nzg_global-3), 8)
     W2mean = W2mean/Real( (nxg_global-3)*( nz_global-2), 8)

     UVmean = UVmean/Real( (nxg_global-3)*(nzg_global-3), 8)

     ! first grid point above IB wall y = 2
     jIB = 2
     Do jj = 2, nyg-1
        If ( yg(jj) > yIB ) Then
           jIB = jj
           Exit
        End If
     End Do

     ! top-channel lower wall = IB wall side
     ! top-channel upper wall = physical top wall
     dUdy_wall_b = ( Umean(jIB+1) - Umean(jIB) ) / &
                   ( yg(jIB+1)  - yg(jIB)  )

     dUdy_wall_t = ( Umean(nyg)  - Umean(nyg-1) ) / &
                   ( yg(nyg)    - yg(nyg-1) )

     dWdy_wall_b = ( Wmean(jIB+1) - Wmean(jIB) ) / &
                   ( yg(jIB+1)  - yg(jIB)  )

     dWdy_wall_t = ( Wmean(nyg)  - Wmean(nyg-1) ) / &
                   ( yg(nyg)    - yg(nyg-1) )

     ! Reynolds stress contribution
     UV_wall_b = UVmean(jIB)
     UV_wall_t = UVmean(nyg)

     VW_wall_b = 0.d0
     VW_wall_t = 0.d0

     Ly_top = yg(nyg) - yIB

     utau = ( ( UV_wall_t - UV_wall_b - nu*dUdy_wall_t + nu*dUdy_wall_b ) / Ly_top )**0.5d0
     !wtau = ( ( VW_wall_t - VW_wall_b - nu*dWdy_wall_t + nu*dWdy_wall_b ) / Ly_top )**0.5d0
wtau = Sqrt(Max(0.0d0, &
       (VW_wall_t - VW_wall_b - nu*dWdy_wall_t + nu*dWdy_wall_b) / Ly_top))
     Retau_u = utau*Ly_top/2.d0/nu
     Retau_w = wtau*Ly_top/2.d0/nu

     Call compute_mean_mass_flow_U(U, Qflow_x)
     Call compute_mean_mass_flow_V(V, Qflow_y)
     Call compute_mean_mass_flow_V(W, Qflow_z)

     Call output_statistics

     If ( Any( Isnan(U) ) ) Stop 'Error: NaNs!'
     If ( Any( Isnan(V) ) ) Stop 'Error: NaNs!'
     If ( Any( Isnan(W) ) ) Stop 'Error: NaNs!'

  End If

End Subroutine compute_statistics_topchannel

Subroutine compute_statistics_topchannel_updated

  Integer(Int32) :: jj
  Integer(Int32) :: jIB

  Real(Int64) :: dUdy_wall_b, dUdy_wall_t
  Real(Int64) :: dWdy_wall_b, dWdy_wall_t

  Real(Int64) :: UV_wall_b, UV_wall_t
  Real(Int64) :: VW_wall_b, VW_wall_t

  Real(Int64) :: Ly_top, yIB

  ! Pressure-gradient-equivalent friction quantities
  Real(Int64) :: Gx_inst
  Real(Int64) :: h_ref

  yIB = 2.d0

  pressure_computed = .False.

  If (Mod(istep,nstats) == 0 .Or. istep == 1) Then

     pressure_computed = .True.

     Call interpolate_x( &
          U, term_1(2:nxg-1,1:nyg,1:nzg) )

     Call interpolate_y( &
          V, term_2(1:nxg,2:nyg-1,1:nzg) )

     term_2(:,1,:)   = -term_2(:,2,:)
     term_2(:,nyg,:) = -term_2(:,nyg-1,:)

     !------------------------------------------------------!
     ! Local statistics.                                    !
     !------------------------------------------------------!

     If (myid < nprocs-1) Then

        Do jj = 1,nyg

           Umean(jj) = Sum( &
                term_1(2:nxg-2,jj,2:nzg-1) )

           Vmean(jj) = Sum( &
                term_2(2:nxg-2,jj,2:nzg-1) )

           Wmean(jj) = Sum( &
                W(2:nxg-2,jj,2:nz-1) )

           U2mean(jj) = Sum( &
                term_1(2:nxg-2,jj,2:nzg-1)**2d0 )

           V2mean(jj) = Sum( &
                term_2(2:nxg-2,jj,2:nzg-1)**2d0 )

           W2mean(jj) = Sum( &
                W(2:nxg-2,jj,2:nz-1)**2d0 )

           UVmean(jj) = Sum( &
                term_1(2:nxg-2,jj,2:nzg-1) * &
                term_2(2:nxg-2,jj,2:nzg-1) )

        End Do

     Else

        Do jj = 1,nyg

           Umean(jj) = Sum( &
                term_1(2:nxg-2,jj,2:nzg-2) )

           Vmean(jj) = Sum( &
                term_2(2:nxg-2,jj,2:nzg-2) )

           Wmean(jj) = Sum( &
                W(2:nxg-2,jj,2:nz-1) )

           U2mean(jj) = Sum( &
                term_1(2:nxg-2,jj,2:nzg-2)**2d0 )

           V2mean(jj) = Sum( &
                term_2(2:nxg-2,jj,2:nzg-2)**2d0 )

           W2mean(jj) = Sum( &
                W(2:nxg-2,jj,2:nz-1)**2d0 )

           UVmean(jj) = Sum( &
                term_1(2:nxg-2,jj,2:nzg-2) * &
                term_2(2:nxg-2,jj,2:nzg-2) )

        End Do

     End If

     !------------------------------------------------------!
     ! Reduce statistics.                                   !
     !------------------------------------------------------!

     If (myid == 0) Then

        Call MPI_Reduce( &
             MPI_IN_PLACE, Umean, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             MPI_IN_PLACE, Vmean, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             MPI_IN_PLACE, Wmean, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             MPI_IN_PLACE, U2mean, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             MPI_IN_PLACE, V2mean, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             MPI_IN_PLACE, W2mean, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             MPI_IN_PLACE, UVmean, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

     Else

        Call MPI_Reduce( &
             Umean, 0, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             Vmean, 0, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             Wmean, 0, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             U2mean, 0, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             V2mean, 0, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             W2mean, 0, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

        Call MPI_Reduce( &
             UVmean, 0, nyg, MPI_real8, MPI_sum, &
             0, MPI_COMM_WORLD, ierr )

     End If

     !------------------------------------------------------!
     ! Normalize plane-averaged statistics.                 !
     !------------------------------------------------------!

     Umean = Umean / Real( &
          (nxg_global-3)*(nzg_global-3), Int64 )

     Vmean = Vmean / Real( &
          (nxg_global-3)*(nzg_global-3), Int64 )

     Wmean = Wmean / Real( &
          (nxg_global-3)*(nz_global-2), Int64 )

     U2mean = U2mean / Real( &
          (nxg_global-3)*(nzg_global-3), Int64 )

     V2mean = V2mean / Real( &
          (nxg_global-3)*(nzg_global-3), Int64 )

     W2mean = W2mean / Real( &
          (nxg_global-3)*(nz_global-2), Int64 )

     UVmean = UVmean / Real( &
          (nxg_global-3)*(nzg_global-3), Int64 )

     !------------------------------------------------------!
     ! First grid point above nominal IB-wall location.     !
     !                                                       !
     ! This fixed-y location is retained only for the old   !
     ! spanwise diagnostic below. It is no longer used to   !
     ! calculate streamwise friction velocity.              !
     !------------------------------------------------------!

     jIB = 2

     Do jj = 2,nyg-1

        If (yg(jj) > yIB) Then
           jIB = jj
           Exit
        End If

     End Do

     !------------------------------------------------------!
     ! Existing fixed-plane gradients.                      !
     !                                                       !
     ! These are not used for the new streamwise utau.      !
     !------------------------------------------------------!

     dUdy_wall_b = &
          (Umean(jIB+1) - Umean(jIB)) / &
          (yg(jIB+1) - yg(jIB))

     dUdy_wall_t = &
          (Umean(nyg) - Umean(nyg-1)) / &
          (yg(nyg) - yg(nyg-1))

     dWdy_wall_b = &
          (Wmean(jIB+1) - Wmean(jIB)) / &
          (yg(jIB+1) - yg(jIB))

     dWdy_wall_t = &
          (Wmean(nyg) - Wmean(nyg-1)) / &
          (yg(nyg) - yg(nyg-1))

     UV_wall_b = UVmean(jIB)
     UV_wall_t = UVmean(nyg)

     VW_wall_b = 0.d0
     VW_wall_t = 0.d0

     Ly_top = yg(nyg) - yIB

     !------------------------------------------------------!
     ! Streamwise pressure-gradient-equivalent friction     !
     ! velocity for the deforming top channel.              !
     !                                                       !
     ! Kim-Choi reference channel:                          !
     !   nominal lower wall: y = 2                          !
     !   upper wall:         y = 4                          !
     !   full height:        2                              !
     !   reference half-height h = 1                        !
     !                                                       !
     ! dPdx_RK is the pressure-gradient forcing used during !
     ! the RK stages. dPdx is the final constant-flow-rate  !
     ! correction divided by dt.                            !
     !------------------------------------------------------!

     h_ref = 1.d0

     If (x_mass_cte == 1) Then

        Gx_inst =  dPdx

     Else

        Gx_inst = dPdx

     End If

     utau = Sqrt( &
          Max(0.d0, Gx_inst*h_ref) )

     Retau_u = utau*h_ref/nu

     !------------------------------------------------------!
     ! Existing spanwise diagnostic.                        !
     !                                                       !
     ! This still uses the nominal fixed-y approximation.   !
     !------------------------------------------------------!

     wtau = Sqrt( &
          Max(0.d0, &
          (VW_wall_t - VW_wall_b                  &
          - nu*dWdy_wall_t + nu*dWdy_wall_b)      &
          / Ly_top) )

     Retau_w = wtau*Ly_top/(2.d0*nu)

     !------------------------------------------------------!
     ! Mean mass-flow rates.                                !
     !------------------------------------------------------!

     Call compute_mean_mass_flow_U(U, Qflow_x)
     Call compute_mean_mass_flow_V(V, Qflow_y)
     Call compute_mean_mass_flow_W(W, Qflow_z)

     Call output_statistics

     If (Any(Isnan(U))) Stop 'Error: NaNs in U!'
     If (Any(Isnan(V))) Stop 'Error: NaNs in V!'
     If (Any(Isnan(W))) Stop 'Error: NaNs in W!'

  End If

End Subroutine compute_statistics_topchannel_updated
End Module statistics

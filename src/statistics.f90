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
       wtau = ( ( VW_wall_t - VW_wall_b - nu*dWdy_wall_t + nu*dWdy_wall_b )/Ly )**0.5d0

       ! friction Reynolds number
       Retau_u = utau*(y(ny)-y(1))/2d0/nu
       Retau_w = wtau*(y(ny)-y(1))/2d0/nu

       ! mean mass flow in x
       Call compute_mean_mass_flow_U(U, Qflow_x)
       Call compute_mean_mass_flow_V(V, Qflow_y)
       Call compute_mean_mass_flow_V(W, Qflow_z)

       ! write statistics
       !Call output_statistics
       !Calculat mean shear stress of IB points
       CALL project_force_to_tangent(fb, tangents_1, fb_t1)
       Call compute_global_mean(fb_t1,tau_w)
       if (myid .eq. 0) Then
         !WRITE(*,*) 'write output stats'
         tau_w_log(store_index,1)=t+REAL(nstep_init)*dt
         tau_w_log(store_index,2)=tau_w 
         dPdx_log(store_index,1)=dPdx
         !WRITE(*,*) 't,tau_w=',tau_w_log(store_index,1),tau_w_log(store_index,2)
         if (store_index .eq. 1000 .or. istep .eq. nsteps) then
           Call output_response
           store_index=1
         elseif ( istep==1) then
          Call output_response
          !store_index=store_index+1
          store_index=1
         elseif (Any( Isnan(U)) .or. Any( Isnan(V)) .or. Any( Isnan(W))) then
            Call output_response
         else
           store_index=store_index+1
         end if

       end if

       ! Sanity check 
       If ( Any( Isnan(U) ) ) Stop 'Error: NaNs!'
       If ( Any( Isnan(V) ) ) Stop 'Error: NaNs!'
       If ( Any( Isnan(W) ) ) Stop 'Error: NaNs!'
       
    End If
    
  End Subroutine compute_statistics

  Subroutine project_force_to_tangent(f_, tangents_, fbt_)

   Implicit None
   Real(Int64), Contiguous,Intent(In)  :: f_(:)
   Real(Int64), Contiguous,Intent(In)  :: tangents_(:)
 
   ! Scalar projection at each IB point
   Real(Int64), Contiguous,Intent(Out) :: fbt_(:)
 
   Integer :: i, ix, iy, iz
 
   Do i = 1, nb
 
      ix = 3*(i-1) + 1
      iy = 3*(i-1) + 2
      iz = 3*(i-1) + 3
 
      ! Dot product: fb · tangents_1
      fbt_(i) = f_(ix) * tangents_(ix) &
               + f_(iy) * tangents_(iy) &
               + f_(iz) * tangents_(iz)
 
   End Do
 
 End Subroutine project_force_to_tangent

 Subroutine compute_global_mean(fbt_, mean_)

   Use mpi
   Implicit None
 
   Real(Int64), Contiguous, Intent(In)  :: fbt_(:)
   Real(Int64),               Intent(Out) :: mean_
 
   Integer :: n_global
 
   Real(Int64) :: local_sum
   Real(Int64) :: global_sum
 
   local_sum = Sum(fbt_)
 
   ! Reduce force sum
   If ( myid == 0 ) Then
 
      Call MPI_Reduce(MPI_IN_PLACE, global_sum, 1, MPI_real8, &
                      MPI_sum, 0, MPI_COMM_WORLD, ierr)
   Else
      global_sum = local_sum
 
      Call MPI_Reduce(global_sum, 0, 1, MPI_real8, &
                      MPI_sum, 0, MPI_COMM_WORLD, ierr)
 
   End If
 
   ! Root computes mean
   If ( myid == 0 ) Then
      mean_ = global_sum / Real(nb, Int64)
   End If
 
 End Subroutine compute_global_mean

End Module statistics

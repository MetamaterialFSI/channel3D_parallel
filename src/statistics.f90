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
    Real   (Int64) ::   UV_wall_b,   UV_wall_t, tmp_t

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

       ! Mean Reynolds stress at the walls
       UV_wall_b = 0d0 
       UV_wall_t = 0d0 

       ! friction velocity
       utau = ( ( UV_wall_t - UV_wall_b - nu*dUdy_wall_t + nu*dUdy_wall_b )/Ly )**0.5d0

       ! friction Reynolds number
       Retau = utau*(y(ny)-y(1))/2d0/nu

       ! mean mass flow in x
       Call compute_mean_mass_flow_U(U,Qflow_x)
       Call compute_mean_mass_flow_V(V,Qflow_y)

       ! write statistics
       Call output_statistics

       ! write statistics of speed
       if (myid .eq. 0) Then
         IB_geo =IB_geo/3
         IB_op =IB_op/3
         non_IB_proj=non_IB_proj/3
         E_1st=E_1st/3
         IB_force=IB_force/3
         R_1st=R_1st/3
         D_1st=D_1st/3
         IB_possion=IB_possion/3
         proj_1st=proj_1st/3
         grad_1st=grad_1st/3
         proj_2nd=proj_2nd/3
         apply_bc = apply_bc/3
         !WRITE(*,*) 'write output stats'
         tmp_t=t+REAL(nstep_init)*dt
         !WRITE(*,*) 'tmp_t=',tmp_t
         time_matrix(store_index,1)=tmp_t
         time_matrix(store_index,2)=IB_geo 
         time_matrix(store_index,3)=IB_op 
         time_matrix(store_index,4)=non_IB_proj 
         time_matrix(store_index,5)=E_1st 
         time_matrix(store_index,6)=IB_force 
         time_matrix(store_index,7)=R_1st 
         time_matrix(store_index,8)=D_1st 
         time_matrix(store_index,9)=IB_possion 
         time_matrix(store_index,10)=proj_1st 
         time_matrix(store_index,11)=grad_1st 
         time_matrix(store_index,12)=proj_2nd 
         time_matrix(store_index,13)=apply_bc 
         time_matrix(store_index,14)=RK1_iter 
         time_matrix(store_index,15)=RK2_iter 
         time_matrix(store_index,16)=RK3_iter 
         !WRITE(*,*) 't,tau_w=',tau_w_log(store_index,1),tau_w_log(store_index,2)
         if (store_index .eq. 1000 .or. istep .eq. nsteps) then
           Call output_time
           store_index=1
         elseif ( istep==1) then
          Call output_time
          !store_index=store_index+1
          store_index=1
         else
           store_index=store_index+1
         end if
         IB_geo =0.d0
         IB_op =0.d0
         non_IB_proj=0.d0
         E_1st=0.d0
         IB_force=0.d0 
         R_1st=0.d0
         D_1st=0.d0
         IB_possion=0.d0
         proj_1st=0.d0
         grad_1st=0.d0
         proj_2nd=0.d0
         apply_bc=0.d0
         RK1_iter=0
         RK2_iter=0
         RK3_iter=0
       end if

       ! Sanity check 
       If ( Any( Isnan(U) ) ) Stop 'Error: NaNs!'
       If ( Any( Isnan(V) ) ) Stop 'Error: NaNs!'
       If ( Any( Isnan(W) ) ) Stop 'Error: NaNs!'
       
    End If
    
  End Subroutine compute_statistics

End Module statistics

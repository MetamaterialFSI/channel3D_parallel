MODULE controls

  ! Control types (stored in ctrl%ctrl_type as a character string):
  ! 'spanwise_const_gauss_x': SIMO feedback control (2D pressure - streamwise periodic Gaussian-type forcing)
  !   Input file format (control_input.NNN.inp):
  !   spanwise_const_gauss_x
  !   nx_local
  !   sigma, error
  !   amp
  !   order of actuator
  !   K_A  (ord_ctrl x ord_ctrl matrix, one row per line)
  !   K_B  (ord_ctrl x 1 vector)
  !   K_C  (1 x ord_ctrl vector)
  !   K_D  (1 x 1 scalar)

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  USE global
  Use mpi
  Use ifport
  Use pressure
  Use boundary_conditions
  Use immersed_boundary_operators

  IMPLICIT NONE

  !  for bodies
  INTEGER, PARAMETER :: maxcontrols = 999 ! a large number
  TYPE(control_t) :: ctrl(maxcontrols)

CONTAINS
! *****************************************************************************************

  SUBROUTINE controls_setup_control
    IMPLICIT NONE
    LOGICAL :: readinput, existctrl, readlogs
    INTEGER :: i,j
    Character(200) :: fname
    CHARACTER(3) :: file_num
    Character(8)   :: ext
    REAL(KIND(0.D0)), ALLOCATABLE :: xm_local(:)
    REAL(KIND(0.D0)) :: lambda

    ! look for controls in input directory
    readinput = .TRUE.
    n_control = 0
    Select case (trim(body_type))
    case ('moving_bottom_wall')
      if (myid .eq. 0) then
        DO WHILE (readinput)
          WRITE(file_num,"(I3.3)")  n_control+1
          INQUIRE(file="./control_input."//file_num//".inp",exist=readinput)
          WRITE(*,*) 'read control input...'
          IF (readinput) THEN
            n_control=n_control+1
            ctrl(n_control)%id = n_control
            ! in each control input file
            OPEN(unit=80,file="./control_input."//file_num//".inp",form='formatted',status='old')
            ! read control type keyword into ctrl_type field
            READ(80,*) ctrl(n_control)%ctrl_type
            ctrl(n_control)%xbs = 0.D0 ! dummy number
            ctrl(n_control)%zbs = 0.D0 ! dummy number
            ctrl(n_control)%width = 0.D0 ! dummy number

            IF (TRIM(ctrl(n_control)%ctrl_type) .eq. 'spanwise_const_gauss_x') Then
              WRITE(*,*) '=> Control no.', n_control, &
              "is 2D SISO feedback controller with streamwise-periodic Gaussian-RPMs"
              ! read actuator parameters
              ALLOCATE(ctrl(n_control)%actuators(1))
              READ(80,*) ctrl(n_control)%actuators(1)%nx_local
              ctrl(n_control)%num_act=(nxb)/ctrl(n_control)%actuators(1)%nx_local
              WRITE(*,*) 'nx_local=',ctrl(n_control)%actuators(1)%nx_local,'num_act=',ctrl(n_control)%num_act
              
              do i = 1, 1
                READ(80,*) ctrl(n_control)%actuators(i)%sigma,ctrl(n_control)%actuators(i)%error_actuator
                READ(80,*) ctrl(n_control)%actuators(i)%amp
                WRITE(*,*) 'sigma=',ctrl(n_control)%actuators(i)%sigma,'error=',ctrl(n_control)%actuators(i)%error_actuator
                WRITE(*,*) 'amp=',ctrl(n_control)%actuators(i)%amp
              end do
              ! compute local zero-flux Gaussian shape function
              ALLOCATE(ctrl(n_control)%zero_func_local(ctrl(n_control)%actuators(1)%nx_local))
              ALLOCATE(ctrl(n_control)%dfdx_local(ctrl(n_control)%actuators(1)%nx_local))
              ALLOCATE(xm_local(ctrl(n_control)%actuators(1)%nx_local))
              ctrl(n_control)%zero_func_local=0.d0
              ctrl(n_control)%dfdx_local=0.d0
              xm_local=0.d0
              do i=1, ctrl(n_control)%actuators(1)%nx_local
                xm_local(i)=(i-0.5d0)*dx-ctrl(n_control)%actuators(1)%nx_local/2.0d0*dx
              end do
              lambda=ctrl(n_control)%actuators(1)%nx_local*dx*2.0D0/3.0D0
              do i=1, ctrl(n_control)%actuators(1)%nx_local
                ctrl(n_control)%zero_func_local(i)=COS(2*pi/lambda*xm_local(i))
                ctrl(n_control)%zero_func_local(i)=ctrl(n_control)%zero_func_local(i)*EXP(-xm_local(i)**2.d0/(2*ctrl(n_control)%actuators(1)%sigma**2.0d0))
                ctrl(n_control)%zero_func_local(i)=ctrl(n_control)%zero_func_local(i)-ctrl(n_control)%actuators(1)%error_actuator
                ! compute dfdx for local profile for tangential and normal vector
                ctrl(n_control)%dfdx_local(i)=EXP(-xm_local(i)**2.d0/(2*ctrl(n_control)%actuators(1)%sigma**2.0d0))
                ctrl(n_control)%dfdx_local(i)=ctrl(n_control)%dfdx_local(i)*(-2*pi/lambda*SIN(2*pi/lambda*xm_local(i))&
                -xm_local(i)/ctrl(n_control)%actuators(1)%sigma**2.0d0*COS(2*pi/lambda*xm_local(i)))
                
              end do
              
              ! read controller state-space parameters (K_A, K_B, K_C, K_D)
              READ(80,*) ctrl(n_control)%ord_ctrl
              WRITE(*,*) 'ord_ctrl=',ctrl(n_control)%ord_ctrl
              ctrl(n_control)%count_u=1
              ctrl(n_control)%count_e=1
              ALLOCATE(ctrl(n_control)%K_A(ctrl(n_control)%ord_ctrl,ctrl(n_control)%ord_ctrl))
              ALLOCATE(ctrl(n_control)%K_B(ctrl(n_control)%ord_ctrl,1)) 
              ALLOCATE(ctrl(n_control)%K_C(2,ctrl(n_control)%ord_ctrl))
              ALLOCATE(ctrl(n_control)%K_D(2,1))
              ctrl(n_control)%K_A=0.d0
              ctrl(n_control)%K_B=0.d0
              ctrl(n_control)%K_C=0.d0
              ctrl(n_control)%K_D=0.d0
              do i = 1, ctrl(n_control)%ord_ctrl
                READ(80,*) ctrl(n_control)%K_A(i,:)
                WRITE(*,*)'K_A',ctrl(n_control)%K_A(i,:)
              end do
              READ(80,*) ctrl(n_control)%K_B(1:ctrl(n_control)%ord_ctrl,1)
              READ(80,*) ctrl(n_control)%K_C(:,1:ctrl(n_control)%ord_ctrl)
              READ(80,*) ctrl(n_control)%K_D(:,1)
              WRITE(*,*)'K_B',ctrl(n_control)%K_B(:,:)
              WRITE(*,*)'K_C',ctrl(n_control)%K_C(:,:)
              WRITE(*,*)'K_D',ctrl(n_control)%K_D(:,:)
              ! allocate history arrays: various columns = various actuators
              ALLOCATE(ctrl(n_control)%y_i_2D(1000,ctrl(n_control)%num_act))
              ALLOCATE(ctrl(n_control)%u_i_2D(1000,ctrl(n_control)%num_act))
              ALLOCATE(ctrl(n_control)%e_i_2D(1000,ctrl(n_control)%num_act))
              ALLOCATE(ctrl(n_control)%t_i(1000))
              ctrl(n_control)%e_i_2D=0.d0
              ctrl(n_control)%y_i_2D=0.d0
              ctrl(n_control)%u_i_2D=0.d0
              ctrl(n_control)%t_i=0.d0
              ! allocate internal state arrays: various columns = internal state for various actuators
              ALLOCATE(ctrl(n_control)%tmp_xn(ctrl(n_control)%ord_ctrl,ctrl(n_control)%num_act))
              ALLOCATE(ctrl(n_control)%tmp_xn1(ctrl(n_control)%ord_ctrl,ctrl(n_control)%num_act))
              ctrl(n_control)%tmp_xn=0.d0
              ctrl(n_control)%tmp_xn1=0.d0
              ! read controller logs for internal states (if restart file exists)
              INQUIRE(file=logs_in,exist=readlogs)
              IF (readlogs) then
                OPEN(unit=81,file=logs_in,form='formatted',status='old')
                READ(81,*) dPdx
                READ(81,*) ctrl(n_control)%u_i_2D(1,:)
                READ(81,*) ctrl(n_control)%y_i_2D(1,:)
                WRITE(*,*) 'initial u_i', ctrl(n_control)%u_i_2D(1,:)
                WRITE(*,*) 'initial y_i', ctrl(n_control)%y_i_2D(1,:)
                Do i=1, ctrl(n_control)%num_act
                  READ(81,*) ctrl(n_control)%tmp_xn(1:ctrl(n_control)%ord_ctrl,i)
                  WRITE(*,*)'tmp_Xn (:,i)',ctrl(n_control)%tmp_xn(:,i)
                end do
              end if
            ELSE
              WRITE(*,*) '=> Unknown ctrl_type: ', TRIM(ctrl(n_control)%ctrl_type), '. Skipping.'
              n_control = n_control - 1
            End if

            ! allocate output values array
            IF (n_control .ge. 1) THEN
              ALLOCATE( ctrl(n_control)%values(1:nsteps) )
              ctrl(n_control)%values = 0.d0
            END IF
            WRITE(*,*) 'read control input done...'
          END if
        END DO
        WRITE(*,*) '=> Import all controls. There are',n_control,'control variables.'
      
      end if
      ! broadcast ctrl_type string and all type-specific parameters to all MPI ranks
      Call Mpi_bcast (  n_control,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
      DO i=1,n_control
        Call Mpi_bcast (  ctrl(i)%ctrl_type,LEN(ctrl(i)%ctrl_type),MPI_character,0,MPI_COMM_WORLD,ierr )
        Select case(trim(ctrl(i)%ctrl_type))
          case ('spanwise_const_gauss_x')
            Call Mpi_bcast (  ctrl(i)%num_act,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
            IF (myid /=0) then
              ALLOCATE(ctrl(i)%actuators(1))
              ALLOCATE(ctrl(i)%u_i_2D(1000,ctrl(i)%num_act))
              ALLOCATE(ctrl(i)%y_i_2D(1000,ctrl(i)%num_act))
            end if
            Call Mpi_bcast (  dPdx,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
            Call Mpi_bcast (  ctrl(i)%u_i_2D,1000*ctrl(i)%num_act,MPI_real8,0,MPI_COMM_WORLD,ierr )
            Call Mpi_bcast (  ctrl(i)%y_i_2D,1000*ctrl(i)%num_act,MPI_real8,0,MPI_COMM_WORLD,ierr )
            Call Mpi_bcast (  ctrl(i)%actuators(1)%nx_local,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
            if ( myid /=0 ) then
              ALLOCATE(ctrl(i)%zero_func_local(ctrl(i)%actuators(1)%nx_local))
              ALLOCATE(ctrl(i)%dfdx_local(ctrl(i)%actuators(1)%nx_local))
            end if 
            Call Mpi_bcast (  ctrl(i)%zero_func_local,ctrl(i)%actuators(1)%nx_local,MPI_real8,0,MPI_COMM_WORLD,ierr )
            Call Mpi_bcast (  ctrl(i)%dfdx_local,ctrl(i)%actuators(1)%nx_local,MPI_real8,0,MPI_COMM_WORLD,ierr )
            do j = 1, 1
              Call Mpi_bcast (  ctrl(i)%actuators(j)%amp,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
              Call Mpi_bcast (  ctrl(i)%actuators(j)%sigma,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
              Call Mpi_bcast (  ctrl(i)%actuators(j)%error_actuator,1,MPI_real8,0,MPI_COMM_WORLD,ierr )
            end do
        END SELECT
      end do
      !-------------------------Done--------------------------------!
      Call Mpi_barrier(MPI_COMM_WORLD,ierr)
    END select

  END SUBROUTINE controls_setup_control

! ******************************************************************************
  SUBROUTINE controls_destroy_control
    IMPLICIT NONE
    INTEGER :: i
    DO i = 1, n_control
      DEALLOCATE( ctrl(i)%values )
    END DO
  END SUBROUTINE controls_destroy_control

! ******************************************************************************
!  Top-level loop: dispatch sensing and controller based on ctrl_type
! ******************************************************************************
  Subroutine controls_sensing_and_compute_control
    IMPLICIT NONE
    INTEGER :: i_ctrl
    do i_ctrl = 1, n_control
      Call controls_sensing(i_ctrl)
      Call controls_controller(i_ctrl)
    end do
  end Subroutine controls_sensing_and_compute_control

! ******************************************************************************
!  Sensing: compute wall-normal force integrated over each actuator patch
! ******************************************************************************
  Subroutine controls_sensing(i_ctrl)

    INTEGER, INTENT(in) :: i_ctrl
  
    Integer  (Int32) :: i
  
    Real(Kind(0.d0)), Allocatable :: P_avg_2D(:)
    Real(Kind(0.d0)), Allocatable :: total_fb_n(:)
    Real(Kind(0.d0)), Allocatable :: fbn_2D(:,:)
  
    Select case (trim(ctrl(i_ctrl)%ctrl_type))
      case ('spanwise_const_gauss_x')
      !WRITE(*,*) 'myid',myid, 'sensing...'
      
      ! ----------------------------------------------------------
      ! Compute pressure first 
      ! ----------------------------------------------------------
      !Call compute_pressure
  
      call compute_interior_pressure(p_interior,rhs_p,fb, normals)
      ! ----------------------------------------------------------
      ! Root rank computes actuator sensing
      ! ----------------------------------------------------------
      If ( myid == 0 ) Then
  
         Allocate(P_avg_2D(ctrl(i_ctrl)%num_act))
         Allocate(fbn_2D(nxb,nzb))
  
         P_avg_2D = 0.d0
         total_fb_n=p_interior*sb
         !write(*,*) 'total_fb(nb)=',total_fb_n(nb)
  
         ! reshape global vector into 2D field
         fbn_2D = Reshape(total_fb_n, (/ nxb, nzb /))
  
         ! integrate force over each actuator patch
         Do i = 1, ctrl(i_ctrl)%num_act
  
            P_avg_2D(i) = Sum( &
                 fbn_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                          i   *ctrl(i_ctrl)%actuators(1)%nx_local , : ) )
  
         End Do
         !WRITE(*,*) 'sensing: count_e',ctrl(i_ctrl)%count_e
         IF (ctrl(i_ctrl)%count_e > SIZE(ctrl(i_ctrl)%e_i_2D,1)) THEN
          WRITE(*,*) 'ERROR: count_e overflow'
          WRITE(*,*) 'count_e=', ctrl(i_ctrl)%count_e
          WRITE(*,*) 'size=', SIZE(ctrl(i_ctrl)%e_i_2D,1)
          STOP
         END IF
         ctrl(i_ctrl)%e_i_2D(ctrl(i_ctrl)%count_e,:) = P_avg_2D
  
         ctrl(i_ctrl)%count_e = ctrl(i_ctrl)%count_e + 1
         !WRITE(*,*) 'sensing: count_e',ctrl(i_ctrl)%count_e
         Deallocate(P_avg_2D,fbn_2D,total_fb_n)
      End If
    end select
  
  End Subroutine controls_sensing
! ******************************************************************************
!  Controller: advance state-space equations and compute actuator output
!  State-space form:  x_{n+1} = K_A * x_n + K_B * e_n
!                     u_n     = K_C * x_n + K_D * e_n
! ******************************************************************************
SUBROUTINE controls_controller(i_ctrl)
  IMPLICIT NONE
  REAL(KIND(0.D0)):: tmp_un(2,1),tmp_e_2D(1,1)
  REAL(KIND(0.D0)), ALLOCATABLE :: tmp_xn(:,:),tmp_xn1(:,:)
  INTEGER, INTENT(in) :: i_ctrl
  INTEGER:: i

  IF (TRIM(ctrl(i_ctrl)%ctrl_type) .eq. 'spanwise_const_gauss_x') THEN
    ! one SISO state-space controller per actuator patch
    if (myid .eq. 0) THEN
      !WRITE(*,*) 'myid',myid, 'compute controller output...'
      ALLOCATE(tmp_xn1(ctrl(i_ctrl)%ord_ctrl,1))
      ALLOCATE(tmp_xn(ctrl(i_ctrl)%ord_ctrl,1))
      do i = 1, ctrl(i_ctrl)%num_act
        tmp_xn(:,1)=ctrl(i_ctrl)%tmp_xn(1:ctrl(i_ctrl)%ord_ctrl,i)
        tmp_e_2D=ctrl(i_ctrl)%e_i_2D(ctrl(i_ctrl)%count_e-1,i)
        tmp_xn1=matmul(ctrl(i_ctrl)%K_A,tmp_xn)
        tmp_xn1=tmp_xn1+ctrl(i_ctrl)%K_B*tmp_e_2D(1,1)
        tmp_un=matmul(ctrl(i_ctrl)%K_C,tmp_xn)
        ctrl(i_ctrl)%u_i_2D(ctrl(i_ctrl)%count_u,i)=tmp_un(1,1)+ctrl(i_ctrl)%K_D(1,1)*tmp_e_2D(1,1)
        ctrl(i_ctrl)%y_i_2D(ctrl(i_ctrl)%count_u,i)=tmp_un(2,1)+ctrl(i_ctrl)%K_D(2,1)*tmp_e_2D(1,1)
        ctrl(i_ctrl)%tmp_xn(:,i)=tmp_xn1(:,1)
      end do
      ctrl(i_ctrl)%count_u=ctrl(i_ctrl)%count_u+1
      DEALLOCATE(tmp_xn1,tmp_xn)
    END IF
    Call Mpi_bcast (  ctrl(i_ctrl)%count_u,1,MPI_integer,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  ctrl(i_ctrl)%u_i_2D,1000*ctrl(i_ctrl)%num_act,MPI_real8,0,MPI_COMM_WORLD,ierr )
    Call Mpi_bcast (  ctrl(i_ctrl)%y_i_2D,1000*ctrl(i_ctrl)%num_act,MPI_real8,0,MPI_COMM_WORLD,ierr )
  END IF

END SUBROUTINE controls_controller
! ******************************************************************************
!  Actuator: prescribe wall-normal IB displacement
!
!  yb(x,z) = zero_func_local(x) * u_i_2D(t,patch) * amp
!
!  Construct full 2D IB wall position field, reshape to 1D,
!  then distribute local portion back to each MPI rank.
! ******************************************************************************
SUBROUTINE controls_actuating(i_ctrl)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: i_ctrl

  Integer(Int32) :: i,k
  REAL(KIND(0.D0)) :: amp

  REAL(KIND(0.D0)), Allocatable :: yb_2D(:,:), ub_2D(:,:), slope_2D(:,:), slope_1D(:)

  IF ( TRIM(ctrl(i_ctrl)%ctrl_type) .eq. 'spanwise_const_gauss_x' ) THEN

    amp = ctrl(i_ctrl)%actuators(1)%amp
    !WRITE(*,*) 'myid',myid,'controls_actuating'
    ! WRITE(*,*) 'myid',myid,'(nxb,nzb)=',nxb,nzb
    ! WRITE(*,*) 'myid',myid,'nb=',nb

    ! ----------------------------------------------------------
    ! Root constructs full IB wall position field
    ! ----------------------------------------------------------
    Allocate(yb_2D(nxb,nzb))
    Allocate(ub_2D(nxb,nzb))
    Allocate(slope_2D(nxb,nzb))
    Allocate(slope_1D(nb))

    yb_2D = 0.d0
    ub_2D = 0.d0
    slope_2D=0.d0
    slope_1D=0.d0

    ! build actuator forcing in x-z plane
    Do k = 1, nzb

      Do i = 1, ctrl(i_ctrl)%num_act

        If ( istep .le. 1 .or. ctrl(i_ctrl)%count_u .eq. 1) Then

            If ( istep .le. 1) Then

              ! displacement
              yb_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                        i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                        ctrl(i_ctrl)%zero_func_local * &
                        ctrl(i_ctrl)%y_i_2D(1,i) * amp

              ! velocity
              ub_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                        i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                        ctrl(i_ctrl)%zero_func_local * &
                        ctrl(i_ctrl)%u_i_2D(1,i) * amp

              ! slope / normal information
              slope_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                          i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                          ctrl(i_ctrl)%dfdx_local * &
                          ctrl(i_ctrl)%y_i_2D(1,i) * amp

            Else

              yb_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                        i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                        ctrl(i_ctrl)%zero_func_local * &
                        ctrl(i_ctrl)%y_i_2D(1000,i) * amp

              ub_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                        i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                        ctrl(i_ctrl)%zero_func_local * &
                        ctrl(i_ctrl)%u_i_2D(1000,i) * amp

              slope_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                          i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                          ctrl(i_ctrl)%dfdx_local * &
                          ctrl(i_ctrl)%y_i_2D(1000,i) * amp

            End If

        Else

            yb_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                    i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                    ctrl(i_ctrl)%zero_func_local * &
                    ctrl(i_ctrl)%y_i_2D(ctrl(i_ctrl)%count_u-1,i) * amp

            ub_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                    i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                    ctrl(i_ctrl)%zero_func_local * &
                    ctrl(i_ctrl)%u_i_2D(ctrl(i_ctrl)%count_u-1,i) * amp

            slope_2D( (i-1)*ctrl(i_ctrl)%actuators(1)%nx_local + 1 : &
                        i   *ctrl(i_ctrl)%actuators(1)%nx_local , k ) = &
                        ctrl(i_ctrl)%dfdx_local * &
                        ctrl(i_ctrl)%y_i_2D(ctrl(i_ctrl)%count_u-1,i) * amp

        End If

      End Do

    End Do
    ! reshape 2D field into global IB vector
    !WRITE(*,*) 'myid',myid,'finish 2D matrix'
    yb(:) = Reshape(yb_2D, (/ nb /))
    ub(nb+1:2*nb) = Reshape(ub_2D, (/ nb /))
    !assigne normal and tangential vector based on the slope vector
    slope_1D=reshape(slope_2D, (/ nb /) )
    normals(1:nb)=slope_1D
    normals(nb+1:2*nb)=-1
    normals(2*nb+1:3*nb)=0.d0
    tangents_1(1:nb)=1
    tangents_1(nb+1:2*nb)=slope_1D
    tangents_1(2*nb+1:3*nb)=0.d0
    tangents_2=0.d0
    tangents_2(2*nb+1:3*nb)=-1
    !WRITE(*,*) 'myid',myid,'finish controls_actuating'


    ! ----------------------------------------------------------
    ! Clean up
    ! ----------------------------------------------------------
    Deallocate(yb_2D, ub_2D, slope_1D, slope_2D)


  END IF

END SUBROUTINE controls_actuating

END MODULE controls


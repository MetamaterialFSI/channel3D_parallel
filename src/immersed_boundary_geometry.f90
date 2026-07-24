Module immersed_boundary_geometry

  Use iso_fortran_env, Only : error_unit, Int32, Int64
  Use global

  ! prevent implicit typing
  Implicit None

  ! declarations
Contains

  Subroutine compute_nb
    Integer(Int32) :: nxb1, nxb2
    Real   (Int64) :: r, r1, r2

    Select Case (trim(body_type))
      Case ('none') ! No IB
        nbodies = 0
        nxb = 0
        nzb = 0
        nb = 0
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case ('center_wall') ! Static planar IB wall centered in y
        nbodies = 1
        nb = nxb * nzb
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case ('double_cylinders_z') ! Double concentric cylinders with axis parallel to z
        nbodies = 2
        r1 = body_param_1
        r2 = body_param_2
        dzb = real(Lzp / nzb, 8)
        dxb = dzb * dx / dz
        nxb1 = int(2 * 3.14159 * r1 / dxb)
        nxb2 = int(2 * 3.14159 * r2 / dxb)
        nxb = nxb1 + nxb2
        nb = nxb * nzb

      Case ('single_cylinder_z')
        nbodies = 1
        r = 1d0
        dzb = real(Lzp / nzb, 8)
        dxb = dzb * dx / dz
        nxb = int(2 * 3.14159 * r / dxb)
        nb = nxb * nzb

      Case ('center_standing_wave_x') ! Center wall undergoing standing wave motion
        nbodies = 1
        nb = nxb * nzb
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)
        
      Case ('standing_wave_x') ! Top and bottom wall undergoing standing wave motion
        nbodies = 2
        nb = 2 * nxb * nzb
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case ('traveling_wave_x') ! Top and bottom wall undergoing traveling wave motion in the x-direction
        nbodies = 2
        nb = 2 * nxb * nzb
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case ('traveling_wave_z') ! Top and bottom wall undergoing traveling wave motion in the z-direction
        nbodies = 2
        nb = 2 * nxb * nzb
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case DEFAULT
        If (myid == 0) Then
            Write(*,*) 'Error: No corresponding body type found for: ', trim(body_type)
        End If
        Error Stop 'Invalid body type'
    End Select


  End Subroutine compute_nb

  Subroutine setup_IB_geometry
    Integer(Int32) :: i, j, k, l, nxb1, nxb2
    Real   (Int64) :: a1, a2, r, r1, r2, xc, yc, theta, dsb, dsb1, dsb2, phi, amp, slope_bottom, slope_top, phase, damp_dt, tau

    Select Case (trim(body_type))
      Case ('none') ! No IB
        moving_body = .False.
        nb_start = 0
        nb_end = -1

      Case ('center_wall') ! Static planar wall centered at y = 1
        If ( grid_type /= 0 ) Stop 'Error: body type is incompatible with grid type'
        moving_body = .False.
        moving_z_flag = .False.

        ub = 0d0
        ! Reference points are the center of the domain
        y_ref_index = ny_global / 2 ! automatically rounds down

        ! Scalar arrays. Arrange such that the points treated by one partition are contiguous
        ! (i.e., fall between an nb_start and nb_end)
        nb_start = nb + 1  ! Initialize to an invalid value (beyond the max index)
        nb_end = 1         ! Initialize to the lowest possible index
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + (j-1) * nxb
            xb(k) = (real(i,8) - 0.5d0) * dxb
            yb(k) = 0.5d0 * Ly_channel
            zb(k) = (real(j,8) - 0.5d0) * dzb
          End Do
          If (zb((j-1) * nxb + 1) >= z(1) .and. nb_start > (j-1) * nxb + 1) then
            nb_start = (j-1) * nxb + 1
          End If
          If (zb(j * nxb) < z(nz-1) .and. nb_end < j * nxb) then
            nb_end = j * nxb
          End If
        End Do
        
        sb = dxb * dzb

        ! Vector arrays
        Do k=1,nb
          tangents_1(k) = 1d0
          tangents_2(2*nb + k) = -1d0
          normals(nb + k) = -1d0
        End Do

      Case ('double_cylinders_z') ! Double rotating cylinders
        If ( grid_type /= 0 ) Stop 'Error: body type is incompatible with grid type'
        moving_body = .False. 
        moving_z_flag = .False.
        nb_start = nb + 1  ! Initialize to an invalid value (beyond the max index)
        nb_end = 1         ! Initialize to the lowest possible index
        y_ref_index = 1 ! The grid has to be uniform for this case, so it doesn't matter what y_ref_index is

        r1 = body_param_1
        r2 = body_param_2
        xc = 0.5d0 * Lxp
        yc = 0.5d0 * Ly_channel
        nxb1 = Int(2 * 3.14159 * r1 / dxb)
        nxb2 = Int(2 * 3.14159 * r2 / dxb)
        nxb = nxb1 + nxb2
        dsb1 = 2 * 3.14159 * r1 / nxb1
        dsb2 = 2 * 3.14159 * r2 / nxb2
        ub = 0d0
        Do j=1,nzb
          ! Inner cylinder
          Do i = 1, nxb1
            theta = Real(i - 1, 8) * dsb1 / r1
            xb(i + (j - 1) * nxb) = r1 * cos(theta) + xc
            yb(i + (j - 1) * nxb) = r1 * sin(theta) + yc
            zb(i + (j - 1) * nxb) = (Real(j,8) - 0.5d0) * dzb
            ub(i + (j - 1) * nxb) = -sin(theta) * body_param_3
            ub(nb + i + (j - 1) * nxb) = cos(theta) * body_param_3
            sb(i + (j - 1) * nxb) = dsb1 * dzb
            normals(i + (j - 1) * nxb) = -cos(theta)
            normals(nb + i + (j - 1) * nxb) = -sin(theta)
            tangents_1(i + (j - 1) * nxb) = sin(theta)
            tangents_1(nb + i + (j - 1) * nxb) = -cos(theta)
          End Do

          ! Outer cylinder
          Do i = 1, nxb2
            theta = Real(i - 1, 8) * dsb2 / r2
            xb(i + (j - 1) * nxb + nxb1) = r2 * cos(theta) + xc
            yb(i + (j - 1) * nxb + nxb1) = r2 * sin(theta) + yc
            zb(i + (j - 1) * nxb + nxb1) = (Real(j,8) - 0.5d0) * dzb
            sb(i + (j - 1) * nxb + nxb1) = dsb2 * dzb
            normals(i + (j - 1) * nxb + nxb1) = cos(theta)
            normals(nb + i + (j - 1) * nxb + nxb1) = sin(theta)
            tangents_1(i + (j - 1) * nxb + nxb1) = -sin(theta)
            tangents_1(nb + i + (j - 1) * nxb + nxb1) = cos(theta)
          End Do

          Do i = 1, nb
            tangents_2(2 * nb + i) = 1d0
          End Do

          If (zb((j-1) * nxb + 1) >= z(1) .and. nb_start > (j-1) * nxb + 1) then
            nb_start = (j-1) * nxb + 1
          End If
          If (zb(j * nxb) < z(nz-1) .and. nb_end < j * nxb) then
            nb_end = j * nxb
          End If
        End Do

      Case ('single_cylinder_z')
        If ( grid_type /= 0 ) Stop 'Error: body type is incompatible with grid type'
        moving_body = .True. 
        moving_z_flag = .False.
        nb_start = nb + 1  ! Initialize to an invalid value (beyond the max index)
        nb_end = 1         ! Initialize to the lowest possible index
        y_ref_index = 1 ! The grid has to be uniform for this case, so it doesn't matter what y_ref_index is

        if (t < body_ramp_up_time .and. body_ramp_up_time > 0d0) then
          tau = t / body_ramp_up_time
          amp = body_param_2 * 0.5d0 * (1d0 - cos(pi * tau))
          ! integrate ramped velocity: integral of 0.5*(1-cos(pi*t/t))dt
          yc = 0.5d0 * Ly_channel + body_param_2 * 0.5d0 * &
               (t - (body_ramp_up_time / pi) * sin(pi * tau))
        else
          amp = body_param_2
          ! ramp displacement at t=t, then constant velocity after
          yc = 0.5d0 * Ly_channel + body_param_2 * 0.5d0 * body_ramp_up_time &
               + body_param_2 * (t - body_ramp_up_time)
        end if

        ! arc length between circle points equals dxb
        r = body_param_1
        xc = 0.5d0 * Lxp
        dsb = 2 * 3.14159 * r / nxb
        ub = 0d0

        Do j = 1, nzb
          do i = 1, nxb
            theta = real(i - 1, 8) * dsb / r
            xb(i + (j - 1) * nxb) = r * cos(theta) + xc
            yb(i + (j - 1) * nxb) = r * sin(theta) + yc
            zb(i + (j - 1) * nxb) = (real(j,8) - 0.5d0) * dzb
            ub(nb + i + (j - 1) * nxb) = amp
            sb(i + (j - 1) * nxb) = dsb * dzb
            normals(i + (j - 1) * nxb) = cos(theta)
            normals(nb + i + (j - 1) * nxb) = sin(theta)
            tangents_1(i + (j - 1) * nxb) = -sin(theta)
            tangents_1(nb + i + (j - 1) * nxb) = cos(theta)
          end do

          If (zb((j-1) * nxb + 1) >= z(1) .and. nb_start > (j-1) * nxb + 1) then
            nb_start = (j-1) * nxb + 1
          End If

          If (zb(j * nxb) < z(nz-1) .and. nb_end < j * nxb) then
            nb_end = j * nxb
          End If

        end do
        do k = 1, nb
          tangents_2(2*nb + k) = 1d0
        end do

        do k = 1, nb
          if (xb(k) >= Lxp) then
            If ( myid==0 ) Then
              print *, "Error: xb(", k, ") = ", xb(k), " is greater than or equal to Lxp =", Lxp
              stop 1
            End If
          end if
        end do

      Case ('center_standing_wave_x') ! Center wall undergoing standing wave motion in x-direction
        If ( grid_type /= 0 ) Stop 'Error: body type is incompatible with grid type'
        moving_body = .True.
        moving_z_flag = .False.

        ! Scalar arrays. Arrange such that the points treated by one partition are contiguous
        ! (i.e., fall between an nb_start and nb_end)
        nb_start = nb + 1 ! Initialize to an invalid value (beyond the max index)
        nb_end = 1        ! Initialize to the lowest possible index
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + nxb * (j - 1)
            xb(k) = (real(i,8) - 0.5d0) * dxb
            yb(k) = 0.5d0 * Ly_channel + body_param_1 * sin(2d0 * pi * body_param_3 * xb(k) / Lxp) * cos(body_param_2 * t)
            zb(k) = (real(j,8) - 0.5d0) * dzb


            If (zb(k) >= z(1) .and. nb_start > k) then
              nb_start = k
            End If
            If (zb(k) < z(nz-1) .and. nb_end < k + nxb) then
              nb_end = k
            End If
          End Do
        End Do
        y_ref_index = ny_global / 2 ! automatically rounds down

        ! Vector arrays
        ub(1:nb) = 0d0
        ub(2 * nb + 1 : 3 * nb) = 0d0
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + nxb * (j - 1)
            phase = 2d0 * pi * body_param_3 * xb(k) / Lxp
            ub(nb + k)       = -body_param_1 * body_param_2 * sin(phase) * sin(body_param_2 * t)

            slope_bottom = body_param_1 * 2d0 * pi * body_param_3 / Lxp * cos(phase) * cos(body_param_2 * t)

            ! Tangent orthogonal to z
            tangents_1(k)            = 1 / sqrt( 1 + slope_bottom ** 2 )
            tangents_1(nb + k)       = slope_bottom / sqrt( 1 + slope_bottom ** 2 )

            ! Tangent parallel to z
            tangents_2(2 * nb + k      ) = 1d0

            normals(k)            = tangents_1(nb + k)
            normals(nb + k)       = -tangents_1(k)
          End Do
        End Do
        ! For now, make a naive sb calculation that assumes no variation in z
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + nxb * (j - 1)
            If (i .eq. 1) Then
              a1 = sqrt((xb(k) - (xb(nxb + nxb * (j - 1)) - Lxp)) ** 2 + (yb(k) - yb(nxb + nxb * (j - 1))) ** 2)
              a2 = sqrt((xb(k) - (xb(k + 1))                        ) ** 2 + (yb(k) - yb(k + 1)                  ) ** 2)
            Else If (i .eq. nxb) Then
              a1 = sqrt((xb(k) - (xb(k - 1))                      ) ** 2   + (yb(k) - yb(k - 1)                  ) ** 2)
              a2 = sqrt((xb(k) - (xb(1 + nxb * (j - 1)) + Lxp)) ** 2   + (yb(k) - yb(1 + nxb * (j - 1))) ** 2)
            Else
              a1 = sqrt((xb(k) - (xb(k - 1))) ** 2                     + (yb(k) - yb(k - 1)) ** 2)
              a2 = sqrt((xb(k) - (xb(k + 1))) ** 2                     + (yb(k) - yb(k + 1)) ** 2)
            End If
            sb(k)       = dzb * (0.5d0 * a1 + 0.5d0 * a2)
          End Do
        End Do

      Case ('standing_wave_x') ! Top and bottom wall undergoing standing wave motion in x-direction
        If ( body_param_1 > min_buffer_width ) Stop 'Error: IB amplitude is bigger than the minimum buffer width'
        moving_body = .True.
        moving_z_flag = .False.

        ! Scalar arrays. Arrange such that the points treated by one partition are contiguous
        ! (i.e., fall between an nb_start and nb_end)
        nb_start = nb + 1 ! Initialize to an invalid value (beyond the max index)
        nb_end = 1        ! Initialize to the lowest possible index
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + 2 * nxb * (j - 1)
            xb(k)       = (real(i,8) - 0.5d0) * dxb
            xb(k + nxb) = (real(i,8) - 0.5d0) * dxb
            yb(k)       =              body_param_1 * sin(2d0 * pi * body_param_3 * xb(k) / Lxp) * cos(body_param_2 * t)
            yb(k + nxb) = Ly_channel + body_param_1 * sin(2d0 * pi * body_param_3 * xb(k) / Lxp) * cos(body_param_2 * t)
            zb(k)       = (real(j,8) - 0.5d0) * dzb
            zb(k + nxb) = (real(j,8) - 0.5d0) * dzb

            y_ref_index(k) = 1
            y_ref_index(k + nxb) = ny_global

            If (yb(k) < y(1 + suppy)) Then
              If ( myid==0 ) Then
                write(*,'(A,I3,A,F13.6,A,F13.6)') "yb(", k, ") = ", yb(k), &
                  " is smaller than y(1 + suppy) = ", y(1 + suppy)
                write(*,*) "Error: body points support exceeds grid dimensions"
              End If
              Stop
            End If
            If (yb(k + nxb) > y(ny_global - suppy)) Then
              If ( myid==0 ) Then
                write(*,'(A,I3,A,F13.6,A,F13.6)') "yb(", k + nxb, ") = ", yb(k + nxb), &
                  " is greater than y(ny_global - suppy) = ", y(ny_global - suppy)
                write(*,*) "Error: body points support exceeds grid dimensions"
              End If
              Stop
            End If

            If (zb(k) >= z(1) .and. nb_start > k) then
              nb_start = k
            End If
            If (zb(k + nxb) < z(nz-1) .and. nb_end < k + nxb) then
              nb_end = k + nxb
            End If
          End Do
        End Do
        ! Vector arrays
        ub(1:nb) = 0d0
        ub(2 * nb + 1 : 3 * nb) = 0d0
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + 2 * nxb * (j - 1)
            phase = 2d0 * pi * body_param_3 * xb(k) / Lxp
            ub(nb + k)       = -body_param_1 * body_param_2 * sin(phase) * sin(body_param_2 * t)
            ub(nb + k + nxb) = -body_param_1 * body_param_2 * sin(phase) * sin(body_param_2 * t)

            slope_bottom = body_param_1 * 2d0 * pi * body_param_3 / Lxp * cos(phase) * cos(body_param_2 * t)
            slope_top    = body_param_1 * 2d0 * pi * body_param_3 / Lxp * cos(phase) * cos(body_param_2 * t)

            ! Tangent orthogonal to z
            tangents_1(k)            = 1 / sqrt( 1 + slope_bottom ** 2 )
            tangents_1(k + nxb)      = 1 / sqrt( 1 + slope_top    ** 2 )
            tangents_1(nb + k)       = slope_bottom / sqrt( 1 + slope_bottom ** 2 )
            tangents_1(nb + k + nxb) = slope_top    / sqrt( 1 + slope_top    ** 2 )

            ! Tangent parallel to z
            tangents_2(2 * nb + k      ) = 1d0
            tangents_2(2 * nb + k + nxb) = 1d0

            normals(k)            = tangents_1(nb + k)
            normals(k + nxb)      = -tangents_1(nb + k + nxb)
            normals(nb + k)       = -tangents_1(k)
            normals(nb + k + nxb) = tangents_1(k + nxb)
          End Do
        End Do
        ! For now, make a naive sb calculation that assumes no variation in z
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + 2 * nxb * (j - 1)
            If (i .eq. 1) Then
              a1 = sqrt((xb(k) - (xb(nxb + 2 * nxb * (j - 1)) - Lxp)) ** 2 + (yb(k) - yb(nxb + 2 * nxb * (j - 1))) ** 2)
              a2 = sqrt((xb(k) - (xb(k + 1))                        ) ** 2 + (yb(k) - yb(k + 1)                  ) ** 2)
            Else If (i .eq. nxb) Then
              a1 = sqrt((xb(k) - (xb(k - 1))                      ) ** 2   + (yb(k) - yb(k - 1)                  ) ** 2)
              a2 = sqrt((xb(k) - (xb(1 + 2 * nxb * (j - 1)) + Lxp)) ** 2   + (yb(k) - yb(1 + 2 * nxb * (j - 1))) ** 2)
            Else
              a1 = sqrt((xb(k) - (xb(k - 1))) ** 2                     + (yb(k) - yb(k - 1)) ** 2)
              a2 = sqrt((xb(k) - (xb(k + 1))) ** 2                     + (yb(k) - yb(k + 1)) ** 2)
            End If
            sb(k)       = dzb * (0.5d0 * a1 + 0.5d0 * a2)
            sb(k + nxb) = dzb * (0.5d0 * a1 + 0.5d0 * a2) ! Assumes that the top and bottom wall undergo the same motion!
          End Do
        End Do

      Case ('traveling_wave_x') ! Top and bottom wall undergoing traveling wave motion
        If ( body_param_1 / (body_param_2 * body_param_3) > min_buffer_width ) Stop 'Error: IB amplitude is bigger than the minimum buffer width'
        moving_body = .True.
        moving_z_flag = .False.

        if (t < body_ramp_up_time .and. body_ramp_up_time > 0d0) then
          tau = t / body_ramp_up_time
          amp = body_param_1 * 0.5d0 * (1d0 - cos(pi * tau))
          damp_dt = body_param_1 * 0.5d0 * (pi / body_ramp_up_time) * sin(pi * tau)
        else
          amp = body_param_1
          damp_dt = 0d0
        end if
        ! Phase difference between top and bottom wave motion
        phi = pi

        ! Scalar arrays. Arrange such that the points treated by one partition are contiguous
        ! (i.e., fall between an nb_start and nb_end)
        nb_start = nb + 1 ! Initialize to an invalid value (beyond the max index)
        nb_end = 1        ! Initialize to the lowest possible index
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + 2 * nxb * (j - 1)
            xb(k)       = (real(i,8) - 0.5d0) * dxb
            xb(k + nxb) = (real(i,8) - 0.5d0) * dxb
            phase = body_param_3 * (xb(k) - body_param_2 * t)
            yb(k)       =              amp / (body_param_2 * body_param_3) * sin(phase)
            yb(k + nxb) = Ly_channel + amp / (body_param_2 * body_param_3) * sin(phase + phi)
            zb(k)       = (real(j,8) - 0.5d0) * dzb
            zb(k + nxb) = (real(j,8) - 0.5d0) * dzb

            y_ref_index(k) = 1
            y_ref_index(k + nxb) = ny_global

            If (yb(k) < y(1 + suppy)) Then
              If ( myid==0 ) Then
                write(*,'(A,I3,A,F13.6,A,F13.6)') "yb(", k, ") = ", yb(k), &
                  " is smaller than y(1 + suppy) = ", y(1 + suppy)
                write(*,*) "Error: body points support exceeds grid dimensions"
              End If
              Stop
            End If
            If (yb(k + nxb) > y(ny_global - suppy)) Then
              If ( myid==0 ) Then
                write(*,'(A,I3,A,F13.6,A,F13.6)') "yb(", k + nxb, ") = ", yb(k + nxb), &
                  " is greater than y(ny_global - suppy) = ", y(ny_global - suppy)
                write(*,*) "Error: body points support exceeds grid dimensions"
              End If
              Stop
            End If

          End Do
          
          If (zb((j-1) * 2 * nxb + 1) >= z(1) .and. nb_start > (j-1) * 2 * nxb + 1) then
            nb_start = (j-1) * 2 * nxb + 1
          End If
          If (zb(j * 2 * nxb) < z(nz-1) .and. nb_end < j * 2 * nxb) then
            nb_end = j * 2 * nxb
          End If

        End Do
        ! Vector arrays
        ub(1:nb) = 0d0
        ub(2 * nb + 1 : 3 * nb) = 0d0
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + 2 * nxb * (j - 1)
            phase = body_param_3 * (xb(k) - body_param_2 * t)

            ub(nb + k) = - amp * cos(phase) + (damp_dt / (body_param_2 * body_param_3)) * sin(phase)
            ub(nb + k + nxb) = - amp * cos(phase + phi) + (damp_dt / (body_param_2 * body_param_3)) * sin(phase + phi)

            slope_bottom = amp / body_param_2 * cos(phase)
            slope_top    = amp / body_param_2 * cos(phase + phi)

            ! Tangent orthogonal to z
            tangents_1(k)            = 1 / sqrt( 1 + slope_bottom ** 2 )
            tangents_1(k + nxb)      = 1 / sqrt( 1 + slope_top    ** 2 )
            tangents_1(nb + k)       = slope_bottom / sqrt( 1 + slope_bottom ** 2 )
            tangents_1(nb + k + nxb) = slope_top    / sqrt( 1 + slope_top    ** 2 )

            ! Tangent parallel to z
            tangents_2(2 * nb + k      ) = 1d0
            tangents_2(2 * nb + k + nxb) = 1d0

            normals(k)            = tangents_1(nb + k)
            normals(k + nxb)      = -tangents_1(nb + k + nxb)
            normals(nb + k)       = -tangents_1(k)
            normals(nb + k + nxb) = tangents_1(k + nxb)
          End Do
        End Do
        ! For now, make a naive sb calculation that assumes no variation in z
        Do j = 1, nzb
          Do i = 1, nxb
            ! Bottom wall
            k = i + 2 * nxb * (j - 1)
            If (i .eq. 1) Then
              a1 = sqrt((xb(k) - (xb(nxb + 2 * nxb * (j - 1)) - Lxp)) ** 2 + (yb(k) - yb(nxb + 2 * nxb * (j - 1))) ** 2)
              a2 = sqrt((xb(k) - (xb(k + 1))                        ) ** 2 + (yb(k) - yb(k + 1)                  ) ** 2)
            Else If (i .eq. nxb) Then
              a1 = sqrt((xb(k) - (xb(k - 1))                      ) ** 2   + (yb(k) - yb(k - 1)                  ) ** 2)
              a2 = sqrt((xb(k) - (xb(1 + 2 * nxb * (j - 1)) + Lxp)) ** 2   + (yb(k) - yb(1 + 2 * nxb * (j - 1))) ** 2)
            Else
              a1 = sqrt((xb(k) - (xb(k - 1))) ** 2                     + (yb(k) - yb(k - 1)) ** 2)
              a2 = sqrt((xb(k) - (xb(k + 1))) ** 2                     + (yb(k) - yb(k + 1)) ** 2)
            End If
            sb(k) = dzb * (0.5d0 * a1 + 0.5d0 * a2)
            ! Top wall
            k = i + nxb + 2 * nxb * (j - 1)
            If (i .eq. 1) Then
              a1 = sqrt((xb(k) - (xb(2 * nxb * (j)) - Lxp)) ** 2               + (yb(k) - yb(2 * nxb * (j))) ** 2)
              a2 = sqrt((xb(k) - (xb(k + 1))                        ) ** 2     + (yb(k) - yb(k + 1)                  ) ** 2)
            Else If (i .eq. nxb) Then
              a1 = sqrt((xb(k) - (xb(k - 1))                      ) ** 2       + (yb(k) - yb(k - 1)                  ) ** 2)
              a2 = sqrt((xb(k) - (xb(nxb + 1 + 2 * nxb * (j - 1)) + Lxp)) ** 2 + (yb(k) - yb(nxb + 1 + 2 * nxb * (j - 1))) ** 2)
            Else
              a1 = sqrt((xb(k) - (xb(k - 1))) ** 2                     + (yb(k) - yb(k - 1)) ** 2)
              a2 = sqrt((xb(k) - (xb(k + 1))) ** 2                     + (yb(k) - yb(k + 1)) ** 2)
            End If
            sb(k) = dzb * (0.5d0 * a1 + 0.5d0 * a2)
          End Do
        End Do

      Case ('traveling_wave_z') ! Top and bottom wall undergoing traveling wave motion
        Write(*,*) 'traveling wave in z-direction not yet implemented'
        Stop 

    End Select

  End Subroutine setup_IB_geometry

  Subroutine remove_mean_per_body(f_)
    Implicit None
    Real(Int64), Dimension(nb), Intent(InOut) :: f_
    Real(Int64), Dimension(nbodies) :: favg_
    Real(Int64) :: r1, r2, xc, yc
    Integer(Int32) :: j
    Integer(Int32) :: nxb1, nxb2

    favg_ = 0d0

    Select Case (trim(body_type))

      Case ('center_wall', 'single_cylinder_z', 'center_standing_wave_x')
        f_ = f_ - sum(f_) / size(f_)

      Case ('double_cylinders_z')
        r1 = body_param_1
        r2 = body_param_2
        xc = 0.5d0 * Lxp
        yc = 0.5d0 * Ly_channel
        nxb1 = int(2 * 3.14159 * r1 / dxb)
        nxb2 = int(2 * 3.14159 * r2 / dxb)
        
        Do j = 1, nzb
          ! Inner cylinder
          favg_(1) = favg_(1) + sum(f_((j-1)*nxb + 1 : (j-1)*nxb + nxb1))

          ! Outer cylinder
          favg_(2) = favg_(2) + sum(f_((j-1)*nxb + nxb1 + 1 : j*nxb))
        End Do

        favg_(1) = favg_(1) / (nxb1 * nzb)
        favg_(2) = favg_(2) / (nxb2 * nzb) 

        Do j = 1, nzb
          ! Inner cylinder
          f_((j-1)*nxb + 1 : (j-1)*nxb + nxb1) = &
               f_((j-1)*nxb + 1 : (j-1)*nxb + nxb1) - favg_(1)

          ! Outer cylinder
          f_((j-1)*nxb + nxb1 + 1 : j*nxb) = &
               f_((j-1)*nxb + nxb1 + 1 : j*nxb) - favg_(2)
        End Do

      Case ('standing_wave_x', 'traveling_wave_x')
        Do j = 1, nzb
          ! Bottom wall
          favg_(1) = favg_(1) + sum(f_((j-1)*2*nxb + 1       : (j-1)*2*nxb + nxb  ))

          ! Top wall
          favg_(2) = favg_(2) + sum(f_((j-1)*2*nxb + nxb + 1 : (j-1)*2*nxb + 2*nxb))
        End Do

        favg_(1) = favg_(1) / (nxb * nzb)
        favg_(2) = favg_(2) / (nxb * nzb) 

        Do j = 1, nzb
          ! Bottom wall
          f_((j-1)*2*nxb + 1       : (j-1)*2*nxb + nxb  ) = f_((j-1)*2*nxb + 1       : (j-1)*2*nxb + nxb  ) - favg_(1)

          ! Top wall
          f_((j-1)*2*nxb + nxb + 1 : (j-1)*2*nxb + 2*nxb) = f_((j-1)*2*nxb + nxb + 1 : (j-1)*2*nxb + 2*nxb) - favg_(2)
        End Do


      Case Default
        stop "subroutine remove_mean_per_body not implemented for current case" 

    End Select

  End Subroutine

  Subroutine exact_masked_fields(U_, V_, W_, P_, Hu_exterior_, Hv_exterior_, Hw_exterior_, Hc_exterior_)
    Implicit None
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: U_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: V_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: W_(:, :, :)
    Real   (Int64), Contiguous, INTENT(INOUT)  :: P_(2:, 2:, 2:)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: Hu_exterior_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: Hv_exterior_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: Hw_exterior_(:, :, :)
    Real   (Int64), CONTIGUOUS, INTENT(INOUT)  :: Hc_exterior_(2:, 2:, 2:)

    ! Local variables
    Integer :: i, j, k
    Real(Int64) :: r1, r2, xc, yc, om, r_ext_smooth_start, r_ext_smooth_end, r_int_smooth_1_start, r_int_smooth_1_end
    Real(Int64) :: r_int_smooth_2_start, r_int_smooth_2_end
    Real(Int64) :: xloc, yloc, r, theta, vel_theta_exterior, vel_theta_interior
    Real(Int64) :: p_exterior, p_interior, p_jump

    Select Case (trim(body_type))

      Case ('double_cylinders_z')

        p_jump = 0d0

        r1 = body_param_1
        r2 = body_param_2
        xc = 0.5d0 * Lxp
        yc = 0.5d0 * Ly_channel
        om = body_param_3 / r1 ! angular velocity

        r_ext_smooth_start = r1 + 0.2 * (r2 - r1)
        r_ext_smooth_end = r1 + 0.8 * (r2 - r1)
        r_int_smooth_1_start = 0.3 * (r1)
        r_int_smooth_1_end = 0.8 * (r1)
        r_int_smooth_2_start = 1.2 * (r2)
        r_int_smooth_2_end = 1.3 * (r2)

        ! U velocities (x-faces)
        Do k = 1, nzg
          Do j = 1, nyg
            Do i = 1, nx
              xloc = x(i)  - xc      ! x-face
              yloc = yg(j) - yc      ! y-center
              r = sqrt(xloc ** 2 + yloc ** 2)
              theta = atan2(yloc, xloc)

              vel_theta_exterior = 0
              vel_theta_interior = 0
              If (r < r_ext_smooth_start) Then
                vel_theta_exterior = -om * r
              ElseIf (r >= r_ext_smooth_start .and. r <= r_ext_smooth_end) Then
                vel_theta_exterior = -om * r &
                  * 0.5 * (1 + cos(pi * (r - r_ext_smooth_start) / (r_ext_smooth_end - r_ext_smooth_start)))
              End If
                
              If (r >= r_int_smooth_1_start .and. r < r_int_smooth_1_end) Then
                vel_theta_interior = -om * r1 ** 2 / (r2 ** 2 - r1 ** 2) * (r2 ** 2 / r - r) &
                  * 0.5 * (1 + cos(pi * (r - r_int_smooth_1_end) / (r_int_smooth_1_start - r_int_smooth_1_end)))
              ElseIf (r >= r_int_smooth_1_end .and. r <= r_int_smooth_2_start) Then
                vel_theta_interior = -om * r1 ** 2 / (r2 ** 2 - r1 ** 2) * (r2 ** 2 / r - r)
              ElseIf (r >= r_int_smooth_2_start .and. r <= r_int_smooth_2_end) Then
                vel_theta_interior = -om * r1 ** 2 / (r2 ** 2 - r1 ** 2) * (r2 ** 2 / r - r) &
                  * 0.5 * (1 + cos(pi * (r - r_int_smooth_2_start) / (r_int_smooth_2_end - r_int_smooth_2_start)))
              End If


              U_(i,j,k) = vel_theta_exterior * sin(theta) * Hu_exterior_(i,j,k) + &
                vel_theta_interior * sin(theta) * (1 - Hu_exterior_(i,j,k))
            End Do
          End Do
        End Do

        ! V velocities (y-faces)
        Do k = 1, nzg
          Do j = 1, ny
            Do i = 1, nxg
              xloc = xg(i) - xc      ! x-center
              yloc = y(j)  - yc      ! y-face
              r = sqrt(xloc**2 + yloc**2)
              theta = atan2(yloc, xloc)

              vel_theta_exterior = 0
              vel_theta_interior = 0
              If (r < r_ext_smooth_start) Then
                vel_theta_exterior = -om * r
              ElseIf (r >= r_ext_smooth_start .and. r <= r_ext_smooth_end) Then
                vel_theta_exterior = -om * r &
                  * 0.5 * (1 + cos(pi * (r - r_ext_smooth_start) / (r_ext_smooth_end - r_ext_smooth_start)))
              End If
                
              If (r >= r_int_smooth_1_start .and. r < r_int_smooth_1_end) Then
                vel_theta_interior = -om * r1 ** 2 / (r2 ** 2 - r1 ** 2) * (r2 ** 2 / r - r) &
                  * 0.5 * (1 + cos(pi * (r - r_int_smooth_1_end) / (r_int_smooth_1_start - r_int_smooth_1_end)))
              ElseIf (r >= r_int_smooth_1_end .and. r <= r_int_smooth_2_start) Then
                vel_theta_interior = -om * r1 ** 2 / (r2 ** 2 - r1 ** 2) * (r2 ** 2 / r - r)
              ElseIf (r >= r_int_smooth_2_start .and. r <= r_int_smooth_2_end) Then
                vel_theta_interior = -om * r1 ** 2 / (r2 ** 2 - r1 ** 2) * (r2 ** 2 / r - r) &
                  * 0.5 * (1 + cos(pi * (r - r_int_smooth_2_start) / (r_int_smooth_2_end - r_int_smooth_2_start)))
              End If

              V_(i,j,k) = -vel_theta_exterior * cos(theta) * Hv_exterior_(i,j,k) - &
                vel_theta_interior * cos(theta) * (1 - Hv_exterior_(i,j,k))
            End Do
          End Do
        End Do

        ! W velocities (z-faces, always zero here)
        Do k = 2, nz - 1
          Do j = 2, nyg - 1
            Do i = 2, nxg - 1
              W_(i,j,k) = 0.d0
            End Do
          End Do
        End Do

        ! Pressure (cell centers)
        Do k = lbound(P_, 3), ubound(P_, 3)
          Do j = lbound(P_, 2), ubound(P_, 2)
            Do i = lbound(P_, 1), ubound(P_, 1)
              xloc = xg(i) - xc
              yloc = yg(j) - yc
              r = sqrt(xloc**2 + yloc**2)

              p_exterior = 0d0
              p_interior = 0d0

              If (r < r_ext_smooth_start) Then
                p_exterior = 0.5d0 * om**2 * r**2
              ElseIf (r >= r_ext_smooth_start .and. r <= r_ext_smooth_end) Then
                p_exterior = 0.5d0 * om**2 * r**2 &
                  * 0.5d0 * (1d0 + cos(pi * (r - r_ext_smooth_start) / (r_ext_smooth_end - r_ext_smooth_start)))
              End If

              If (r > 0d0) Then
                If (r >= r_int_smooth_1_start .and. r < r_int_smooth_1_end) Then
                  p_interior = (om**2 * r1**4 / (r2**2 - r1**2)**2 * &
                    (-0.5d0 * r2**4 / r**2 - 2d0 * r2**2 * log(r) + 0.5d0 * r**2 + 2d0 * r2**2 * log(r2)) + p_jump) &
                    * 0.5d0 * (1d0 + cos(pi * (r - r_int_smooth_1_end) / (r_int_smooth_1_start - r_int_smooth_1_end)))
                ElseIf (r >= r_int_smooth_1_end .and. r <= r_int_smooth_2_start) Then
                  p_interior = om**2 * r1**4 / (r2**2 - r1**2)**2 * &
                    (-0.5d0 * r2**4 / r**2 - 2d0 * r2**2 * log(r) + 0.5d0 * r**2 + 2d0 * r2**2 * log(r2)) + p_jump
                ElseIf (r >= r_int_smooth_2_start .and. r <= r_int_smooth_2_end) Then
                  p_interior = (om**2 * r1**4 / (r2**2 - r1**2)**2 * &
                    (-0.5d0 * r2**4 / r**2 - 2d0 * r2**2 * log(r) + 0.5d0 * r**2 + 2d0 * r2**2 * log(r2)) + p_jump) &
                    * 0.5d0 * (1d0 + cos(pi * (r - r_int_smooth_2_start) / (r_int_smooth_2_end - r_int_smooth_2_start)))
                End If
              End If

              P_(i,j,k) = p_exterior * Hc_exterior_(i,j,k) + p_interior * (1d0 - Hc_exterior_(i,j,k))
            End Do
          End Do
        End Do

    End Select

  End Subroutine exact_masked_fields

End Module immersed_boundary_geometry

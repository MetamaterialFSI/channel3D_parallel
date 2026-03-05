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
        nxb = 0
        nzb = 0
        nb = 0
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case ('center_wall') ! Static planar IB wall centered at y = 1
        nb = nxb * nzb
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case ('double_cylinders_z') ! Double concentric cylinders with axis parallel to z
        r1 = body_param_1
        r2 = body_param_2
        dzb = real(Lzp / nzb, 8)
        dxb = dzb * dx / dz
        nxb1 = int(2 * 3.14159 * r1 / dxb)
        nxb2 = int(2 * 3.14159 * r2 / dxb)
        nxb = nxb1 + nxb2
        nb = nxb * nzb

      Case ('standing_wave') ! Top and bottom wall undergoing standing wave motion
        nb = 2 * nxb * nzb
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case ('traveling_wave_x') ! Top and bottom wall undergoing traveling wave motion in the x-direction
        nb = 2 * nxb * nzb
        dxb = real(Lxp / nxb, 8)
        dzb = real(Lzp / nzb, 8)

      Case ('traveling_wave_z') ! Top and bottom wall undergoing traveling wave motion in the z-direction
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
    Real   (Int64) :: a1, a2, r1, r2, xc, yc, theta, dsb1, dsb2, phi, amp

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
        nb_end = 0         ! Initialize to the lowest possible index
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + (j-1) * nxb
            xb(k) = (real(i,8) - 0.5d0) * dxb
            yb(k) = 2.0d0
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
        moving_body = .False. ! should be False, but set to True for speed test
        moving_z_flag = .False.
        nb_start = nb + 1  ! Initialize to an invalid value (beyond the max index)
        nb_end = 0         ! Initialize to the lowest possible index
        y_ref_index = 1 ! The grid has to be uniform for this case, so it doesn't matter what y_ref_index is

        r1 = body_param_1
        r2 = body_param_2
        xc = 1.0d0
        yc = 1.0d0
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

      Case ('standing_wave') ! Top and bottom wall undergoing standing wave motion in x-direction
        If ( body_param_1 > min_buffer_width ) Stop 'Error: IB amplitude is bigger than the minimum buffer width'
        moving_body = .True.
        moving_z_flag = .False.

        ! Scalar arrays. Arrange such that the points treated by one partition are contiguous
        ! (i.e., fall between an nb_start and nb_end)
        nb_start = nb + 1 ! Initialize to an invalid value (beyond the max index)
        nb_end = 0        ! Initialize to the lowest possible index
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + 2 * nxb * (j - 1)
            xb(k)       = (real(i,8) - 0.5d0) * dxb
            xb(k + nxb) = (real(i,8) - 0.5d0) * dxb
            yb(k)       =       body_param_1 * sin(2d0 * pi * body_param_3 * xb(k) / Lxp) * cos(body_param_2 * t)
            yb(k + nxb) = 2d0 + body_param_1 * sin(2d0 * pi * body_param_3 * xb(k) / Lxp) * cos(body_param_2 * t)
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
            ub(nb + k)       = -body_param_1 * body_param_2 * sin(2d0 * pi * body_param_3 * xb(k) / Lxp) * sin(body_param_2 * t)
            ub(nb + k + nxb) = -body_param_1 * body_param_2 * sin(2d0 * pi * body_param_3 * xb(k) / Lxp) * sin(body_param_2 * t)

            tangents_1(nb + k)       = body_param_1 * 2d0 * pi * body_param_3 / Lxp &
              * cos(2d0 * pi * body_param_3 * xb(k) / Lxp) * cos(body_param_2 * t)
            tangents_1(nb + k + nxb) = body_param_1 * 2d0 * pi * body_param_3 / Lxp &
              * cos(2d0 * pi * body_param_3 * xb(k) / Lxp) * cos(body_param_2 * t)
            ! scale to unit vectors
            tangents_1(k)            = 1 / sqrt( 1 + tangents_1(nb + k) ** 2 )
            tangents_1(k + nxb)      = 1 / sqrt( 1 + tangents_1(nb + k + nxb) ** 2 )
            tangents_1(nb + k)       = tangents_1(nb + k) / sqrt( 1 + tangents_1(nb + k) ** 2 )
            tangents_1(nb + k + nxb) = tangents_1(nb + k + nxb) / sqrt( 1 + tangents_1(nb + k + nxb) ** 2 )

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
        moving_body = .True.

        If (t < body_ramp_up_time .and. body_ramp_up_time > 0) Then
          amp = body_param_1 * t / body_ramp_up_time
        Else
          amp = body_param_1
        End If

        ! Phase difference between top and bottom wave motion
        phi = pi

        ! Scalar arrays. Arrange such that the points treated by one partition are contiguous
        ! (i.e., fall between an nb_start and nb_end)
        nb_start = nb + 1 ! Initialize to an invalid value (beyond the max index)
        nb_end = 0        ! Initialize to the lowest possible index
        Do j = 1, nzb
          Do i = 1, nxb
            k = i + 2 * nxb * (j - 1)
            xb(k)       = (real(i,8) - 0.5d0) * dxb
            xb(k + nxb) = (real(i,8) - 0.5d0) * dxb
            yb(k)       =       amp / (body_param_2 * body_param_3) * sin(body_param_3 * (xb(k) - body_param_2 * t))
            yb(k + nxb) = 2d0 + amp / (body_param_2 * body_param_3) * sin(body_param_3 * (xb(k) - body_param_2 * t) + phi)
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
            ub(nb + k)       = amp * cos(body_param_3 * (xb(k) - body_param_2 * t))
            ub(nb + k + nxb) = amp * cos(body_param_3 * (xb(k) - body_param_2 * t) + phi)

            tangents_1(nb + k)       = amp * body_param_3 &
              * cos(body_param_3 * (xb(k) - body_param_2 * t))
            tangents_1(nb + k + nxb) = amp * body_param_3 &
              * cos(body_param_3 * (xb(k) - body_param_2 * t) + phi)
            ! scale to unit vectors
            tangents_1(k)            = 1 / sqrt( 1 + tangents_1(nb + k) ** 2 )
            tangents_1(k + nxb)      = 1 / sqrt( 1 + tangents_1(nb + k + nxb) ** 2 )
            tangents_1(nb + k)       = tangents_1(nb + k) / sqrt( 1 + tangents_1(nb + k) ** 2 )
            tangents_1(nb + k + nxb) = tangents_1(nb + k + nxb) / sqrt( 1 + tangents_1(nb + k + nxb) ** 2 )

            tangents_2(2 * nb + k      ) = 1d0
            tangents_2(2 * nb + k + nxb) = -1d0

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
              a1 = sqrt((xb(k) - (xb(1 + 2 * nxb * (j - 1)) + Lxp)) ** 2   + (yb(k) - yb(1 + 2 * nxb * (j - 1))) ** 2)
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

End Module immersed_boundary_geometry

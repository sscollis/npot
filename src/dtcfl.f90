!=============================================================================!
        subroutine dtcfl(nx, ny, u, v, c, m1l, m2l, n1l, n2l, detJ, &
                         dxi, deta, cflmax, loctime, iter, dtl)
!  
!  Compute the CFL and time step
!  
!=============================================================================!
        use constants
        implicit none

        integer :: nx, ny
        real :: u(nx,ny), v(nx,ny), c(nx,ny)
        real :: m1l(nx,ny), m2l(nx,ny), n1l(nx,ny), n2l(nx,ny), detJ(nx,ny)
        real :: dtl(nx,ny)
        real :: dxi, deta, cflmax
        integer :: loctime, iter
        
        real :: u1l, u2l, cl, dm, dn, cfl, cfll, dtlmax, delt
        integer :: i, j, ji, il, jl
        
        real :: cfla, cflc, cflal, cflcl
        integer :: ia, ja, ic, jc

        real :: vol, sxj, syj, sxi, syi
!=============================================================================!

!.... new way of computing CFL

        cfla = 0.0
        cflc = 0.0

        do j = 1, ny
          do i = 1, nx
          
            sxj =  n2l(i,j) * dxi / detJ(i,j)
            syj = -m2l(i,j) * dxi / detJ(i,j)
            
            sxi = -n1l(i,j) * deta / detJ(i,j)
            syi =  m1l(i,j) * deta / detJ(i,j)

            vol = dxi * deta / detJ(i,j)
            
            u1l = u(i,j)
            u2l = v(i,j)
            
            cl = c(i,j)
            
!.... acoustic CFL

            cflal = (abs( u1l * sxi + u2l * syi ) + &
                     abs( u1l * sxj + u2l * syj ) + &
                     cl * sqrt( sxi**2 + two * sxi * syi + syi**2 + &
                                sxj**2 + two * sxj * syj + syj**2 ) ) / vol

!.... use the acoustic CFL to scale the local time step

            dtl(i,j) = cflal  ! one + sqrt(detJ)
            
            if ( cflal .gt. cfla ) then
              cfla = cflal
              ia = i
              ja = j
            end if

!.... convective CFL

            cflcl = (abs( u1l * sxi + u2l * syi ) + &
                     abs( u1l * sxj + u2l * syj ) ) / vol

            if ( cflcl .gt. cflc) then
              cflc = cflcl
              ic = i
              jc = j
            end if
                                
          end do        ! loop on i
        end do          ! loop on j

        cfl = cfla              ! use the acoustic CFL
        
        delt = cflmax / cfl
        if (loctime.eq.0) then
          dtl  = delt
        else
          dtl  = cflmax / dtl
!         dtl  = delt * (one + sqrt(m1l(ja,ia) * n2l(ja,ia) - &
!                m2l(ja,ia) * n1l(ja,ia))) / dtl
        end if
                   
        if (.false.) then
          write(*,"  ('Maximum Acoustic CFL/dt   = ',1pe13.6, &
                & ' at (',i3,',',i3,')')") cfla, ia, ja
          write(*,"  ('Maximum Convective CFL/dt = ',1pe13.6, &
                & ' at (',i3,',',i3,')')") cflc, ic, jc
          write(*,"(/,'Delt = ',1pe13.6,'  CFL  = ',1pe13.6)") &
                delt, cfl*delt
          write(*,  "('CFLa = ',1pe13.6,'  CFLc = ',1pe13.6,/)") &
                cfla*delt, cflc*delt
        end if

!.... find the maximum dt

        if (loctime.eq.1) then
          dtlmax = zero
          do j = 1, ny
            do i = 1, nx
              if (dtl(i,j).gt.dtlmax) then
                ia = i
                ja = j
                dtlmax = dtl(i,j)
              end if
            end do
          end do
          write(80,"(i5,1x,1pe13.6,1x,2(1pe13.6,1x,i5,1x,i5,1x))") &
                iter, delt, dtlmax, ia, ja
        end if
        
        return
        end

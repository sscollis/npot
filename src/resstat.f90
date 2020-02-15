!============================================================================!
        subroutine resstat( )
!============================================================================!
        use global
        use stencil
        implicit none
!============================================================================!

        resmax = zero
        resmin = infty
        norm = zero
        im = 0; jm = 0; in = 0; jn = 0
!!$omp parallel do private(i,j,resl) reduction(+: norm) reduction(max: resmax)
        do j = 1, ny
          do i = 1, nx
!           resl = res(i,j)**2 / (sigma * omega)**2
            resl = res(i,j)**2 / (dtl(i,j) * sigma * omega)**2
            norm = norm + resl
!           resmax = max(resl,resmax)
            if (sqrt(resl) .gt. resmax) then
              im = i
              jm = j
              resmax = sqrt(resl)
            end if
            if (sqrt(resl) .lt. resmin .and. sqrt(resl) .gt. zero) then
              in = i
              jn = j
              resmin = sqrt(resl)
            end if
          end do
        end do
        norm = sqrt( norm / real(nx*ny) )
!       resmax = sqrt( resmax )
!       resmin = sqrt( resmin )

        return
        end

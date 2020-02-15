!=============================================================================!
        subroutine smoother(rl, vl, eps_e, nx, ny) 
!  
!       Fourth order explicit smoother for the NS equations.
!       From: Computational Fluid Mechanics and Heat Transfer
!             D.A. Anderson, J. C. Tannehill, Rl H. Pletcher
!             Page:  450, 495 
!
!       Written: 6-9-95
!
!       Revised: 7-9-96
!
!       Notes:   Using symmetry on left and bottom boundaries
!=============================================================================!
        use stencil
        use constants
        implicit none

        integer :: nx, ny
        real :: rl(nx,ny), vl(nx,ny), eps_e

        real :: buff(nx,ny)
        integer :: i, j
!=============================================================================!

!.... fourth-order dissipation with reduced order near the boundaries

        do j = 1, ny
          do i = 1, nx
            buff(i,j) = -1.0
          end do
        end do

!.... \eta direction
        
        do i = 1, nx

          j = 3
          rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                      vl(i,j-2) - &
                                4.0 * vl(i,j-1) + &
                                6.0 * vl(i,j  ) - &
                                4.0 * vl(i,j+1) + &
                                1.0 * vl(i,j+2) )

          j = 2
          rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                -1.0 * vl(i,j-1) + &
                                 2.0 * vl(i,j  ) - &
                                 1.0 * vl(i,j+1) )
        end do

        !$omp parallel do private(i,j)
        do j = 4, ny-3
          do i = 1, nx
            rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                fa1 * vl(i,j-3) + &
                                fa2 * vl(i,j-2) + &
                                fa3 * vl(i,j-1) + &
                                fa4 * vl(i,j  ) + &
                                fa5 * vl(i,j+1) + &
                                fa6 * vl(i,j+2) + &
                                fa7 * vl(i,j+3) )
          end do
        end do
        
        do i = 1, nx

          j = ny-2
          rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                      vl(i,j-2) - &
                                4.0 * vl(i,j-1) + &
                                6.0 * vl(i,j  ) - &
                                4.0 * vl(i,j+1) + &
                                1.0 * vl(i,j+2) )

          j = ny-1
          rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                -1.0 * vl(i,j-1) + &
                                 2.0 * vl(i,j  ) - &
                                 1.0 * vl(i,j+1) )
        end do
                                          
!.... \xi direction

        !$omp parallel do private(i,j)
        do j = 1, ny

           i = 3
           rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                      vl(i-2,j) - &
                                4.0 * vl(i-1,j) + &
                                6.0 * vl(i  ,j) - &
                                4.0 * vl(i+1,j) + &
                                1.0 * vl(i+2,j) )

           i = 2
           rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                -1.0 * vl(i-1,j) + &
                                 2.0 * vl(i,  j) - &
                                 1.0 * vl(i+1,j) )

           do i = 4, nx-3
             rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                 fa1 * vl(i-3,j) + &
                                 fa2 * vl(i-2,j) + &
                                 fa3 * vl(i-1,j) + &
                                 fa4 * vl(i  ,j) + &
                                 fa5 * vl(i+1,j) + &
                                 fa6 * vl(i+2,j) + &
                                 fa7 * vl(i+3,j) )
           end do

           i = nx-2
           rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                      vl(i-2,j) - &
                                4.0 * vl(i-1,j) + &
                                6.0 * vl(i  ,j) - &
                                4.0 * vl(i+1,j) + &
                                1.0 * vl(i+2,j) )

           i = nx-1
           rl(i,j) = rl(i,j) + eps_e * buff(i,j) * ( &
                                -1.0 * vl(i-1,j) + &
                                 2.0 * vl(i,  j) - &
                                 1.0 * vl(i+1,j) )
        end do

        return
        end


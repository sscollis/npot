function field(il,jl)
  use global
  implicit none

  real :: field
  integer :: il, jl

  field = xy(1,il,jl) + pt5*circ/pi*(pi-atan2(xy(2,il,jl),-(xy(1,il,jl)-pt5)))
! field = pt5*circ/pi*(pi-atan2(xy(2,il,jl),-(xy(1,il,jl)-pt5)))

  if (jl .eq. 1) then
    if (il.ge.1 .and. il.le.na) then
      field = xy(1,il,jl) + pt5*circ/pi*(pi-atan2(xy(2,il,jl), &
            -(xy(1,il,jl)-pt5))+two*pi)
!     field = pt5*circ/pi*(pi-atan2(xy(2,il,jl),-(xy(1,il,jl)-pt5))+two*pi)
    end if
  end if

end function field
  
!============================================================================!
        subroutine initial
!============================================================================!
        use global
        implicit none

        real :: xil, etal, phiu, phil
        integer :: iu, il

        real, external :: field
!============================================================================!

!.... The conformal mapping for the parabolic cylinder is given by

!....    x = 1/2 * ( xi^2 - eta^2 ) + 1/2
!....    y = xi * eta

        if (restart .eq. 0) then             ! uniform flow in the +x direction
          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              phi(i,j) = field(i,j)
            end do
          end do
          j = 1
          do i = 1, na
            iu = nx - i + 1 
            phil = field(i,j)
            phiu = field(iu,j)
            write(55,"(5(1pe13.6,1x))") xy(1,i,1), phil, phiu, phiu-phil
          end do
          !call flush(55)
        else if (restart .eq. 1) then        ! restart from a previous solution
          open(20,file='restart.dat',form='unformatted',status='old',err=10)
          read(20) phi
          close(20)
        else if (restart .eq. 2) then        ! attempt to add mach correction
          beta = sqrt(one - Ma**2)
          !$omp parallel do private(i,etal,xil,u,v,t,rho)
          do j = 1, ny
            do i = 1, nx
              etal = sqrt(pt5 - xy(1,i,j) + pt5*sqrt((two*xy(1,i,j)-one)**2 + &
                     four * beta**2 * xy(2,i,j)**2))
              xil = beta * xy(2,i,j) / etal
              phi(i,j) = xy(1,i,j) + (pt5*(xil**2-etal**2 )+etal - &
                         xy(1,i,j))/beta**2
              u = one + ((etal**2 - etal + xil**2) / &
                  (xil**2 + etal**2) - one) / beta**2
              v = xil / (xil**2 + etal**2) / beta
              t = one - pt5*gamma1*Ma**2*( u**2 + v**2 - one )
              rho = t**(one/gamma1)
            end do
          end do
        else if (restart .eq. 3) then        ! works for pcyl
          !$omp parallel do private(i,etal,xil,u,v,t,rho)
          do j = 1, ny
            do i = 1, nx
              etal = sqrt(pt5-xy(1,i,j)+pt5*sqrt((two*xy(1,i,j)-one)**2 + &
                     four*xy(2,i,j)**2))
              xil = xy(2,i,j) / etal
              phi(i,j) = pt5 * ( xil**2 - etal**2 ) + etal
              phis(i,j) = phi(i,j)
              u = (etal**2 - etal + xil**2) / (xil**2 + etal**2)
              v = xil / (xil**2 + etal**2)
              t = one - pt5*gamma1*Ma**2*( u**2 + v**2 - one )
              rho = t**(one/gamma1)
            end do
          end do
        else if (restart .eq. 4) then
          !$omp parallel do private(i,etal,xil,u,v,t,rho)
          do j = 1, ny
            do i = 1, nx
              phi(i,j) = xy(1,i,j) + pt5 * log( xy(1,i,j)**2 + xy(2,i,j)**2 )
            end do
          end do
        end if

        return
10      call error('initial$','Unable to open restart file.$')
        end

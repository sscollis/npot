!============================================================================!
        subroutine itrbc( )
!============================================================================!
        use global
        use stencil
        implicit none

!.... Riemann boundary variables

        real, allocatable :: etab(:), xib(:), ub(:), vb(:), cb(:)
        real, allocatable :: bn1(:), bn2(:), bn(:), vn(:)
        real, allocatable :: rhob(:), tb(:), pb(:), phib(:)

        integer :: ib
        real, external :: field
!============================================================================!

!.... top-bottom periodic boundary condition

        if (yper) then
          phi(:,ny) = phi(:,1)
        end if

!.... hard boundary condition on top side

        if (top.eq.0) then
          if (restart .eq. 3) then
            allocate( etab(nx), xib(nx) )
            etab = sqrt(pt5 - xy(1,:,ny) + &
                   pt5*sqrt((two*xy(1,:,ny)-one)**2 + four * xy(2,:,ny)**2))
            xib = xy(2,:,ny) / etab
            phi(1:bx,ny) = pt5 * ( xib(1:bx)**2 - etab(1:bx)**2 ) + &
                           etab(1:bx)
            deallocate( etab, xib )
          else if (restart.eq.4) then
            phi(1:bx,ny) = xy(1,1:bx,ny) + pt5 * log( xy(1,1:bx,ny)**2 + &
                           xy(2,1:bx,ny)**2 )
          else 
            do i = 1, bx
              phi(i,ny) = field(i,ny)
            end do
          end if
        end if

!.... Reimann on top boundary

!!$     if (top.eq.1) then
!!$       call grad(ndof, nx, ny, phi, g1phi, g2phi, dxi, deta, optx, opty, &
!!$                 xper, yper, lsym, rsym, bsym, tsym, carp)
!!$
!!$       u   = m1 * g1phi + n1 * g2phi 
!!$       v   = m2 * g1phi + n2 * g2phi
!!$       t   = one + gamma1 * Ma**2 * pt5 * ( one - u**2 - v**2 )
!!$       c   = sqrt(t) / Ma
!!$       rho = t ** ( one / gamma1 )
!!$       p   = rho * t / (gamma * Ma**2) 
!!$
!!$       allocate( bn(nx), bn1(nx), bn2(nx), uint(nx), vint(nx), cint(nx) )
!!$       allocate( etab(nx), xib(nx), ub(nx), vb(nx), cb(nx), phib(nx) )
!!$       etab = sqrt(pt5 - xy(1,:,ny) + pt5*sqrt((two*xy(1,:,ny)-one)**2 + &
!!$              four * xy(2,:,ny)**2))
!!$       xib  = xy(2,:,ny) / etab
!!$       ub   = (etab**2 - etab + xib**2) / (xib**2 + etab**2)
!!$       vb   = xib / (xib**2 + etab**2)
!!$       cb   = sqrt(one - pt5*gamma1*Ma**2*( ub**2 + vb**2 - one )) / Ma
!!$
!!$       bn(:)  = one / sqrt( n1(:,ny)**2 + n2(:,ny)**2 )
!!$       bn1(:) = n1(:,ny) * bn
!!$       bn2(:) = n2(:,ny) * bn
!!$
!!$       uint = two * u(:,ny-1) - u(:,ny-2)
!!$       vint = two * v(:,ny-1) - v(:,ny-2)
!!$       cint = two * c(:,ny-1) - c(:,ny-2)
!!$       phib = ( pt5 * deta * bn * &
!!$              ( bn1 * ub + bn2 * vb - two * cb / gamma1 + &
!!$                bn1 * uint + bn2 * vint + two * cint / gamma1 ) + &
!!$                gc2 * phi(:,ny-1) + gc3 * phi(:,ny-2) + &
!!$                gc4 * phi(:,ny-3) + gc5 * phi(:,ny-4) ) / (-gc1)
!!$       phi(1:bx,ny) = phib(1:bx)
!!$       deallocate( etab, xib, ub, vb, cb, phib )
!!$       deallocate( bn, bn1, bn2, uint, vint, cint )
!!$     end if

!.... bottom boundary for inviscid wall

        j = 1
        if (wakecut) then

           call grad(ndof, nx, ny, phi, gphi, dxi, deta, optx, opty, &
                 xper, yper, lsym, rsym, bsym, tsym, carp)

!!$            do i = 1, na
!!$           ib = nx - i + 1
!!$           phi(i,1) = ( -m(4,i,1)/deta*(gc2 * phi(i,2) + gc3 * phi(i,3) + &
!!$                                         gc4 * phi(i,4) + gc5 * phi(i,5) ) + &
!!$                           m(4,ib,1)/deta * ( gc1 * (-circ) + &
!!$                                         gc2 * phi(ib,2) + gc3 * phi(ib,3) + &
!!$                                         gc4 * phi(ib,4) + gc5 * phi(ib,5))&
!!$                   +(-m(2,i,1)*gphi(1,i,1) + m(2,ib,1)*gphi(1,ib,1))  &
!!$                   )/(m(4,i,1)/deta*gc1 - m(4,ib,1)/deta*gc1)
!!$         end do
        
!!$            do i = 1, na
!!$           ib = nx - i + 1
!!$           phi(i,1) = ( -m(4,i,1)/deta*(gk2 * phi(i,2) + gk3 * phi(i,3)) + &
!!$                           m(4,ib,1)/deta * ( gk1 * (-circ) + &
!!$                                         gk2 * phi(ib,2) + gk3 * phi(ib,3) )&
!!$                   +(-m(2,i,1)*gphi(1,i,1) + m(2,ib,1)*gphi(1,ib,1))  &
!!$                   )/(m(4,i,1)/deta*gk1 - m(4,ib,1)/deta*gk1)
!!$         end do
!!$
            do i = 1, na
              ib = nx - i + 1
              phi(i,1) = ( -m(4,i,1) * ( gc2 * phi(i,2) + gc3 * phi(i,3) + &
                                         gc4 * phi(i,4) + gc5 * phi(i,5) ) + &
                           m(4,ib,1) * ( gc1 * (-circ) + &
                                         gc2 * phi(ib,2) + gc3 * phi(ib,3) + &
                                         gc4 * phi(ib,4) + gc5 * phi(ib,5))&
                   +(-m(2,i,1)*gphi(1,i,1) + m(2,ib,1)*gphi(1,ib,1))*deta   &
                   )/(m(4,i,1)*gc1 - m(4,ib,1)*gc1)
            end do

!!$
!!$         phi(na+1:nb-1,1) = -( gk2 * phi(na+1:nb-1,2) + &
!!$                               gk3 * phi(na+1:nb-1,3) ) / gk1

            phi(na+1:nb-1,1) = -( gc2 * phi(na+1:nb-1,2) + &
                                  gc3 * phi(na+1:nb-1,3) + &
                                  gc4 * phi(na+1:nb-1,4) + &
                                  gc5 * phi(na+1:nb-1,5) ) / gc1

            do i = nb, nx
              phi(i,1) = phi(na-(i-nb),1) - circ
            end do

         else
          if (bottom.eq.0) then
            phi(1:bx,1) = -( gc2 * phi(1:bx,2) + &
                             gc3 * phi(1:bx,3) + &
                             gc4 * phi(1:bx,4) + &
                             gc5 * phi(1:bx,5) ) / gc1
          end if
        end if

!.... left-right periodic boundary condition

        if (xper) then
          phi(nx,:) = phi(1,:)
        end if

!.... hard boundary condition on left side

        if (left.eq.0) then
          do j = 1, by
            phi(1,j) = field(1,j)
          end do
        end if

!.... Extrapolation of left boundary

        if (left.eq.2) then
          phi(1,1:by) = two * phi(2,1:by) - phi(3,1:by)
        end if

!.... hacked wake test

        if (right.eq.-2) then
          do j = 1, ny
            phi(nx,j) = phi(1,j)
          end do
        end if

!.... hard boundary condition on right side

        if (right.eq.0) then
          if (restart .eq. 3) then
            allocate( etab(ny), xib(ny) )
            etab = sqrt(pt5 - xy(1,nx,:) + pt5*sqrt((two*xy(1,nx,:)-one)**2 + &
                   four * xy(2,nx,:)**2))
            xib = xy(2,nx,:) / etab
            phi(nx,1:by) = pt5 * ( xib(1:by)**2 - etab(1:by)**2 ) + &
                           etab(1:by)
            deallocate( etab, xib )
          else if (restart.eq.4) then
            phi(nx,1:by) = xy(1,nx,1:by) + pt5 * log( xy(1,nx,1:by)**2 + &
                           xy(2,nx,1:by)**2 )
          else
            do j = 1, by
              phi(nx,j) = field(nx,j)
            end do
          end if
        end if

!.... Reimann on right boundary

!!$     if (right.eq.1) then
!!$       call grad(ndof, nx, ny, phi, g1phi, g2phi, dxi, deta, optx, opty, &
!!$                 xper, yper, lsym, rsym, bsym, tsym, carp)
!!$
!!$       u   = m1 * g1phi + n1 * g2phi 
!!$       v   = m2 * g1phi + n2 * g2phi
!!$       t   = one + gamma1 * Ma**2 * pt5 * ( one - u**2 - v**2 )
!!$       c   = sqrt(t) / Ma
!!$       rho = t ** ( one / gamma1 )
!!$       p   = rho * t / (gamma * Ma**2)
!!$
!!$       allocate( bn(ny), bn1(ny), bn2(ny), uint(ny), vint(ny), cint(ny) )
!!$       allocate( etab(ny), xib(ny), ub(ny), vb(ny), cb(ny), phib(ny) )
!!$       etab = sqrt(pt5 - xy(1,nx,:) + pt5*sqrt((two*xy(1,nx,:)-one)**2 + &
!!$              four * xy(2,nx,:)**2))
!!$       xib  = xy(2,nx,:) / etab
!!$       ub   = (etab**2 - etab + xib**2) / (xib**2 + etab**2)
!!$       vb   = xib / (xib**2 + etab**2)
!!$       cb   = sqrt(one - pt5*gamma1*Ma**2*( ub**2 + vb**2 - one )) / Ma
!!$
!!$       bn(:)  = one / sqrt( m1(nx,:)**2 + m2(nx,:)**2 )
!!$       bn1(:) = m1(nx,:) * bn
!!$       bn2(:) = m2(nx,:) * bn
!!$       uint = two * u(nx-1,:) - u(nx-2,:)
!!$       vint = two * v(nx-1,:) - v(nx-2,:)
!!$       cint = two * c(nx-1,:) - c(nx-2,:)
!!$       phib = ( pt5 * dxi * bn * &
!!$              ( bn1 * ub + bn2 * vb - two * cb / gamma1 + &
!!$                bn1 * uint + bn2 * vint + two * cint / gamma1 ) + &
!!$                gc2 * phi(nx-1,:) + gc3 * phi(nx-2,:) + &
!!$                gc4 * phi(nx-3,:) + gc5 * phi(nx-4,:) ) / (-gc1)
!!$       phi(nx,1:by) = phib(1:by)
!!$       deallocate( etab, xib, ub, vb, cb, phib )
!!$       deallocate( bn, bn1, bn2, uint, vint, cint)
!!$     end if

!.... Extrapolation of right boundary

        if (right.eq.2) then
          phi(nx,1:by) = two * phi(nx-1,1:by) - phi(nx-2,1:by)
        end if

        return
        end

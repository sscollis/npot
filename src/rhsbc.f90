!============================================================================!
        subroutine rhsbc( )
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
          res(:,ny) = zero
        end if

!.... hold the initial values of phi

        if (top.eq.0) then
          do i = 1, bx
            res(i,ny) = phi(i,ny)-field(i,ny)
          end do
        end if
 
        if (top.eq.3.or.top.eq.4) then
          do i = 1, bx
            res(i,ny) = zero
          end do
        end if

!.... Reimann on top boundary

!!$     if (top.eq.1) then
!!$       call grad(ndof, nx, ny, phi, g1phi, g2phi, dxi, deta, optx, opty, &
!!$                 xper, yper, lsym, rsym, bsym, tsym, carp)
!!$       
!!$       u = m1 * g1phi + n1 * g2phi 
!!$       v = m2 * g1phi + n2 * g2phi
!!$       t = one + gamma1 * Ma**2 * pt5 * ( one - u**2 - v**2 )
!!$       c = sqrt(t) / Ma
!!$
!!$       allocate( bn(nx), bn1(nx), bn2(nx), uint(nx), vint(nx), cint(nx) )
!!$       allocate( etab(nx), xib(nx), ub(nx), vb(nx), cb(nx), phib(nx) )
!!$       etab = sqrt(pt5 - x(:,ny) + pt5*sqrt((two*x(:,ny)-one)**2 + &
!!$              four * y(:,ny)**2))
!!$       xib  = y(:,ny) / etab
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
!!$       phib = gc1 * phi(:,ny) + ( pt5 * deta * bn * &
!!$              ( bn1 * ub + bn2 * vb  - two * cb / gamma1 + &
!!$                bn1 * uint + bn2 * vint + two * cint / gamma1 ) + &
!!$              gc2 * phi(:,ny-1) + gc3 * phi(:,ny-2) + &
!!$              gc4 * phi(:,ny-3) + gc5 * phi(:,ny-4) )
!!$       res(1:bx,ny) = phib(1:bx)
!!$       deallocate( etab, xib, ub, vb, cb, phib )
!!$       deallocate( bn, bn1, bn2, uint, vint, cint )
!!$     end if
          
!.... bottom boundary for inviscid wall

        if (wakecut) then

!.... left portion (Continuity of first derivative)

          call grad(ndof, nx, ny, phi, gphi, dxi, deta, optx, opty, &
                    xper, yper, lsym, rsym, bsym, tsym, carp)
          j = 1
          do i = 1, na
            ib = nx - i + 1
!!$          res(i,1) = m(4,i,j)*(phi(i,2)-phi(i,1)) - &
!!$                     m(4,ib,j)*(phi(ib,2)-phi(ib,1))

!!$          res(i,1) =  m(4,i,1)*(gk1*phi(i,1) +gk2*phi(i,2) +gk3*phi(i,3)) - &
!!$                     m(4,ib,1)*(gk1*phi(ib,1)+gk2*phi(ib,2)+gk3*phi(ib,3))+&
!!$            (m(2,i,1)*gphi(1,i,1) - m(2,ib,1)*gphi(1,ib,1))*deta

             res(i,1) = m(4,i,1) * ( gc1 * phi(i,1) + &
                                     gc2 * phi(i,2) + &
                                     gc3 * phi(i,3) + &
                                     gc4 * phi(i,4) + &
                                     gc5 * phi(i,5) ) - &
                        m(4,ib,1)* ( gc1 * phi(ib,1) + &
                                     gc2 * phi(ib,2) + &
                                     gc3 * phi(ib,3) + &
                                     gc4 * phi(ib,4) + &
                                     gc5 * phi(ib,5) )  + &
                      (m(2,i,1)*gphi(1,i,1) - m(2,ib,1)*gphi(1,ib,1))*deta
          end do

!.... center wall portion

          res(na+1:nb-1,1) = gc1 * phi(na+1:nb-1,1) + &
                             gc2 * phi(na+1:nb-1,2) + &
                             gc3 * phi(na+1:nb-1,3) + &
                             gc4 * phi(na+1:nb-1,4) + &
                             gc5 * phi(na+1:nb-1,5)

!!$       res(na+1:nb-1,1) = gk1 * phi(na+1:nb-1,1) + &
!!$                          gk2 * phi(na+1:nb-1,2) + &
!!$                          gk3 * phi(na+1:nb-1,3)

!.... right portion (continuity or jump condition)

          do i = nb, nx
            res(i,1) = phi(i,1) - phi(na-(i-nb),1) + circ
          end do

        else
          if (bottom.eq.0) then
            res(1:bx,1) = gc1 * phi(1:bx,1) + gc2 * phi(1:bx,2) + &
                          gc3 * phi(1:bx,3) + gc4 * phi(1:bx,4) + &
                          gc5 * phi(1:bx,5)
          end if
        end if

!.... left-right periodic boundary condition

        if (xper) then
          res(nx,:) = zero
        end if

!.... hold the initial values of phi

        if (left.eq.0) then
          res(1,1:by) = zero
        end if

        if (left.eq.2) then
          res(1,1:by) = (phi(1,1:by) - ( two * phi(2,1:by) - &
                                         phi(3,1:by) ) )
        end if

!.... hold the initial values of phi

        if (right.eq.0) then
          res(nx,1:by) = zero
        end if
        
        if (right.eq.3.or.right.eq.4) then
          res(nx,1:by) = zero
        end if

!.... Reimann on right boundary

!!$     if (right.eq.1) then
!!$       call grad(ndof, nx, ny, phi, g1phi, g2phi, dxi, deta, optx, opty, &
!!$                 xper, yper, lsym, rsym, bsym, tsym, carp)
!!$       
!!$       u = m1 * g1phi + n1 * g2phi 
!!$       v = m2 * g1phi + n2 * g2phi
!!$       t = one + gamma1 * Ma**2 * pt5 * ( one - u**2 - v**2 )
!!$       c = sqrt(t) / Ma
!!$
!!$       allocate( bn(ny), bn1(ny), bn2(ny), uint(ny), vint(ny), cint(ny) )
!!$       allocate( etab(ny), xib(ny), ub(ny), vb(ny), cb(ny), phib(ny) )
!!$       etab = sqrt(pt5 - x(nx,:) + pt5*sqrt((two*x(nx,:)-one)**2 + &
!!$              four * y(nx,:)**2))
!!$       xib  = y(nx,:) / etab
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
!!$       phib = gc1 * phi(nx,:) + ( pt5 * dxi * bn * &
!!$              ( bn1 * ub + bn2 * vb - two * cb / gamma1 + &
!!$                bn1 * uint + bn2 * vint + two * cint / gamma1 ) + &
!!$               gc2 * phi(nx-1,:) + gc3 * phi(nx-2,:) + &
!!$               gc4 * phi(nx-3,:) + gc5 * phi(nx-4,:) )
!!$       res(nx,1:by) = phib(1:by)
!!$       deallocate( etab, xib, ub, vb, cb, phib )
!!$       deallocate( bn, bn1, bn2, uint, vint, cint)
!!$     end if

!..... Extrapolation on right boundary

        if (right.eq.2) then
          res(nx,1:by) = (phi(nx,1:by) - ( two * phi(nx-1,1:by) - &
                                           phi(nx-2,1:by) ) )
        end if

        return
        end


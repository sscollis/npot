!=============================================================================!
        subroutine grad2 (ndof, nx, ny, v, gv, ggv, dx, dy, &
                          optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)
!
!  Take the second derivative of a 2-D field.
!  Updated to sixth-order accurate differencing on the interior with
!  the option for optimized fourth-order differencing.
!
!  Added g1v as input to improve efficiency
!
!  Revised: 6-28-95
!=============================================================================!
        use stencil
        implicit none
        
        integer :: ndof, nx, ny, optx, opty
        logical :: xper, yper
        logical :: lsym, rsym, bsym, tsym, carp
        real    :: v(ndof,nx,ny), gv(2,ndof,nx,ny), ggv(3,ndof,nx,ny)
        !$sgi distribute v(*,*,block), gv(*,*,*,block), ggv(*,*,*,block)
        real    :: dx, dy

        real, parameter :: zero = 0.0, one = 1.0, pt5 = 0.5, two = 2.0
        real :: dxinv, dyinv, dxsinv, dysinv
        real :: a, b, c, w
        real :: gx1, gx2, gx3, gx4, gx5, gx6
        real :: gy1, gy2, gy3, gy4, gy5, gy6
        real :: dx1, dx2, dx3, dx4, dx5, dx6, dx7
        real :: dy1, dy2, dy3, dy4, dy5, dy6, dy7
        
        integer :: i, j, idof
        real, parameter :: eps = 1.0e-12
        real :: isign
!=============================================================================!

        dxinv  = one / dx
        dyinv  = one / dy
        dxsinv = one / dx**2
        dysinv = one / dy**2

!.... seven point stencil in x

        if (optx.eq.0) then
          c = 1.0 / 60.0
        else if (optx.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gx1 =  -c * dxinv
        gx2 =  -b * dxinv
        gx3 =  -a * dxinv
        gx4 =   a * dxinv
        gx5 =   b * dxinv
        gx6 =   c * dxinv

        if (optx.eq.0) then
          c = 1.0 / 90.0
        else if (optx.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = -(3.0*w**2 - 16.0*Sin(w/2.0)**2 + Sin(w)**2) / &
               (12.0*(-15.0*Sin(w/2.0)**2 + 6.0*Sin(w)**2 -  &
               Sin(3.0*w/2.0)**2))
        end if
        
        a = (4.0 + 45.0 * c) / 3.0
        b = -(1.0 + 72.0 * c) / 12.0

        dx1 =  c * dxsinv
        dx2 =  b * dxsinv
        dx3 =  a * dxsinv
        dx4 = -2.0 * ( a + b + c ) * dxsinv
        dx5 =  a * dxsinv
        dx6 =  b * dxsinv
        dx7 =  c * dxsinv

!.... seven point stencil in y

        if (opty.eq.0) then
          c = 1.0 / 60.0
        else if (opty.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = (w/2.0 - 2.0*Sin(w)/3.0 + Sin(2.0*w)/12.0) / &
              (5.0*Sin(w) - 4.0*Sin(2.0*w) + Sin(3.0*w))
        end if
        
        a = (2.0 + 15.0 * c) / 3.0
        b = -(1.0 + 48.0 * c) / 12.0

        gy1 =  -c * dyinv
        gy2 =  -b * dyinv
        gy3 =  -a * dyinv
        gy4 =   a * dyinv
        gy5 =   b * dyinv
        gy6 =   c * dyinv

        if (opty.eq.0) then
          c = 1.0 / 90.0
        else if (opty.eq.-1) then
          c = 0.0
        else
          w = 2.0 * 3.1415926535897932385e+0 / 12.0
          c = -(3.0*w**2 - 16.0*Sin(w/2.0)**2 + Sin(w)**2) / &
               (12.0*(-15.0*Sin(w/2.0)**2 + 6.0*Sin(w)**2 -  &
               Sin(3.0*w/2.0)**2))
        end if
        
        a = (4.0 + 45.0 * c) / 3.0
        b = -(1.0 + 72.0 * c) / 12.0

        dy1 =  c * dysinv
        dy2 =  b * dysinv
        dy3 =  a * dysinv
        dy4 = -2.0 * ( a + b + c ) * dysinv
        dy5 =  a * dysinv
        dy6 =  b * dysinv
        dy7 =  c * dysinv

!=============================================================================!
!.... compute the second derivative in x
!=============================================================================!

        if (xper) then

!$doacross
!$omp parallel do
          do j = 1, ny

          ggv(1,:,1,j)       = ( dx1 * v(:,nx-3,j)   + &
                                dx2 * v(:,nx-2,j)   + &
                                dx3 * v(:,nx-1,j)   + &
                                dx4 * v(:,1,j)      + &
                                dx5 * v(:,2,j)      + &
                                dx6 * v(:,3,j)      + &
                                dx7 * v(:,4,j)      ) 
  
          ggv(1,:,2,j)       = ( dx1 * v(:,nx-2,j)   + &
                                dx2 * v(:,nx-1,j)   + &
                                dx3 * v(:,1,j)      + &
                                dx4 * v(:,2,j)      + &
                                dx5 * v(:,3,j)      + &
                                dx6 * v(:,4,j)      + &
                                dx7 * v(:,5,j)      ) 
  
          ggv(1,:,3,j)       = ( dx1 * v(:,nx-1,j)   + &
                                dx2 * v(:,1,j)      + &
                                dx3 * v(:,2,j)      + &
                                dx4 * v(:,3,j)      + &
                                dx5 * v(:,4,j)      + &
                                dx6 * v(:,5,j)      + &
                                dx7 * v(:,6,j)      ) 
  
          ggv(1,:,nx-2,j)    = ( dx1 * v(:,nx-5,j)   + &
                                dx2 * v(:,nx-4,j)   + &
                                dx3 * v(:,nx-3,j)   + &
                                dx4 * v(:,nx-2,j)   + &
                                dx5 * v(:,nx-1,j)   + &
                                dx6 * v(:,1,j)      + &
                                dx7 * v(:,2,j)      ) 
  
          ggv(1,:,nx-1,j)    = ( dx1 * v(:,nx-4,j)   + &
                                dx2 * v(:,nx-3,j)   + &
                                dx3 * v(:,nx-2,j)   + &
                                dx4 * v(:,nx-1,j)   + &
                                dx5 * v(:,1,j)      + &
                                dx6 * v(:,2,j)      + &
                                dx7 * v(:,3,j)   ) 
  
          ggv(1,:,nx,j) = ggv(1,:,1,j)

          end do
          
        else
        
!$doacross
!$omp parallel do
          do j = 1, ny

          ggv(1,:,1,j)       = ( dd1 * v(:,1,j) + &
                                dd2 * v(:,2,j) + &
                                dd3 * v(:,3,j) + &
                                dd4 * v(:,4,j) + &
                                dd5 * v(:,5,j) ) * dxsinv
  
          ggv(1,:,2,j)       = ( db1 * v(:,1,j) + &
                                db2 * v(:,2,j) + &
                                db3 * v(:,3,j) + &
                                db4 * v(:,4,j) + &
                                db5 * v(:,5,j) ) * dxsinv
  
          ggv(1,:,3,j)       = ( da1 * v(:,1,j) + &
                                da2 * v(:,2,j) + &
                                da3 * v(:,3,j) + &
                                da4 * v(:,4,j) + &
                                da5 * v(:,5,j) ) * dxsinv
  
          ggv(1,:,nx-2,j)    = ( da1 * v(:,nx-4,j) + &
                                da2 * v(:,nx-3,j) + &
                                da3 * v(:,nx-2,j) + &
                                da4 * v(:,nx-1,j) + &
                                da5 * v(:,nx,j)   ) * dxsinv
  
          ggv(1,:,nx-1,j)    =   ( db1 * v(:,nx,j)   + &
                                  db2 * v(:,nx-1,j) + &
                                  db3 * v(:,nx-2,j) + &
                                  db4 * v(:,nx-3,j) + &
                                  db5 * v(:,nx-4,j) ) * dxsinv
  
          ggv(1,:,nx,j)      =   ( dd1 * v(:,nx,j)   + &
                                  dd2 * v(:,nx-1,j) + &
                                  dd3 * v(:,nx-2,j) + &
                                  dd4 * v(:,nx-3,j) + &
                                  dd5 * v(:,nx-4,j) ) * dxsinv

          end do

        end if
        
!.... implement symmetry conditions

        if (lsym) then
!$doacross local( idof, isign )
!$omp parallel do private( idof, isign )
          do j = 1, ny
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
            
            ggv(1,idof,1,j)    = ( isign * dx1 * v(idof,4,j)    + &
                                  isign * dx2 * v(idof,3,j)     + &
                                  isign * dx3 * v(idof,2,j)     + &
                                  dx4 * v(idof,1,j)             + &
                                  dx5 * v(idof,2,j)             + &
                                  dx6 * v(idof,3,j)             + &
                                  dx7 * v(idof,4,j)             ) 
    
            ggv(1,idof,2,j)    = ( isign * dx1 * v(idof,3,j)    + &
                                  isign * dx2 * v(idof,2,j)     + &
                                  dx3 * v(idof,1,j)             + &
                                  dx4 * v(idof,2,j)             + &
                                  dx5 * v(idof,3,j)             + &
                                  dx6 * v(idof,4,j)             + &
                                  dx7 * v(idof,5,j)             ) 
    
            ggv(1,idof,3,j)    = ( isign * dx1 * v(idof,2,j)    + &
                                  dx2 * v(idof,1,j)             + &
                                  dx3 * v(idof,2,j)             + &
                                  dx4 * v(idof,3,j)             + &
                                  dx5 * v(idof,4,j)             + &
                                  dx6 * v(idof,5,j)             + &
                                  dx7 * v(idof,6,j)             ) 
          end do
          end do

        end if

        if (rsym) then
!$doacross local( idof, isign )
!$omp parallel do private( idof, isign )
          do j = 1, ny
          do idof = 1, ndof
            if (idof .eq. 3) then
              isign = -one
            else
              isign = one
            end if
  
            ggv(1,idof,nx-2,j) = ( dx1 * v(idof,nx-5,j)         + &
                                  dx2 * v(idof,nx-4,j)          + &
                                  dx3 * v(idof,nx-3,j)          + &
                                  dx4 * v(idof,nx-2,j)          + &
                                  dx5 * v(idof,nx-1,j)          + &
                                  dx6 * v(idof,nx,j)            + &
                                  isign * dx7 * v(idof,nx-1,j)  ) 
    
            ggv(1,idof,nx-1,j) = ( dx1 * v(idof,nx-4,j)         + &
                                  dx2 * v(idof,nx-3,j)          + &
                                  dx3 * v(idof,nx-2,j)          + &
                                  dx4 * v(idof,nx-1,j)          + &
                                  dx5 * v(idof,nx,j)            + &
                                  isign * dx6 * v(idof,nx-1,j)  + &
                                  isign * dx7 * v(idof,nx-2,j)  ) 
    
            ggv(1,idof,nx,j)   = ( dx1 * v(idof,nx-3,j)         + &
                                  dx2 * v(idof,nx-2,j)          + &
                                  dx3 * v(idof,nx-1,j)          + &
                                  dx4 * v(idof,nx,j)            + &
                                  isign * dx5 * v(idof,nx-1,j)  + &
                                  isign * dx6 * v(idof,nx-2,j)  + &
                                  isign * dx7 * v(idof,nx-3,j)  ) 
          end do
          end do

        end if
        
!.... interior

!$doacross local(i,idof)
!$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 4, nx-3
            do idof = 1, ndof
              ggv(1,idof,i,j) = ( dx1 * v(idof,i-3,j)   + &
                                 dx2 * v(idof,i-2,j)   + &
                                 dx3 * v(idof,i-1,j)   + &
                                 dx4 * v(idof,i,j)     + &
                                 dx5 * v(idof,i+1,j)   + &
                                 dx6 * v(idof,i+2,j)   + &
                                 dx7 * v(idof,i+3,j)   )
            end do
          end do
        end do

!=============================================================================!
!.... compute the second derivative in y
!=============================================================================!

        if (yper) then
        
!$doacross
!$omp parallel do
          do i = 1, nx

          ggv(3,:,i,1)       = ( dy1 * v(:,i,ny-3)   + &
                                dy2 * v(:,i,ny-2)   + &
                                dy3 * v(:,i,ny-1)   + &
                                dy4 * v(:,i,1)      + &
                                dy5 * v(:,i,2)      + &
                                dy6 * v(:,i,3)      + &
                                dy7 * v(:,i,4)      ) 
  
          ggv(3,:,i,2)       = ( dy1 * v(:,i,ny-2)   + &
                                dy2 * v(:,i,ny-1)   + &
                                dy3 * v(:,i,1)      + &
                                dy4 * v(:,i,2)      + &
                                dy5 * v(:,i,3)      + &
                                dy6 * v(:,i,4)      + &
                                dy7 * v(:,i,5)      ) 
  
          ggv(3,:,i,3)       = ( dy1 * v(:,i,ny-1)   + &
                                dy2 * v(:,i,1)      + &
                                dy3 * v(:,i,2)      + &
                                dy4 * v(:,i,3)      + &
                                dy5 * v(:,i,4)      + &
                                dy6 * v(:,i,5)      + &
                                dy7 * v(:,i,6)      ) 
  
          ggv(3,:,i,ny-2)    = ( dy1 * v(:,i,ny-5)   + &
                                dy2 * v(:,i,ny-4)   + &
                                dy3 * v(:,i,ny-3)   + &
                                dy4 * v(:,i,ny-2)   + &
                                dy5 * v(:,i,ny-1)   + &
                                dy6 * v(:,i,1)      + &
                                dy7 * v(:,i,2)      ) 
  
          ggv(3,:,i,ny-1)    = ( dy1 * v(:,i,ny-4)   + &
                                dy2 * v(:,i,ny-3)   + &
                                dy3 * v(:,i,ny-2)   + &
                                dy4 * v(:,i,ny-1)   + &
                                dy5 * v(:,i,1)      + &
                                dy6 * v(:,i,2)      + &
                                dy7 * v(:,i,3)   ) 
  
          ggv(3,:,i,ny) = ggv(3,:,i,1)

          end do

        else if (bsym .or. tsym) then

          call error('grad$','Bsym and Tsym flags not implemented$')

        else
        
!$doacross
!$omp parallel do
          do i = 1, nx

          ggv(3,:,i,1)      =  ( dd1 * v(:,i,1) + &
                                 dd2 * v(:,i,2) + &
                                 dd3 * v(:,i,3) + &
                                 dd4 * v(:,i,4) + &
                                 dd5 * v(:,i,5) ) * dysinv
!!$
!!$       ggv(3,:,i,1)      =  ( dk1 * v(:,i,1) + &
!!$                              dk2 * v(:,i,2) + &
!!$                              dk3 * v(:,i,3) ) * dysinv

          ggv(3,:,i,2)      =  ( db1 * v(:,i,1) + &
                                 db2 * v(:,i,2) + &
                                 db3 * v(:,i,3) + &
                                 db4 * v(:,i,4) + &
                                 db5 * v(:,i,5) ) * dysinv
  
          ggv(3,:,i,3)       = ( da1 * v(:,i,1) + &
                                 da2 * v(:,i,2) + &
                                 da3 * v(:,i,3) + &
                                 da4 * v(:,i,4) + &
                                 da5 * v(:,i,5) ) * dysinv
  
          ggv(3,:,i,ny-2)    = ( da1 * v(:,i,ny-4) + &
                                 da2 * v(:,i,ny-3) + &
                                 da3 * v(:,i,ny-2) + &
                                 da4 * v(:,i,ny-1) + &
                                 da5 * v(:,i,ny)   ) * dysinv
  
          ggv(3,:,i,ny-1)   =  ( db1 * v(:,i,ny)   + &
                                 db2 * v(:,i,ny-1) + &
                                 db3 * v(:,i,ny-2) + &
                                 db4 * v(:,i,ny-3) + &
                                 db5 * v(:,i,ny-4) ) * dysinv
  
          ggv(3,:,i,ny)     =  ( dd1 * v(:,i,ny)   + &
                                 dd2 * v(:,i,ny-1) + &
                                 dd3 * v(:,i,ny-2) + &
                                 dd4 * v(:,i,ny-3) + &
                                 dd5 * v(:,i,ny-4) ) * dysinv

          end do
                               
        endif

!.... interior

!$doacross local(j,idof)
!$omp parallel do private(j,idof)
        do i = 1, nx
          do j = 4, ny-3
            do idof = 1, ndof
              ggv(3,idof,i,j) = ( dy1 * v(idof,i,j-3)   + &
                                  dy2 * v(idof,i,j-2)   + &
                                  dy3 * v(idof,i,j-1)   + &
                                  dy4 * v(idof,i,j  )   + &
                                  dy5 * v(idof,i,j+1)   + &
                                  dy6 * v(idof,i,j+2)   + &
                                  dy7 * v(idof,i,j+3)   ) 
            end do
          end do
        end do

!
!=============================================================================!
!.... compute the cross derivative
!=============================================================================!

        if (yper) then
        
!$doacross
!$omp parallel do
          do i = 1, nx

          ggv(2,:,i,1)       = ( gy1 * gv(1,:,i,ny-3)   + &
                                gy2 * gv(1,:,i,ny-2)    + &
                                gy3 * gv(1,:,i,ny-1)    + &
                                gy4 * gv(1,:,i,2)       + &
                                gy5 * gv(1,:,i,3)       + &
                                gy6 * gv(1,:,i,4)       ) 
  
          ggv(2,:,i,2)       = ( gy1 * gv(1,:,i,ny-2)   + &
                                gy2 * gv(1,:,i,ny-1)    + &
                                gy3 * gv(1,:,i,1)       + &
                                gy4 * gv(1,:,i,3)       + &
                                gy5 * gv(1,:,i,4)       + &
                                gy6 * gv(1,:,i,5)       ) 
  
          ggv(2,:,i,3)       = ( gy1 * gv(1,:,i,ny-1)   + &
                                gy2 * gv(1,:,i,1)       + &
                                gy3 * gv(1,:,i,2)       + &
                                gy4 * gv(1,:,i,4)       + &
                                gy5 * gv(1,:,i,5)       + &
                                gy6 * gv(1,:,i,6)       ) 
  
          ggv(2,:,i,ny-2)    = ( gy1 * gv(1,:,i,ny-5)   + &
                                gy2 * gv(1,:,i,ny-4)    + &
                                gy3 * gv(1,:,i,ny-3)    + &
                                gy4 * gv(1,:,i,ny-1)    + &
                                gy5 * gv(1,:,i,1)       + &
                                gy6 * gv(1,:,i,2)       ) 
  
          ggv(2,:,i,ny-1)    = ( gy1 * gv(1,:,i,ny-4)   + &
                                gy2 * gv(1,:,i,ny-3)    + &
                                gy3 * gv(1,:,i,ny-2)    + &
                                gy4 * gv(1,:,i,1)       + &
                                gy5 * gv(1,:,i,2)       + &
                                gy6 * gv(1,:,i,3)       ) 
  
          ggv(2,:,i,ny) = ggv(2,:,i,1)

          end do

        else if (bsym .or. tsym) then

          call error('grad$','Bsym and Tsym flags not implemented$')

        else if (carp) then  ! Implement Carpenter's boundary stencil

!$doacross
!$omp parallel do
          do i = 1, nx

          ggv(2,:,i,1)    = ( gg1 * gv(1,:,i,1) + &
                             gg2 * gv(1,:,i,2)          + &
                             gg3 * gv(1,:,i,3)          + &
                             gg4 * gv(1,:,i,4)          + &
                             gg5 * gv(1,:,i,5)          + &
                             gg6 * gv(1,:,i,6)  ) * dyinv
  
          ggv(2,:,i,2)    = ( gh1 * gv(1,:,i,1)         + &
                             gh2 * gv(1,:,i,2)          + &
                             gh3 * gv(1,:,i,3)          + &
                             gh4 * gv(1,:,i,4)          + &
                             gh5 * gv(1,:,i,5)          + &
                             gh6 * gv(1,:,i,6)  ) * dyinv

          ggv(2,:,i,3)    = ( gi1 * gv(1,:,i,1)         + &
                             gi2 * gv(1,:,i,2)          + &
                             gi3 * gv(1,:,i,3)          + &
                             gi4 * gv(1,:,i,4)          + &
                             gi5 * gv(1,:,i,5)          + &
                             gi6 * gv(1,:,i,6)  ) * dyinv

          ggv(2,:,i,4)    = ( gj1 * gv(1,:,i,1)         + &
                             gj2 * gv(1,:,i,2)          + &
                             gj3 * gv(1,:,i,3)          + &
                             gj4 * gv(1,:,i,4)          + &
                             gj5 * gv(1,:,i,5)          + &
                             gj6 * gv(1,:,i,6)  ) * dyinv

          ggv(2,:,i,ny-3) =-( gj1 * gv(1,:,i,ny)    + &
                             gj2 * gv(1,:,i,ny-1)  + &
                             gj3 * gv(1,:,i,ny-2)  + &
                             gj4 * gv(1,:,i,ny-3)  + &
                             gj5 * gv(1,:,i,ny-4)  + &
                             gj6 * gv(1,:,i,ny-5)  ) * dyinv

          ggv(2,:,i,ny-2) =-( gi1 * gv(1,:,i,ny)    + &
                             gi2 * gv(1,:,i,ny-1)  + &
                             gi3 * gv(1,:,i,ny-2)  + &
                             gi4 * gv(1,:,i,ny-3)  + &
                             gi5 * gv(1,:,i,ny-4)  + &
                             gi6 * gv(1,:,i,ny-5)  ) * dyinv

          ggv(2,:,i,ny-1) =-( gh1 * gv(1,:,i,ny)    + &
                             gh2 * gv(1,:,i,ny-1)  + &
                             gh3 * gv(1,:,i,ny-2)  + &
                             gh4 * gv(1,:,i,ny-3)  + &
                             gh5 * gv(1,:,i,ny-4)  + &
                             gh6 * gv(1,:,i,ny-5)  ) * dyinv

          ggv(2,:,i,ny)   =-( gg1 * gv(1,:,i,ny)    + &
                             gg2 * gv(1,:,i,ny-1)  + &
                             gg3 * gv(1,:,i,ny-2)  + &
                             gg4 * gv(1,:,i,ny-3)  + &
                             gg5 * gv(1,:,i,ny-4)  + &
                             gg6 * gv(1,:,i,ny-5)  ) * dyinv

          end do

        else   ! normal boundary difference

!$doacross
!$omp parallel do
          do i = 1, nx

          ggv(2,:,i,1)     = ( gc1 * gv(1,:,i,1)  + &
                              gc2 * gv(1,:,i,2)  + &
                              gc3 * gv(1,:,i,3)  + &
                              gc4 * gv(1,:,i,4)  + &
                              gc5 * gv(1,:,i,5)  ) * dyinv
!!$  
!!$       ggv(2,:,i,1)    = ( gk1 * gv(1,:,i,1)  + &
!!$                           gk2 * gv(1,:,i,2)  + &
!!$                           gk3 * gv(1,:,i,3)  ) * dyinv
 
          ggv(2,:,i,2)     = ( gb1 * gv(1,:,i,1)  + &
                              gb2 * gv(1,:,i,2)  + &
                              gb3 * gv(1,:,i,3)  + &
                              gb4 * gv(1,:,i,4)  + &
                              gb5 * gv(1,:,i,5)  ) * dyinv
  
          ggv(2,:,i,3)     = ( ga1 * gv(1,:,i,1)  + &
                              ga2 * gv(1,:,i,2)  + &
                              ga3 * gv(1,:,i,4)  + &
                              ga4 * gv(1,:,i,5)  ) * dyinv
  
          ggv(2,:,i,ny-2)  = ( ga1 * gv(1,:,i,ny-4)  + &
                              ga2 * gv(1,:,i,ny-3)  + &
                              ga3 * gv(1,:,i,ny-1)  + &
                              ga4 * gv(1,:,i,ny  )  ) * dyinv
  
          ggv(2,:,i,ny-1) = -( gb1 * gv(1,:,i,ny  )  + &
                              gb2 * gv(1,:,i,ny-1)  + &
                              gb3 * gv(1,:,i,ny-2)  + &
                              gb4 * gv(1,:,i,ny-3)  + &
                              gb5 * gv(1,:,i,ny-4)  ) * dyinv
  
          ggv(2,:,i,ny)   = -( gc1 * gv(1,:,i,ny  )  + &
                              gc2 * gv(1,:,i,ny-1)  + &
                              gc3 * gv(1,:,i,ny-2)  + &
                              gc4 * gv(1,:,i,ny-3)  + &
                              gc5 * gv(1,:,i,ny-4)  ) * dyinv

          end do
                            
        end if
        
!.... interior

!$doacross local(j,idof)
!$omp parallel do private(j,idof)
        do i = 1, nx
          do j = 4, ny-3
            do idof = 1, ndof
              ggv(2,idof,i,j) = ( gy1 * gv(1,idof,i,j-3)        + &
                                 gy2 * gv(1,idof,i,j-2) + &
                                 gy3 * gv(1,idof,i,j-1) + &
                                 gy4 * gv(1,idof,i,j+1) + &
                                 gy5 * gv(1,idof,i,j+2) + &
                                 gy6 * gv(1,idof,i,j+3) ) 
            end do
          end do
        end do

!.... implement a filter of roundoff noise

        if (.false.) then

!$doacross local(i,idof)
!$omp parallel do private(i,idof)
        do j = 1, ny
          do i = 1, nx
            do idof = 1, ndof
              if ( abs(ggv(1,idof,i,j)) .lt. eps ) ggv(1,idof,i,j) = zero
              if ( abs(ggv(2,idof,i,j)) .lt. eps ) ggv(2,idof,i,j) = zero
              if ( abs(ggv(3,idof,i,j)) .lt. eps ) ggv(3,idof,i,j) = zero
            end do
          end do
        end do
        
        end if

        return
        end

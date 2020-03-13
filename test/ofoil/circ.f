      program circ
      implicit real*8 (a-h,o-z)
      pi  = 3.1415926535897932385d+0
!      write(*,"('Enter Nx ==> ',$)")
!      read(*,*) nx
      nx = 64
      write(*,"('Enter Radius ==> ',$)")
      read(*,*) r
      x0 = 1.0-r
      y0 = 0.0
      dth = -2.0 * pi / float(nx-1)
      thmin = 2.0 * pi
!     write(10,*) nx
      do i = 1, nx
         th = thmin + float(i-1) * dth
         x  = r * cos(th)
         y  = r * sin(th)
         write(20,10) x+x0, y+y0
 10      format(2(1pe20.13,1x))
      end do
      stop
      end

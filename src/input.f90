!============================================================================!
        subroutine input
!============================================================================!
        use global
        implicit none
!============================================================================!
        write(*,"('Enter Ma ==> ',$)")
        read(*,*) Ma
        write(*,"('Enter sigma, omega ==> ',$)")
        read(*,*) sigma, omega
        write(*,"('Enter niter ==> ',$)")
        read(*,*) niter
        write(*,"('Enter restart ==> ',$)")
        read(*,*) restart
        write(*,"('Enter sweep angle ==> ',$)")
        read(*,*) lambda
        lambda = lambda * pi / 180.0
        write(*,"('Enter ispg, As, xs ==> ',$)") 
        read(*,*) ispg, As, xs
        write(*,"('Enter eps_e ==> ',$)")
        read(*,*) eps_e
!       write(*,"('Enter loctime, cflmax ==> ',$)")
!       read(*,*) loctime, cflmax

!.... setup boundary condition flags

!       left   = 0      ! symmetry
!       right  = 1      ! Reimann
!       right  = 0      ! Hard incompressible
!       bottom = 0      ! wall
!       top    = 1      ! Reimann
!       top    = 0      ! Hard incompressible

        write(*,"('Enter left, right, bottom, top BC flags ==> ',$)")
        read(*,*) left, right, bottom, top
        write(*,"('Enter xper, yper ==> ',$)")
        read(*,*) xper, yper
        write(*,"('Enter lsym, rsym, bsym, tsym ==> ',$)")
        read(*,*) lsym, rsym, bsym, tsym
        write(*,"('Enter carp ==> ',$)")
        read(*,*) carp
        write(*,*)
!       read(*,*) wakecut
!       write(*,*)

!.... get airfoil starting and stopping node for c-grid case

        if (wakecut) then
          write(*,"('Enter na, nb, circ==> ',$)")
          read(*,*) na, nb, circ
        end if
        write(*,*)

!.... turn off boundary conditions if periodic

        if (xper) then
          left  = -1
          right = -1
          lsym  = .false.
          rsym  = .false.
        end if

        if (yper) then
          top = -1
          bottom = -1
          bsym   = .false.
          tsym   = .false.
        end if

!.... check for symmetries

        if (bsym) bottom = -1
        if (tsym) top    = -1
        if (lsym) left   = -1
        if (rsym) right  = -1

        return
        end

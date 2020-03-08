!============================================================================!
        module global
!============================================================================!

!.... constants

        real, parameter :: zero    = 0.0000000000000000000d+00
        real, parameter :: pt25    = 2.5000000000000000000d-01
        real, parameter :: pt33    = 3.3333333333333333333d-01
        real, parameter :: pt5     = 5.0000000000000000000d-01
        real, parameter :: pt66    = 6.6666666666666666666d-01
        real, parameter :: one     = 1.0000000000000000000d+00
        real, parameter :: onept25 = 1.2500000000000000000d+00    
        real, parameter :: onept33 = 1.3333333333333333333d+00
        real, parameter :: onept5  = 1.5000000000000000000d+00
        real, parameter :: two     = 2.0000000000000000000d+00
        real, parameter :: three   = 3.0000000000000000000d+00
        real, parameter :: pi      = 3.1415926535897932385d+00
        real, parameter :: four    = 4.0000000000000000000d+00
        real, parameter :: five    = 5.0000000000000000000d+00
        real, parameter :: infty   = 1.0000000000000000000d+30

!.... flow variables

        real, allocatable :: phi(:,:), q(:,:,:)
        real, allocatable :: gphi(:,:,:), ggphi(:,:,:)
        real, allocatable :: A(:,:,:)
        real, allocatable :: uint(:), vint(:), cint(:)
        real, allocatable :: p(:,:), gp(:,:,:)

!$sgi distribute phi(*,block), q(*,*,block)
!$sgi distribute gphi(*,*,block), ggphi(*,*,block)
!$sgi distribute A(*,*,block)
!$sgi distribute p(*,block), g1p(*,block), g2p(*,block)

        real :: g1phil, g2phil
        real :: g11phil, g12phil, g22phil
        real :: rho, u, v, t, c, cinv, Mx, My

!.... linear system

        real, allocatable :: res(:,:), mat(:,:,:)
!$sgi distribute res(*,block), mat(*,*,block)
                
!.... wake cut treatment

        integer :: na, nb

!.... flags
#ifndef __GFORTRAN__
        integer, external :: iargc
#endif
        logical :: wakecut=.false.
        integer :: type=0
        integer, parameter :: mfile = 3 
        integer :: nfile=0, narg, iarg, ifile(mfile)
        character(80) :: arg 

!.... mesh

        integer :: i, j, k, idof
        integer :: nx, ny, nz, ndof=1, ny2
        real, allocatable :: xy(:,:,:), xi(:), eta(:)
!$sgi distribute xy(*,*,block)
        real :: dxi, deta, dxiinv, detainv, minJ

!.... time stepping

        real :: dtmax, cflmax
        real :: omega, sigma
        real, allocatable :: dtl(:,:)
!$sgi distribute dtl(*,block)
        integer :: loctime, iter, niter
        
!.... metrics

        real, allocatable :: m(:,:,:), detJ(:,:)
!$sgi distribute m(*,*,block), detJ(*,block)

!.... flow parameters

        real :: Re, Pr
        real :: Ma     = 0.1
        real :: gamma  = 1.4
        real :: gamma1 = 0.4
        real :: cv     = 716.5
        real :: cp     = 1003.1
        real :: Rgas   = 286.6
        real :: beta, lambda
        real :: circ   = 0.0
        
!.... boundary conditions and initial conditions

        integer :: restart
        integer :: left, right, bottom, top
        integer :: optx=-1, opty=-1
        logical :: xper=.false., yper=.false.
        logical :: lsym=.false., rsym=.false., bsym=.false., tsym=.false.
        logical :: carp=.false.
        integer :: bx, by

!.... sponge variables

        integer :: ispg, Ns = 4
        real    :: As, xs = 0.0, xt = 1.0
        real, allocatable :: spg(:,:), phis(:,:)
!$sgi distribute spg(*,block), phis(*,block)

!.... smoother variables

        real :: eps_e

!.... residual statistics

        real    :: norm, resl, resmax, resmin, dvmax
        integer :: im, jm, in, jn, io, jo

!.... error handler

        integer :: ier

        end module global

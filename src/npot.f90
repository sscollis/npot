!============================================================================!
        program npot
!============================================================================!
!  
!  Purpose:  Compressible potential flow solver in nonconservative form
!
!  Inputs:
!
!       Ma      => Freestream chordwise Mach number
!       sigma   => Pseudo time-step [0.1:1e5], try different value to
!                  speed convergence
!       omega   => Over-relaxation parameter (set equal to 1)
!       niter   => Number of iterations (I like 500 for each sigma)
!       restart => 0 = fresh start, 1 = restart
!       lambda  => Sweep angle (Only affects output)
!       ispg    => Sponge flag (not helpful:  set to 0)
!       As      => Sponge amplitude (set to 1)
!       xs      => Start of sponge (set to 0)
!       eps_e   => Smoother amplitude (set to 0)
!       left    => Left BC (recommend 0)
!       right   => Right BC (recommend 0)
!       bottom  => Bottom BC (recommend 0)
!       top     => Top BC (recommend 0)
!
!  Input Files:
!
!       grid.dat    => standard Plot3d grid-file
!       metric.dat  => LNS format metric file
!       restart.dat => Read if restart = 1
!
!  Output Files:
!
!       output.q  => Plot3d file for visualization
!       lns.dat     => LNS data-file format
!       *.pot       => Data for LNS potential flow boundary conditions
!
!  Notes:
!
!       There may be a better way of implementing the Riemann boudary
!       condition for the potential-flow equations.  See Shapiro: 
!       "The Dynamics and Thermodynamics of Compressible Fluid Flow, Vol. I"
!
!       The current hard (i.e. BC flag = 0) freestream boundary conditions
!       set the solution to the exact incompressible potential solution.
!       This doesn't work very well, so that the boundaries must be placed
!       very far away from the body (on the order of 1e8 nose radii).
!       I recommend an exponentially streched mesh (see confpc.f)
!
!  Author:   Scott Collis
!  
!  Date:     6-18-97
!
!  Revised:  10-12-00    Switched indices on lns.dat output file
!            10-12-00    Switched all indices
!            10-16-00    Corrected Penta solvers
!            10-18-00    Parallel version was tested and works for std
!                        boundary conditions -- need to implement
!                        periodic solvers
!            10-19-00    Changes metric data structures and gradients to 
!                        improve cache hits
!
!  Parallel environment:
!
!  OMP_NUM_THREADS #
!  OMP_SCHEDULE STATIC
!  OMP_DYNAMIC FALSE     OS cannot change number of threads 
!                        (This seems to lockup sometimes...?)
!  _DSM_MUSTRUN TRUE     Locks a thread to a CPU
!  _DSM_WAIT SPIN        Don't ever surrender a CPU
!
!  _DSM_MUSTRUN TRUE     This appears to be the important parameter!
!
!============================================================================!
        use global
        use stencil
        implicit none

!.... local variables

        real :: tmp
        real, allocatable :: mult(:), fact(:), per(:,:,:), per2(:,:)
!$sgi distribute mult(block), fact(block), per(*,*,block), per2(*,block)

        real, allocatable :: mat2(:,:,:), res2(:,:), mult2(:), fact2(:)
!$sgi distribute mat2(*,*,block), res2(*,block), mult2(block), fact2(block)

        real :: cpu, cpu1
        integer, parameter :: start = 0, end = 1
        logical :: metric_ji = .false.

        integer :: conserve=0
        real, allocatable :: res_uniform(:,:)

        integer :: blksize, numprocs, ia, ib

        real :: g1phit, g2phit, g11phit, g12phit, g22phit
        real :: gx1, gx2, gx3, gx4, gx5, gx6
        real :: gy1, gy2, gy3, gy4, gy5, gy6

        real, external :: field

        integer, external :: OMP_GET_NUM_THREADS
!============================================================================!

!.... initialize global variables (must do on SGI)

        ndof = 1; Ma = 0.1; gamma = 1.4; gamma1 = 0.4; cv = 716.5; 
        cp = 1003.1; Rgas = 286.6; Ns = 4; 
        optx=-1; opty=-1; xper=.false.; yper=.false.; 
        lsym=.false.; rsym=.false.; bsym=.false.; tsym=.false.; carp=.true.
        xs = 0.0; xt = 1.0;

!.... parse argument list

        narg = iargc()
        do iarg = 1, narg
          call getarg(iarg,arg)
          if (arg(1:1) .ne. '-') then
            nfile = nfile + 1
            if (nfile .gt. mfile) then
              write(*,*) '>> Error in argument list, too many file names'
                call exit(1)
           end if
            ifile(nfile) = iarg
          else
            select case (arg(1:3))
            case ('-wc')                ! for nonuniform bottom boundary
              wakecut = .true.
            case ('-nc') 
              carp = .false.
            case ('-ms') 
              metric_ji = .true.
            case ('-c1')
              conserve = 1
            case ('-c2')
              conserve = 2
            case ('-h')
              write(*,"('-----------------------------------------------')")
              write(*,"('Usage:  npot [options] < input.dat ')")
              write(*,"('-----------------------------------------------')")
              write(*,"('   -h:  this help')")
              write(*,"('-----------------------------------------------')")
              write(*,"('   -wc: for wake cut')")
              write(*,"('   -nc: for no carpenters')")
              write(*,"('   -ms: read metric file assuming JI ordering')")
              write(*,"('   -c1: generate conservation residual')")
              write(*,"('   -c2: use conservation residual')")
              write(*,"('-----------------------------------------------')")
              call exit(0)
            case default
              write(*,"('Argument ',i2,' ignored.')") iarg
            end select
          end if
        end do

!.... get user input

        call input

!.... read in the grid file

        open(unit=10,file='grid.dat',form='unformatted',status='old',err=100)
        read(10) nx, ny, nz
        write(*,*) "(nx,ny,nz) = ", nx, ny, nz
        allocate( xy(2,nx,ny), xi(nx), eta(ny), STAT=ier )
!$omp parallel
!$omp single
        numprocs = 1
!$      numprocs = OMP_GET_NUM_THREADS()
        blksize = nx*ny*8/numprocs/(1024)
        write(*,*) 'Running on ',numprocs,' processor(s).'
        write(*,*) 'Block size ', blksize,'KB  or  ',blksize/16,' pages.'
        if (blksize < 16) then
          write(*,*) 'WARNING:  Block size should be greater than 16 KB'
        end if
!$omp end single nowait
!$omp do private(i,j)
        do j = 1, ny
          do i = 1, nx
            xy(:,i,j) = zero
          end do
        end do
!$omp end do
!$omp end parallel 
        if (ier .ne. 0) call error('pot$', 'Insufficient Memory for grid$')
        read(10) (((xy(1,i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((xy(2,i,j), i = 1, nx), j = 1, ny), k = 1, nz), &
                 (((      tmp, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... define the boundary limits (for the boundary corners)

        bx = nx
        by = ny

!.... make the xi grid
        
        dxi = one / float(nx-1)
        dxiinv = one / dxi
        do i = 1, nx
          xi(i) = real(i-1) * dxi
        end do

!.... make the eta grid

        deta = one / float(ny-1)
        detainv = one / deta
        do j = 1, ny
          eta(j) = real(j-1) * deta
        end do

!.... allocate storage for metrics

        allocate (m(10,nx,ny), detJ(nx,ny), STAT=ier )
        if (ier .ne. 0) call error('pot$', 'Insufficient Memory for metrics$')

!$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            m(:,i,j) = zero
            detJ(i,j) = zero
          end do
        end do

!.... read in the metric file

        if (metric_ji) then
          open (unit=10,file='metric.dat',form='unformatted', &
                status='old',err=110)
          read(10) (((m(idof,i,j), j=1,ny), i=1,nx), idof = 1,10)
          close(10)
        else
          open (unit=10,file='metric.dat',form='unformatted', &
                status='old',err=110)
          read(10) (m(idof,:,:), idof = 1, 10)
          close(10)
        end if

!$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            detJ(i,j) = m(1,i,j) * m(4,i,j) - m(3,i,j) * m(2,i,j)
          end do
        end do
!       minJ = minval(detJ)
        
!.... Allocate memory for the flow data

        if (xper .or. yper) then
          allocate( per(2,nx,ny), per2(6,max(nx,ny)), STAT=ier )
!$omp parallel do private(i,j,idof)
          do j = 1, ny
            do i = 1, nx
              do idof = 1, 2
                per(idof,i,j) = zero
              end do
            end do
            do idof = 1, 6
              per2(idof,j) = zero
            end do
          end do
!$omp end parallel do
        end if
        if (ier .ne. 0) then
          call error('npot$', 'Insufficient Memory for flow variables$')
        end if

        allocate( phi(nx,ny), res(nx,ny), mat(5,nx,ny), mult(max(nx,ny)), &
                  fact(max(nx,ny)), gphi(2,nx,ny), ggphi(3,nx,ny), &
                  A(4,nx,ny), dtl(nx,ny), p(nx,ny), gp(2,nx,ny), &
                  phis(nx,ny), spg(nx,ny), STAT=ier )
        if (ier .ne. 0) then
          call error('npot$', 'Insufficient Memory for flow variables$')
        end if

!$omp parallel do private(i,j,idof)
        do j = 1, ny
          mult(j) = zero; fact(j) = zero
          do i = 1, nx
            phi(i,j) = zero; res(i,j) = zero; mat(:,i,j) = zero
            gphi(:,i,j) = zero; ggphi(:,i,j) = zero
            A(:,i,j) = zero; dtl(i,j) = one; p(i,j) = zero; gp(:,i,j) = zero
            phis(i,j) = zero; spg(i,j) = zero
          end do
        end do

150     continue

!.... Set the initial condition

        call initial()

!.... make the sponge

        if (ispg.ne.0) then
!$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              spg(i,j) = zero
              if ( xi(i) .ge. xs ) then
                spg(i,j) = As * ((xi(i)-xs)/(xt-xs))**Ns
              end if
            end do
          end do
        end if

        call get_stencil( optx, dxi, gx1, gx2, gx3, gx4, gx5, gx6 )
        call get_stencil( opty, deta, gy1, gy2, gy3, gy4, gy5, gy6 )

        if (conserve.eq.1) then
          !$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              phi(i,j) = xy(1,i,j) ! field(i,j) ! xy(1,i,j)
            end do
          end do
        else if (conserve.eq.2) then
          allocate( res_uniform(nx,ny) )
          open(20,file='uniform.dat',form='unformatted')
          read(20) res_uniform
          close(20)
        end if

!.... Begin the main loop

        do iter = 1, niter

        call timer(cpu,start,'')

!.... set the boundary conditions

        if (conserve.ne.1) call itrbc()

!.... Compute first derivatives of field in the mapped space

        call timer(cpu1,start,'')
        call grad(ndof, nx, ny, phi, gphi, dxi, deta, &
                  optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)
        call timer(cpu1, end,'Grad')
        
!.... Compute second derivatives of field
        
        call timer(cpu1, start,'')
        call grad2(ndof, nx, ny, phi, gphi, ggphi, dxi, deta, &
                   optx, opty, xper, yper, lsym, rsym, bsym, tsym, carp)
        call timer(cpu1, end,'Grad2')

!.... write out a Plot3d file

        if (.false.) then
        open(unit=20, file='grad.dat', form='unformatted')
        write(20) nx, ny, 1
        write(20) Ma, zero, zero, zero
        write(20) (( gphi(1,i,j), i = 1, nx), j = 1, ny), &
                  (( gphi(2,i,j), i = 1, nx), j = 1, ny), &
                  ((ggphi(1,i,j), i = 1, nx), j = 1, ny), &
                  ((ggphi(2,i,j), i = 1, nx), j = 1, ny), &
                  ((ggphi(3,i,j), i = 1, nx), j = 1, ny)
        close(20)
        end if

!.... transform the gradients to physical space

        call timer(cpu1, start,'')
!$omp parallel 
!$omp do private(i,j,g1phil,g2phil,g11phil,g12phil,g22phil,u,v,t,c,cinv,Mx,My)
        do j = 1, ny
          do i = 1, nx

            g1phil  = gphi(1,i,j) * m(1,i,j) + gphi(2,i,j) * m(3,i,j)
            g2phil  = gphi(1,i,j) * m(2,i,j) + gphi(2,i,j) * m(4,i,j)
            
            g11phil = ggphi(1,i,j)       * m(1,i,j)*m(1,i,j)    + &
                      two * ggphi(2,i,j) * m(1,i,j)*m(3,i,j)    + &
                      ggphi(3,i,j)       * m(3,i,j)*m(3,i,j)    + &
                      gphi(1,i,j)        * m(5,i,j)             + &
                      gphi(2,i,j)        * m(8,i,j)
            
            g12phil = ggphi(1,i,j)       * m(1,i,j)*m(2,i,j)    + &
                      ggphi(2,i,j)       * m(1,i,j)*m(4,i,j)    + &
                      ggphi(2,i,j)       * m(2,i,j)*m(3,i,j)    + &
                      ggphi(3,i,j)       * m(3,i,j)*m(4,i,j)    + &
                      gphi(1,i,j)        * m(6,i,j)             + &
                      gphi(2,i,j)        * m(9,i,j)
            
            g22phil = ggphi(1,i,j)       * m(2,i,j)*m(2,i,j)    + &
                      two * ggphi(2,i,j) * m(2,i,j)*m(4,i,j)    + &
                      ggphi(3,i,j)       * m(4,i,j)*m(4,i,j)    + &
                      gphi(1,i,j)        * m(7,i,j)             + &
                      gphi(2,i,j)        * m(10,i,j)

!.... compute other flow quantities

            u    = g1phil
            v    = g2phil
            t    = one + gamma1 * Ma**2 * pt5 * ( one - u**2 - v**2 )
            c    = sqrt(t) / Ma
            cinv = one / c
        
            Mx = u * cinv
            My = v * cinv
        
            A(1,i,j) = (one-Mx**2)*m(1,i,j)**2 + (one-My**2)*m(2,i,j)**2 - &
                       two*Mx*My*m(1,i,j)*m(2,i,j)
            A(2,i,j) = (one-Mx**2)*m(3,i,j)**2 + (one-My**2)*m(4,i,j)**2 - &
                       two*Mx*My*m(3,i,j)*m(4,i,j)
            A(3,i,j) = (one-Mx**2)*m(5,i,j) + (one-My**2)*m(7,i,j)   - &
                       two*Mx*My*m(6,i,j)
            A(4,i,j) = (one-Mx**2)*m(8,i,j) + (one-My**2)*m(10,i,j)   - &
                       two*Mx*My*m(9,i,j)

!.... compute the local time-step

!!$         call dtcfl(nx, ny, u, v, c, m(1,:,:), m(2,:,:), &
!!$                       m(3,:,:), m(4,:,:), detJ, &
!!$                    dxi, deta, cflmax, loctime, iter, dtl)

!           dtl = sqrt( maxval(detJ) ) / sqrt(detJ)

!           dtl(i,j) = one

!============================================================================!
!     Form residual
!============================================================================!

            res(i,j) = (one - Mx**2) * g11phil + (one - My**2) * g22phil - &
                       two * Mx * My * g12phil

!.... sponge

            if (ispg.ne.0) &
                 res(i,j) = res(i,j) - spg(i,j) * ( phi(i,j) - phis(i,j) )

!.... put in the time term

            res(i,j) = dtl(i,j) * sigma * omega * res(i,j)

          end do
        end do
!$omp end do
!$omp end parallel
        call timer(cpu1,end,'RHS')

!.... smoother

        if (eps_e.ne.zero) then
          call smoother( res, phi, eps_e, nx, ny )
        end if

!.... save conservative residual

        if (conserve.eq.1) then
          open(20,file='uniform.dat',form='unformatted')
          write(20) res
          close(20)
          call exit(0)
        else if (conserve.eq.2) then
          !$omp parallel do private(i,j)
          do j = 1, ny
            do i = 1, nx
              res(i,j) = res(i,j) - res_uniform(i,j)
            end do
          end do
        end if

!.... enforce BC's on the RHS

        call timer(cpu1,start,'RHSBC')
        call rhsbc()
        call timer(cpu1,end,'RHSBC')

!.... compute the norm of the residual

        call timer(cpu1,start,'Resstat')
        call resstat()
        call timer(cpu1,end,'Resstat')

!.... write out a Plot3d file

        if (.false.) then
        open(unit=20, file='res.dat', form='unformatted')
        write(20) nx, ny, 1
        write(20) Ma, zero, zero, zero
        write(20) ((    phi(i,j), i = 1, nx), j = 1, ny), &
                  ((    res(i,j), i = 1, nx), j = 1, ny), &
                  ((        one, i = 1, nx), j = 1, ny), &
                  ((one, i = 1, nx), j = 1, ny), &
                  ((one, i = 1, nx), j = 1, ny)
        close(20)
        end if

!============================================================================!
!     Form LHS and solve
!============================================================================!

        call timer(cpu1,start,'')
        call lhs1()
        call timer(cpu1,end,'LHS1')
        call timer(cpu1,start,'')
        if (xper) then
          call penta2p_blk( nx, ny, 1, mat, res, per, per2, mult, fact ) 
        else
!         call penta2bc_blk( nx, ny, 1, mat, res, mult, fact )
          call penta2bc( nx, ny, mat, res, mult, fact )
        end if
        call timer(cpu1,end,'Solve1')
          
        if (.true.) then

        call timer(cpu1,start,'')
        call lhs2()
        call timer(cpu1,end,'LHS2')
        call timer(cpu1,start,'')

        if (yper) then
          call penta1p_blk( nx, ny, 1, mat, res, per, per2, mult, fact )
        else
          if (wakecut) then
            allocate( mat2(5,na,2*ny), res2(na,2*ny), &
                      mult2(max(na,2*ny)), fact2(max(na,2*ny)) )
            do j = 1, ny
              do i = 1, na
                do k = 1, 5
                  mat2(k,i,j) = mat(5-k+1,nx-i+1,ny-j+1)
                end do
                res2(i,j) = res(nx-i+1,ny-j+1)
              end do
            end do
            mat2(:,1:na,ny+1:2*ny) = mat(:,1:na,1:ny)
            res2(1:na,ny+1:2*ny) = res(1:na,1:ny)
!           call penta1bc_blk( na, 2*ny, 1, mat2, res2, mult2, fact2 )
            call penta1bc( na, 2*ny, mat2, res2, mult2, fact2 )
            res(1:na,1:ny) = res2(1:na,ny+1:2*ny)
            do j = 1, ny
              do i = 1, na
                res(nx-i+1,ny-j+1) = res2(i,j)
              end do
            end do
            deallocate( mat2, res2, mult2, fact2 )
!           call penta1bc_blk( (nb-1)-(na+1)+1, ny, 1, mat(:,na+1:nb-1,:), &
!                              res(na+1:nb-1,:), mult, fact )
            call penta1bc( (nb-1)-(na+1)+1, ny, mat(:,na+1:nb-1,:), &
                           res(na+1:nb-1,:), mult, fact )
          else
!           call penta1bc_blk( nx, ny, 1, mat, res, mult, fact )
            call penta1bc( nx, ny, mat, res, mult, fact )
          end if
        end if
        call timer(cpu1,end,'Solve2')

        end if
        
!.... update the potential function

!$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            phi(i,j) = phi(i,j) + res(i,j)
          end do
        end do

!.... find the maximum change

        dvmax = zero
        io = 0; jo = 0

!!$omp parallel do private(i,j,resl) reduction(max: dvmax)
        do j = 1, ny
          do i = 1, nx
            resl = abs(res(i,j))
!           dvmax = max(resl,dvmax)
            if (resl .gt. dvmax) then
              io = i
              jo = j
              dvmax = resl
            end if
          end do
        end do

        call timer(cpu,end,'-')
        write(*,20) iter, norm, im, jm, resmax, io, jo, dvmax, cpu
!       write(*,20) iter, norm, im, jm, resmax, in, jn, resmin, io, jo, dvmax 
  20    format(i4,1x,1pe13.6,2(1x,'(',i3,',',i3,')',1x,1pe13.6),1x,1pe10.3)

        end do  ! iter
        
!============================================================================!
!     Output results
!============================================================================!

        call itrbc()

        allocate( q(5,nx,ny) )

!.... Compute first derivatives of field in the mapped space

        call grad(ndof, nx, ny, phi, gphi, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, carp)

!.... transform the gradients to physical space

!$omp parallel do private(i,j,g1phil,g2phil,u,v,t,rho)
        do j = 1, ny
          do i = 1, nx
            g1phil  = gphi(1,i,j) * m(1,i,j) + gphi(2,i,j) * m(3,i,j)
            g2phil  = gphi(1,i,j) * m(2,i,j) + gphi(2,i,j) * m(4,i,j)
            u   = g1phil
            v   = g2phil
            t   = one + gamma1 * Ma**2 * pt5 * (one - u**2 - v**2 )
            rho = t ** ( one / gamma1 )
            p(i,j)   = rho * t / (gamma * Ma**2)
            q(1,i,j) = rho
            q(2,i,j) = u
            q(3,i,j) = v
            q(4,i,j) = phi(i,j)
            q(5,i,j) = p(i,j)
          end do
        end do

        if (wakecut) then
          write(*,"(i4,1x,3(1pe16.8E3,1x))") na, xy(1,na,1), xy(2,na,1), p(na,1)
          write(*,"(i4,1x,3(1pe16.8E3,1x))") na, xy(1,nb,1), xy(2,nb,1), p(nb,1)
        endif

        if (wakecut) then
          open(unit=10, file='cp.dat', form='formatted')
          !write(*,*) "Output cp.dat with na, nb = ", na, nb
          do i = nb,nx
            ia = nx - i + 1
            write(10,"(8(1pe12.4E3,1x))")  xy(1,i,1),  q(4,i,1), q(3,i,1), &
                                           xy(1,ia,1), q(4,ia,1), q(3,ia,1)
!           write(*,"(8(1pe9.2,1x))") m(1,i,1), m(2,i,1), m(3,i,1), m(4,i,1)
          end do
          close(10)
        else
          open(unit=10, file='cp.dat', form='formatted')
          !write(*,*) "Output cp.dat"
          do i = 1,nx
            write(10,"(8(1pe12.4E3,1x))")  xy(1,i,1),  q(4,i,1), q(3,i,1)
          end do
          close(10)
        end if

!.... write out a Plot3d file

        open(unit=20, file='output.q', form='unformatted')
        write(20) nx, ny, 1
        write(20) Ma, zero, zero, zero
        write(20) (((q(idof,i,j), i = 1, nx), j = 1, ny), idof = 1, 5)
        close(20)

        open(unit=20, file='output.f  ', form='unformatted')
        write(20) nx, ny, 1, 5
        write(20) (((q(idof,i,j), i = 1, nx), j = 1, ny), idof = 1, 5)
        close(20)

!.... write computational grid

        open(unit=10, file='comp.xyz', form='unformatted', status='unknown')
        write(10) nx, ny, nz
        write(10) ((( xi(i), i = 1, nx), j = 1, ny), k = 1, nz), &
                  (((eta(j), i = 1, nx), j = 1, ny), k = 1, nz), &
                  ((( zero, i = 1, nx), j = 1, ny), k = 1, nz)
        close(10)

!.... write out a restart file

        open(20,file='restart.dat',form='unformatted')
        write(20) phi
        close(20)
        
!.... write out a LNS file

        deallocate( q )
        ny2 = ny                ! watch out for pg
        allocate( q(5,nx,ny2) )

!$omp parallel do private(i,j,g1phil,g2phil,u,v,t,rho)
        do j = 1, ny2
          do i = 1, nx
            g1phil = gphi(1,i,j) * m(1,i,j) + gphi(2,i,j) * m(3,i,j)
            g2phil = gphi(1,i,j) * m(2,i,j) + gphi(2,i,j) * m(4,i,j)
            u   = g1phil
            v   = g2phil
            t   = one + gamma1 * Ma**2 * pt5 * (one-u**2-v**2 )
            rho = t ** ( one / gamma1 )
            p(i,j)   = rho * t / (gamma * Ma**2)
            q(1,i,j) = rho
            q(2,i,j) = u
            q(3,i,j) = v
            q(4,i,j) = tan(lambda)
            q(5,i,j) = t
          end do
        end do

!.... write out boundary values

        open(20,file='top.pot',form='formatted')
        j = ny2
        do i = 1, nx
          write(20,*) q(1,i,j), q(2,i,j), q(3,i,j), q(4,i,j), q(5,i,j)
        end do
        close(20)

        open(20,file='right.pot',form='formatted')
        i = nx
        do j = 1, ny2
          write(20,*) q(1,i,j), q(2,i,j), q(3,i,j), q(4,i,j), q(5,i,j)
        end do
        close(20)

!.... write out the pressure gradient field

        call grad(ndof, nx, ny, p, gp, dxi, deta, optx, opty, &
                  xper, yper, lsym, rsym, bsym, tsym, carp )

        open(20,file='pg.pot',form='unformatted')
        write(20) gp(1,:,:)
        close(20)

!.... write out the LNS restart file

        open(20,file='lns.dat',form='unformatted')
        write(20) 0, zero, nx, ny2, 1, 5, zero, Ma, zero, gamma, Cv
        write(20) (((q(idof,i,j), idof = 1, 5), i = 1, nx), j = 1, ny2)
        close(20)

        if (ny2.ne.ny) then
        
!.... write the grid file

        open(unit=10, file='grid.sml', form='unformatted', status='unknown')
        write(10) nx, ny2, nz
        write(10) (((xy(1,i,j), i = 1, nx), j = 1, ny2), k = 1, nz), &
                  (((xy(2,i,j), i = 1, nx), j = 1, ny2), k = 1, nz), &
                  (((     zero, i = 1, nx), j = 1, ny2), k = 1, nz)
        close(10)

        end if
        
        stop
10      format(8(1pe13.6,1x))
100     call error('npot$','Unable to open grid file.$')
110     call error('npot$','Unable to open metric file.$')

        end

subroutine get_stencil( optx, dx, gx1, gx2, gx3, gx4, gx5, gx6)

integer :: optx
real :: dx, gx1, gx2, gx3, gx4, gx5, gx6
real :: c, w, a, b, dxinv

dxinv = one / dx

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
return
end subroutine get_stencil


subroutine timer( cpu, flag, string )

  real :: cpu
  integer :: flag
  real*4, external :: second
  character(*) :: string

  if (flag.eq.0) then
    cpu = second()
  else
    cpu = second()-cpu
#ifdef VERBOSE
    if (string(1:1).ne.'-') write(*,10) string, cpu
#endif
  end if
  return
10 format(1x,a,': ',2x,1pe10.4)
end subroutine timer

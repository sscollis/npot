!============================================================================!
        subroutine lhs2
!============================================================================!
        use global
        use stencil
        implicit none

        real :: a1, a2, a3, a4, a5
        real :: b1, b2, b3, b4, b5

        integer :: ib
!============================================================================!

!.... fourth-order first-derivative stencil

        a1 = ga1  * detainv
        a2 = ga2  * detainv
        a3 = zero * detainv
        a4 = ga3  * detainv
        a5 = ga4  * detainv

!.... fourth-order second-derivative stencil

        b1 = da1 * detainv**2
        b2 = da2 * detainv**2
        b3 = da3 * detainv**2
        b4 = da4 * detainv**2
        b5 = da5 * detainv**2

!.... second LHS interior

!$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            mat(1,i,j) = a1 * A(4,i,j) + b1 * A(2,i,j)
            mat(2,i,j) = a2 * A(4,i,j) + b2 * A(2,i,j)
            mat(3,i,j) = a3 * A(4,i,j) + b3 * A(2,i,j)
            mat(4,i,j) = a4 * A(4,i,j) + b4 * A(2,i,j)
            mat(5,i,j) = a5 * A(4,i,j) + b5 * A(2,i,j)
          end do
        end do

!.... implicit damping term

        if (eps_e.ne.zero) then
!$omp parallel do private(i,j)
          do j = 1, ny
            do i = 1, nx
              mat(1,i,j) = mat(1,i,j) - eps_e * fb1
              mat(2,i,j) = mat(2,i,j) - eps_e * fb2
              mat(3,i,j) = mat(3,i,j) - eps_e * fb3
              mat(4,i,j) = mat(4,i,j) - eps_e * fb4
              mat(5,i,j) = mat(5,i,j) - eps_e * fb5
            end do
          end do
        end if

        if (.not. yper) then

!.... bottom boundary nodes
        
        if (bsym) then
          mat(1,:,1) = zero
          mat(2,:,1) = zero
          mat(3,:,1) =      a3 * A(4,:,1) +      b3 * A(2,:,1)
          mat(4,:,1) = (a4+a2) * A(4,:,1) + (b4+b2) * A(2,:,1)
          mat(5,:,1) = (a5+a1) * A(4,:,1) + (b5+b1) * A(2,:,1)

          if (eps_e.ne.zero) then
            mat(5,:,1) = mat(5,:,1) - eps_e * fb1
            mat(4,:,1) = mat(4,:,1) - eps_e * fb2
            mat(3,:,1) = mat(3,:,1) - eps_e * fb3
            mat(4,:,1) = mat(4,:,1) - eps_e * fb4
            mat(5,:,1) = mat(5,:,1) - eps_e * fb5
          end if

          mat(1,:,2) = zero
          mat(2,:,2) =      a2 * A(4,:,2) +      b2 * A(2,:,2)
          mat(3,:,2) = (a3+a1) * A(4,:,2) + (b3+b1) * A(2,:,2)
          mat(4,:,2) =      a4 * A(4,:,2) +      b4 * A(2,:,2)
          mat(5,:,2) =      a5 * A(4,:,2) +      b5 * A(2,:,2)

          if (eps_e.ne.zero) then
            mat(3,:,2) = mat(3,:,2) - eps_e * fb1
            mat(2,:,2) = mat(2,:,2) - eps_e * fb2
            mat(3,:,2) = mat(3,:,2) - eps_e * fb3
            mat(4,:,2) = mat(4,:,2) - eps_e * fb4
            mat(5,:,2) = mat(5,:,2) - eps_e * fb5
          end if
        else if (wakecut) then

!.... left side of wake cut

          j = 1
          mat(3,1:na,j) = gk1 * detainv * A(4,1:na,j) + dk1 * detainv**2 * A(2,1:na,j) 
          mat(4,1:na,j) = gk2 * detainv * A(4,1:na,j) + dk2 * detainv**2 * A(2,1:na,j)
          mat(5,1:na,j) = gk3 * detainv * A(4,1:na,j) + dk3 * detainv**2 * A(2,1:na,j)
          mat(1,1:na,j) = zero
          mat(2,1:na,j) = zero

          j = 2
          mat(2,1:na,j) = gl1 * detainv * A(4,1:na,j) + dl1 * detainv**2 * A(2,1:na,j) 
          mat(3,1:na,j) = gl2 * detainv * A(4,1:na,j) + dl2 * detainv**2 * A(2,1:na,j)
          mat(4,1:na,j) = gl3 * detainv * A(4,1:na,j) + dl3 * detainv**2 * A(2,1:na,j)
          mat(5,1:na,j) = gl4 * detainv * A(4,1:na,j) + dl4 * detainv**2 * A(2,1:na,j)
          mat(1,1:na,j) = zero

          if (eps_e.ne.zero) then
            mat(2,1:na,j) = mat(2,1:na,j) - eps_e * (-one)
            mat(3,1:na,j) = mat(3,1:na,j) - eps_e * (two)
            mat(4,1:na,j) = mat(4,1:na,j) - eps_e * (-one)
          end if

!.... center of wake cut (body)

          j = 1
          mat(3,na+1:nb-1,j) = gc1 * detainv * A(4,na+1:nb-1,j) + dd1 * detainv**2 * A(2,na+1:nb-1,j) 
          mat(4,na+1:nb-1,j) = gc2 * detainv * A(4,na+1:nb-1,j) + dd2 * detainv**2 * A(2,na+1:nb-1,j)
          mat(5,na+1:nb-1,j) = gc3 * detainv * A(4,na+1:nb-1,j) + dd3 * detainv**2 * A(2,na+1:nb-1,j)
          mat(1,na+1:nb-1,j) = gc4 * detainv * A(4,na+1:nb-1,j) + dd4 * detainv**2 * A(2,na+1:nb-1,j)
          mat(2,na+1:nb-1,j) = gc5 * detainv * A(4,na+1:nb-1,j) + dd5 * detainv**2 * A(2,na+1:nb-1,j)

!!$       mat(3,na+1:nb-1,j) = gk1 * detainv * A(4,na+1:nb-1,j) + dk1 * detainv**2 * A(2,na+1:nb-1,j) 
!!$       mat(4,na+1:nb-1,j) = gk2 * detainv * A(4,na+1:nb-1,j) + dk2 * detainv**2 * A(2,na+1:nb-1,j)
!!$       mat(5,na+1:nb-1,j) = gk3 * detainv * A(4,na+1:nb-1,j) + dk3 * detainv**2 * A(2,na+1:nb-1,j)
!!$       mat(1,na+1:nb-1,j) = zero
!!$       mat(2,na+1:nb-1,j) = zero

          j = 2
          mat(2,na+1:nb-1,j) = gb1 * detainv * A(4,na+1:nb-1,j) + db1 * detainv**2 * A(2,na+1:nb-1,j) 
          mat(3,na+1:nb-1,j) = gb2 * detainv * A(4,na+1:nb-1,j) + db2 * detainv**2 * A(2,na+1:nb-1,j)
          mat(4,na+1:nb-1,j) = gb3 * detainv * A(4,na+1:nb-1,j) + db3 * detainv**2 * A(2,na+1:nb-1,j)
          mat(5,na+1:nb-1,j) = gb4 * detainv * A(4,na+1:nb-1,j) + db4 * detainv**2 * A(2,na+1:nb-1,j)
          mat(1,na+1:nb-1,j) = gb5 * detainv * A(4,na+1:nb-1,j) + db5 * detainv**2 * A(2,na+1:nb-1,j)

          if (eps_e.ne.zero) then
            mat(2,na+1:nb-1,j) = mat(2,na+1:nb-1,j) - eps_e * (-one)
            mat(3,na+1:nb-1,j) = mat(3,na+1:nb-1,j) - eps_e * (two)
            mat(4,na+1:nb-1,j) = mat(4,na+1:nb-1,j) - eps_e * (-one)
          end if

!.... right side of wake cut

          j = 1
          mat(3,nb:nx,j) = gk1 * detainv * A(4,nb:nx,j) + dk1 * detainv**2 * A(2,nb:nx,j) 
          mat(4,nb:nx,j) = gk2 * detainv * A(4,nb:nx,j) + dk2 * detainv**2 * A(2,nb:nx,j)
          mat(5,nb:nx,j) = gk3 * detainv * A(4,nb:nx,j) + dk3 * detainv**2 * A(2,nb:nx,j)
          mat(1,nb:nx,j) = zero
          mat(2,nb:nx,j) = zero

          j = 2
          mat(2,nb:nx,j) = gl1 * detainv * A(4,nb:nx,j) + dl1 * detainv**2 * A(2,nb:nx,j) 
          mat(3,nb:nx,j) = gl2 * detainv * A(4,nb:nx,j) + dl2 * detainv**2 * A(2,nb:nx,j)
          mat(4,nb:nx,j) = gl3 * detainv * A(4,nb:nx,j) + dl3 * detainv**2 * A(2,nb:nx,j)
          mat(5,nb:nx,j) = gl4 * detainv * A(4,nb:nx,j) + dl4 * detainv**2 * A(2,nb:nx,j)
          mat(1,nb:nx,j) = zero

          if (eps_e.ne.zero) then
            mat(2,nb:nx,j) = mat(2,nb:nx,j) - eps_e * (-one)
            mat(3,nb:nx,j) = mat(3,nb:nx,j) - eps_e * (two)
            mat(4,nb:nx,j) = mat(4,nb:nx,j) - eps_e * (-one)
          end if
        else
          j = 1
          mat(3,:,j) = gc1 * detainv * A(4,:,j) + dd1 * detainv**2 * A(2,:,j) 
          mat(4,:,j) = gc2 * detainv * A(4,:,j) + dd2 * detainv**2 * A(2,:,j)
          mat(5,:,j) = gc3 * detainv * A(4,:,j) + dd3 * detainv**2 * A(2,:,j)
          mat(1,:,j) = gc4 * detainv * A(4,:,j) + dd4 * detainv**2 * A(2,:,j)
          mat(2,:,j) = gc5 * detainv * A(4,:,j) + dd5 * detainv**2 * A(2,:,j)

          j = 2
          mat(2,:,j) = gb1 * detainv * A(4,:,j) + db1 * detainv**2 * A(2,:,j) 
          mat(3,:,j) = gb2 * detainv * A(4,:,j) + db2 * detainv**2 * A(2,:,j)
          mat(4,:,j) = gb3 * detainv * A(4,:,j) + db3 * detainv**2 * A(2,:,j)
          mat(5,:,j) = gb4 * detainv * A(4,:,j) + db4 * detainv**2 * A(2,:,j)
          mat(1,:,j) = gb5 * detainv * A(4,:,j) + db5 * detainv**2 * A(2,:,j)

          if (eps_e.ne.zero) then
            mat(2,:,j) = mat(2,:,j) - eps_e * (-one)
            mat(3,:,j) = mat(3,:,j) - eps_e * (two)
            mat(4,:,j) = mat(4,:,j) - eps_e * (-one)
          end if
        end if

!.... top boundary nodes

        if (tsym) then
          call error('lhs2$', 'tsym is not currently supported$')
        else
          mat(5,:,ny-1) = -gb5 * detainv * A(4,:,ny-1) + &
                           db5 * detainv**2 * A(2,:,ny-1) 
          mat(1,:,ny-1) = -gb4 * detainv * A(4,:,ny-1) + &
                           db4 * detainv**2 * A(2,:,ny-1)
          mat(2,:,ny-1) = -gb3 * detainv * A(4,:,ny-1) + &
                           db3 * detainv**2 * A(2,:,ny-1)
          mat(3,:,ny-1) = -gb2 * detainv * A(4,:,ny-1) + &
                           db2 * detainv**2 * A(2,:,ny-1)
          mat(4,:,ny-1) = -gb1 * detainv * A(4,:,ny-1) + &
                           db1 * detainv**2 * A(2,:,ny-1)

          if (eps_e.ne.zero) then
            mat(2,:,ny-1) = mat(2,:,ny-1) - eps_e * (-one)
            mat(3,:,ny-1) = mat(3,:,ny-1) - eps_e * (two)
            mat(4,:,ny-1) = mat(4,:,ny-1) - eps_e * (-one)
          end if

          mat(4,:,ny) = -gc5 * detainv * A(4,:,ny) + dd5 * detainv**2 * A(2,:,ny) 
          mat(5,:,ny) = -gc4 * detainv * A(4,:,ny) + dd4 * detainv**2 * A(2,:,ny)
          mat(1,:,ny) = -gc3 * detainv * A(4,:,ny) + dd3 * detainv**2 * A(2,:,ny)
          mat(2,:,ny) = -gc2 * detainv * A(4,:,ny) + dd2 * detainv**2 * A(2,:,ny)
          mat(3,:,ny) = -gc1 * detainv * A(4,:,ny) + dd1 * detainv**2 * A(2,:,ny)
        end if

        end if  ! yper

!.... put in time term

!$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            mat(1,i,j) = -sigma * dtl(i,j) * mat(1,i,j)
            mat(2,i,j) = -sigma * dtl(i,j) * mat(2,i,j)
            mat(3,i,j) =  one - sigma * dtl(i,j) * mat(3,i,j)
            mat(4,i,j) = -sigma * dtl(i,j) * mat(4,i,j)
            mat(5,i,j) = -sigma * dtl(i,j) * mat(5,i,j)
          end do
        end do
        
!.... left boundary condition

        if (left.eq.0 .or. left.eq.2) then
!$omp parallel do private(j)
          do j = 1, by
            mat(:,1,j) = zero
            mat(3,1,j) = one
          end do
        else if (left.ne.-1) then
          call error('lhs2$', 'Non supported left BC flag$')
        end if

!.... right boundary condition

        if (right.eq.0 .or. right.eq.1 .or. right.eq.2) then
!$omp parallel do private(j)
          do j = 1, by
            mat(:,nx,j) = zero
            mat(3,nx,j) = one
          end do
        else if (right.eq.-2) then

        else if (right.ne.-1) then
          call error('lhs2$', 'Non supported right BC flag$')
        end if
        
!.... top boundary condition

        if (top.eq.0) then
          mat(:,1:bx,ny) = zero
          mat(3,1:bx,ny) = -one
        else if (top.eq.1) then
          mat(:,1:bx,ny) = zero
          mat(4,1:bx,ny) = -gc5
          mat(5,1:bx,ny) = -gc4
          mat(1,1:bx,ny) = -gc3
          mat(2,1:bx,ny) = -gc2
          mat(3,1:bx,ny) = -gc1
        else if (top.ne.-1) then
          call error('lhs2$', 'Non supported top BC flag$')
        end if

!.... bottom boundary condition
        
        if (wakecut) then

!.... left portion (continuity of derivative)

          do i = 1, na
            ib = nx - i + 1
            mat(:,i,1) = zero

!!$            mat(1,i,1) =  m(4,ib,1)*gc2
!!$            mat(2,i,1) =  m(4,ib,1)*gc1
!!$            mat(3,i,1) = -m(4,i,1)*gc1
!!$            mat(4,i,1) = -m(4,i,1)*gc2
!!$            mat(5,i,1) = -m(4,i,1)*gc3

!!$            mat(1,i,1) =  m(4,ib,1)*gk2
!!$            mat(2,i,1) =  m(4,ib,1)*gk1
!!$            mat(3,i,1) = -m(4,i,1)*gk1
!!$            mat(4,i,1) = -m(4,i,1)*gk2
!!$            mat(5,i,1) = -m(4,i,1)*gk3

            mat(1,i,1) =  m(4,ib,1)
            mat(2,i,1) = -m(4,ib,1)
            mat(3,i,1) =  m(4,i,1)
            mat(4,i,1) = -m(4,i,1)
            mat(5,i,1) = zero
          end do

!.... center wall portion

          mat(:,na+1:nb-1,1) = zero
          mat(3,na+1:nb-1,1) = -gc1
          mat(4,na+1:nb-1,1) = -gc2
          mat(5,na+1:nb-1,1) = -gc3
          mat(1,na+1:nb-1,1) = -gc4
          mat(2,na+1:nb-1,1) = -gc5

!!$       mat(:,na+1:nb-1,1) = zero
!!$       mat(3,na+1:nb-1,1) = -gk1
!!$       mat(4,na+1:nb-1,1) = -gk2
!!$       mat(5,na+1:nb-1,1) = -gk3
!!$       mat(1,na+1:nb-1,1) = zero
!!$       mat(2,na+1:nb-1,1) = zero

!.... left portion (continuity or jump condition)

          mat(:,nb:nx,1) =  zero
          mat(3,nb:nx,1) = -one
          mat(2,nb:nx,1) =  one
        else          
          if (bottom.eq.0) then
            mat(:,1:bx,1) = zero
            mat(3,1:bx,1) = -gc1
            mat(4,1:bx,1) = -gc2
            mat(5,1:bx,1) = -gc3
            mat(1,1:bx,1) = -gc4
            mat(2,1:bx,1) = -gc5
!!$         mat(:,1:bx,1 = zero
!!$         mat(3,1:bx,1) = -gk1
!!$         mat(4,1:bx,1) = -gk2
!!$         mat(5,1:bx,1) = -gk3
!!$         mat(1,1:bx,1) = zero
!!$         mat(2,1:bx,1) = zero
          else if (bottom.ne.-1) then
            call error('lhs2$', 'Non supported bottom BC flag$')
          end if
        end if

        return
        end

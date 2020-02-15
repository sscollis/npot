!============================================================================!
        subroutine lhs1
!============================================================================!
        use global
        use stencil
        implicit none

        real :: a1, a2, a3, a4, a5
        real :: b1, b2, b3, b4, b5
!============================================================================!

!.... fourth-order first-derivative stencil

        a1 = ga1  * dxiinv
        a2 = ga2  * dxiinv
        a3 = zero * dxiinv
        a4 = ga3  * dxiinv
        a5 = ga4  * dxiinv

!.... fourth-order second-derivative stencil

        b1 = da1 * dxiinv**2
        b2 = da2 * dxiinv**2
        b3 = da3 * dxiinv**2
        b4 = da4 * dxiinv**2
        b5 = da5 * dxiinv**2

!.... first LHS interior

!$omp parallel do private(i,j)
        do j = 1, ny
          do i = 1, nx
            mat(1,i,j) = a1 * A(3,i,j) + b1 * A(1,i,j)
            mat(2,i,j) = a2 * A(3,i,j) + b2 * A(1,i,j)
            mat(3,i,j) = a3 * A(3,i,j) + b3 * A(1,i,j)
            mat(4,i,j) = a4 * A(3,i,j) + b4 * A(1,i,j)
            mat(5,i,j) = a5 * A(3,i,j) + b5 * A(1,i,j)
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

        if (.not. xper) then

!.... left boundary nodes

        if (lsym) then
!$omp parallel do private(i,j)
          do j = 1, ny
            i = 1
            mat(1,i,j) = zero
            mat(2,i,j) = zero
            mat(3,i,j) =      a3 * A(3,i,j) +      b3 * A(1,i,j)
            mat(4,i,j) = (a4+a2) * A(3,i,j) + (b4+b2) * A(1,i,j)
            mat(5,i,j) = (a5+a1) * A(3,i,j) + (b5+b1) * A(1,i,j)
            
            if (eps_e.ne.zero) then
              mat(5,i,j) = mat(5,i,j) - eps_e * fb1
              mat(4,i,j) = mat(4,i,j) - eps_e * fb2
              mat(3,i,j) = mat(3,i,j) - eps_e * fb3
              mat(4,i,j) = mat(4,i,j) - eps_e * fb4
              mat(5,i,j) = mat(5,i,j) - eps_e * fb5
            end if
            
            i = 2
            mat(1,i,j) = zero
            mat(2,i,j) =      a2 * A(3,i,j) +      b2 * A(1,i,j)
            mat(3,i,j) = (a3+a1) * A(3,i,j) + (b3+b1) * A(1,i,j)
            mat(4,i,j) =      a4 * A(3,i,j) +      b4 * A(1,i,j)
            mat(5,i,j) =      a5 * A(3,i,j) +      b5 * A(1,i,j)
            
            if (eps_e.ne.zero) then
              mat(3,i,j) = mat(3,i,j) - eps_e * fb1
              mat(2,i,j) = mat(2,i,j) - eps_e * fb2
              mat(3,i,j) = mat(3,i,j) - eps_e * fb3
              mat(4,i,j) = mat(4,i,j) - eps_e * fb4
              mat(5,i,j) = mat(5,i,j) - eps_e * fb5
            end if
          end do
        else
!$omp parallel do private(i,j)
          do j = 1, ny
            i = 1
            mat(3,i,j) = gc1 * dxiinv * A(3,i,j) + dd1 * dxiinv**2 * A(1,i,j) 
            mat(4,i,j) = gc2 * dxiinv * A(3,i,j) + dd2 * dxiinv**2 * A(1,i,j)
            mat(5,i,j) = gc3 * dxiinv * A(3,i,j) + dd3 * dxiinv**2 * A(1,i,j)
            mat(1,i,j) = gc4 * dxiinv * A(3,i,j) + dd4 * dxiinv**2 * A(1,i,j)
            mat(2,i,j) = gc5 * dxiinv * A(3,i,j) + dd5 * dxiinv**2 * A(1,i,j)
            
            i = 2
            mat(2,i,j) = gb1 * dxiinv * A(3,i,j) + db1 * dxiinv**2 * A(1,i,j)
            mat(3,i,j) = gb2 * dxiinv * A(3,i,j) + db2 * dxiinv**2 * A(1,i,j)
            mat(4,i,j) = gb3 * dxiinv * A(3,i,j) + db3 * dxiinv**2 * A(1,i,j)
            mat(5,i,j) = gb4 * dxiinv * A(3,i,j) + db4 * dxiinv**2 * A(1,i,j)
            mat(1,i,j) = gb5 * dxiinv * A(3,i,j) + db5 * dxiinv**2 * A(1,i,j)
            
            if (eps_e.ne.zero) then
              mat(2,i,j) = mat(2,i,j) - eps_e * (-one)
              mat(3,i,j) = mat(3,i,j) - eps_e * (two)
              mat(4,i,j) = mat(4,i,j) - eps_e * (-one)
            end if
          end do
        end if

!.... right boundary nodes

        if (rsym) then
!$omp parallel do private(i,j)
          do j = 1, ny
            i = nx-1
            mat(1,i,j) = a1 * A(3,i,j) + b1 * A(1,i,j)
            mat(2,i,j) = a2 * A(3,i,j) + b2 * A(1,i,j)
            mat(3,i,j) = (a3+a5) * A(3,i,j) + (b3+b5) * A(1,i,j)
            mat(4,i,j) = a4 * A(3,i,j) + b4 * A(1,i,j)
            mat(5,i,j) = zero

            if (eps_e.ne.zero) then
              mat(1,i,j) = mat(1,i,j) - eps_e * fb1
              mat(2,i,j) = mat(2,i,j) - eps_e * fb2
              mat(3,i,j) = mat(3,i,j) - eps_e * fb3
              mat(4,i,j) = mat(4,i,j) - eps_e * fb4
              mat(3,i,j) = mat(3,i,j) - eps_e * fb5
            end if

            i = nx
            mat(1,i,j) = (a1+a5) * A(3,i,j) + (b1+b5) * A(1,i,j)
            mat(2,i,j) = (a2+a4) * A(3,i,j) + (b2+b4) * A(1,i,j)
            mat(3,i,j) = a3 * A(3,i,j) + b3 * A(1,i,j)
            mat(4,i,j) = zero
            mat(5,i,j) = zero

            if (eps_e.ne.zero) then
              mat(1,i,j) = mat(1,i,j) - eps_e * fb1
              mat(2,i,j) = mat(2,i,j) - eps_e * fb2
              mat(3,i,j) = mat(3,i,j) - eps_e * fb3
              mat(2,i,j) = mat(2,i,j) - eps_e * fb4
              mat(1,i,j) = mat(1,i,j) - eps_e * fb5
            end if
          end do
        else
!$omp parallel do private(i,j)
          do j = 1, ny
            i = nx-1
            mat(5,i,j) = -gb5 * dxiinv * A(3,i,j) + db5 * dxiinv**2 * A(1,i,j)
            mat(1,i,j) = -gb4 * dxiinv * A(3,i,j) + db4 * dxiinv**2 * A(1,i,j)
            mat(2,i,j) = -gb3 * dxiinv * A(3,i,j) + db3 * dxiinv**2 * A(1,i,j)
            mat(3,i,j) = -gb2 * dxiinv * A(3,i,j) + db2 * dxiinv**2 * A(1,i,j)
            mat(4,i,j) = -gb1 * dxiinv * A(3,i,j) + db1 * dxiinv**2 * A(1,i,j)
            
            if (eps_e.ne.zero) then
              mat(2,i,j) = mat(2,i,j) - eps_e * (-one)
              mat(3,i,j) = mat(3,i,j) - eps_e * (two)
              mat(4,i,j) = mat(4,i,j) - eps_e * (-one)
            end if
            
            i = nx
            mat(4,i,j) = -gc5 * dxiinv * A(3,i,j) + dd5 * dxiinv**2 * A(1,i,j) 
            mat(5,i,j) = -gc4 * dxiinv * A(3,i,j) + dd4 * dxiinv**2 * A(1,i,j)
            mat(1,i,j) = -gc3 * dxiinv * A(3,i,j) + dd3 * dxiinv**2 * A(1,i,j)
            mat(2,i,j) = -gc2 * dxiinv * A(3,i,j) + dd2 * dxiinv**2 * A(1,i,j)
            mat(3,i,j) = -gc1 * dxiinv * A(3,i,j) + dd1 * dxiinv**2 * A(1,i,j)
          end do
        end if  ! rsym

        end if  ! xper

!.... wakecut

!!$        if (wakecut) then
!!$          j = 1                        !.... left side
!!$       i = na-1
!!$       mat(5,i,j) = zero
!!$       mat(1,i,j) = -gl4 * dxiinv * A(3,i,j) + dl4 * dxiinv**2 * A(1,i,j)
!!$       mat(2,i,j) = -gl3 * dxiinv * A(3,i,j) + dl3 * dxiinv**2 * A(1,i,j)
!!$       mat(3,i,j) = -gl2 * dxiinv * A(3,i,j) + dl2 * dxiinv**2 * A(1,i,j)
!!$       mat(4,i,j) = -gl1 * dxiinv * A(3,i,j) + dl1 * dxiinv**2 * A(1,i,j)
!!$
!!$       if (eps_e.ne.zero) then
!!$         mat(2,i,j) = mat(2,i,j) - eps_e * (-one)
!!$         mat(3,i,j) = mat(3,i,j) - eps_e * (two)
!!$         mat(4,i,j) = mat(4,i,j) - eps_e * (-one)
!!$       end if
!!$
!!$       i = na
!!$       mat(4,i,j) = zero
!!$       mat(5,i,j) = zero
!!$       mat(1,i,j) = -gk3 * dxiinv * A(3,i,j) + dk3 * dxiinv**2 * A(1,i,j)
!!$       mat(2,i,j) = -gk2 * dxiinv * A(3,i,j) + dk2 * dxiinv**2 * A(1,i,j)
!!$       mat(3,i,j) = -gk1 * dxiinv * A(3,i,j) + dk1 * dxiinv**2 * A(1,i,j)
!!$       
!!$
!!$       i = nb                       !.... right side
!!$       mat(3,i,j) = gk1 * dxiinv * A(3,i,j) + dk1 * dxiinv**2 * A(1,i,j) 
!!$       mat(4,i,j) = gk2 * dxiinv * A(3,i,j) + dk2 * dxiinv**2 * A(1,i,j)
!!$       mat(5,i,j) = gk3 * dxiinv * A(3,i,j) + dk3 * dxiinv**2 * A(1,i,j)
!!$       mat(1,i,j) = zero
!!$       mat(2,i,j) = zero
!!$
!!$       i = nb+1
!!$       mat(2,i,j) = gl1 * dxiinv * A(3,i,j) + dl1 * dxiinv**2 * A(1,i,j)
!!$       mat(3,i,j) = gl2 * dxiinv * A(3,i,j) + dl2 * dxiinv**2 * A(1,i,j)
!!$       mat(4,i,j) = gl3 * dxiinv * A(3,i,j) + dl3 * dxiinv**2 * A(1,i,j)
!!$       mat(5,i,j) = gl4 * dxiinv * A(3,i,j) + dl4 * dxiinv**2 * A(1,i,j)
!!$       mat(1,i,j) = zero
!!$
!!$       if (eps_e.ne.zero) then
!!$         mat(2,i,j) = mat(2,i,j) - eps_e * (-one)
!!$         mat(3,i,j) = mat(3,i,j) - eps_e * (two)
!!$         mat(4,i,j) = mat(4,i,j) - eps_e * (-one)
!!$       end if
!!$        end if

!.... sponge

        if (ispg.ne.0) then
!$omp parallel do private(i,j)
          do j = 1, ny
            do i = 1, nx
              mat(3,i,j) = mat(3,i,j) - spg(i,j)
            end do
          end do
        end if

!.... put in the time term

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

        if (left.eq.0) then
!$omp parallel do private(i,j)
          do j = 1, by
            mat(:,1,j) = zero
            mat(3,1,j) = -one
          end do
        else if (left.eq.2) then
!$omp parallel do private(i,j)
          do j = 1, by
            mat(:,1,j) =  zero
            mat(5,1,j) = -one
            mat(4,1,j) =  two
            mat(3,1,j) = -one
          end do
        else if (left.ne.-1) then
          call error('lhs1$', 'Non supported left BC flag$')
        end if
          
!.... right boundary condition

        if (right.eq.0) then
!$omp parallel do private(i,j)
          do j = 1, by
            mat(:,nx,j) = zero
            mat(3,nx,j) = -one
          end do
        else if (right.eq.1) then
!$omp parallel do private(i,j)
          do j = 1, by
            mat(:,nx,j) = zero
            mat(4,nx,j) = -gc5
            mat(5,nx,j) = -gc4
            mat(1,nx,j) = -gc3
            mat(2,nx,j) = -gc2
            mat(3,nx,j) = -gc1
          end do
        else if (right.eq.2) then
!$omp parallel do private(i,j)
          do j = 1, by
            mat(:,nx,j) =  zero
            mat(1,nx,j) = -one
            mat(2,nx,j) =  two
            mat(3,nx,j) = -one
          end do
        else if (right.eq.-2) then

        else if (right.ne.-1) then
          call error('lhs1$', 'Non supported right BC flag$')
        end if
        
!.... bottom boundary condition

        if (wakecut) then
          mat(:,1:na,1) = zero
          mat(3,1:na,1) = one

          mat(:,na+1:nb-1,1) = zero
          mat(3,na+1:nb-1,1) = one

          mat(:,nb:nx,1) = zero
          mat(3,nb:nx,1) = one
        else
          if (bottom.eq.0) then
            mat(:,1:bx,1) = zero
            mat(3,1:bx,1) = one
          else if (bottom.ne.-1) then
            call error('lhs1$', 'Non supported bottom BC flag$')
          end if
        end if

!.... top boundary condition

        if (top.eq.0 .or. top.eq.1) then
          mat(:,1:bx,ny) = zero
          mat(3,1:bx,ny) = one
        else if (top.ne.-1) then
          call error('lhs1$', 'Non supported top BC flag$')
        end if

        return
        end

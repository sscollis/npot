!=============================================================================!
        subroutine penta2bc( n, nsys, mat, rhs, mult, fact )
!=============================================================================!
!  
!  Solves a pentadiagonal system of equations without pivoting.
!  Vectorized over the second index.
!
!  This version includes the Boundary Treatment
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!
!  inout:
!    mat        : the pentadiagonal matrix to be solved
!    rhs        : the right hand side
!
!  Revised:  3-14-96
!=============================================================================!
        implicit none

        integer n, nsys
        real mat(5,n,nsys), rhs(n,nsys)
        real mult(nsys), fact(nsys)
        
!.... useful constants

        real zero, one
        parameter (zero = 0.0, one = 1.0)

!.... local variables
        
        integer nc, i, j, p, q, r, s, t, iv
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!

!$omp parallel do private(iv,i,p,r,s,t)
        do iv = 1, nsys

!.... do the first row

        i = 1

        p = i + 1
        r = i + 2
        
        fact(iv) = one / mat(3,i,iv)
        
        mult(iv) = -mat(2,p,iv) * fact(iv)
        mat(3,p,iv) = mat(3,p,iv) + mult(iv) * mat(4,i,iv)
        mat(4,p,iv) = mat(4,p,iv) + mult(iv) * mat(5,i,iv)
        mat(5,p,iv) = mat(5,p,iv) + mult(iv) * mat(1,i,iv)
        mat(1,p,iv) = mat(1,p,iv) + mult(iv) * mat(2,i,iv)
        rhs(p,iv) = rhs(p,iv) + mult(iv) * rhs(i,iv)
        
        mult(iv) = -mat(1,r,iv) * fact(iv)
        mat(2,r,iv) = mat(2,r,iv) + mult(iv) * mat(4,i,iv)
        mat(3,r,iv) = mat(3,r,iv) + mult(iv) * mat(5,i,iv)
        mat(4,r,iv) = mat(4,r,iv) + mult(iv) * mat(1,i,iv)
        mat(5,r,iv) = mat(5,r,iv) + mult(iv) * mat(2,i,iv)
        rhs(r,iv) = rhs(r,iv) + mult(iv) * rhs(i,iv)

!.... do the second row

        i = 2
        
        p = i + 1
        r = i + 2
        
        fact(iv) = one / mat(3,i,iv)

        mult(iv) = -mat(2,p,iv) * fact(iv)
        mat(3,p,iv) = mat(3,p,iv) + mult(iv) * mat(4,i,iv)
        mat(4,p,iv) = mat(4,p,iv) + mult(iv) * mat(5,i,iv)
        mat(5,p,iv) = mat(5,p,iv) + mult(iv) * mat(1,i,iv)
        rhs(p,iv) = rhs(p,iv) + mult(iv) * rhs(i,iv)
          
        mult(iv) = -mat(1,r,iv) * fact(iv)
        mat(2,r,iv) = mat(2,r,iv) + mult(iv) * mat(4,i,iv)
        mat(3,r,iv) = mat(3,r,iv) + mult(iv) * mat(5,i,iv)
        mat(4,r,iv) = mat(4,r,iv) + mult(iv) * mat(1,i,iv)
        rhs(r,iv) = rhs(r,iv) + mult(iv) * rhs(i,iv)

!.... do the interior rows

        do i = 3, n - 5

          p = i + 1
          r = i + 2
          
          fact(iv) = one / mat(3,i,iv)

          mult(iv) = -mat(2,p,iv) * fact(iv)
          mat(3,p,iv) = mat(3,p,iv) + mult(iv) * mat(4,i,iv)
          mat(4,p,iv) = mat(4,p,iv) + mult(iv) * mat(5,i,iv)
          rhs(p,iv) = rhs(p,iv) + mult(iv) * rhs(i,iv)

          mult(iv) = -mat(1,r,iv) * fact(iv)
          mat(2,r,iv) = mat(2,r,iv) + mult(iv) * mat(4,i,iv)
          mat(3,r,iv) = mat(3,r,iv) + mult(iv) * mat(5,i,iv)
          rhs(r,iv) = rhs(r,iv) + mult(iv) * rhs(i,iv)
              
        end do          ! end of loop on i

!.... do the (n-4)th row

        i = n - 4

        p = i + 1
        r = i + 2
        s = i + 3
        t = i + 4
        
        fact(iv) = one / mat(3,i,iv)
        
        mult(iv) = -mat(2,p,iv) * fact(iv)
        mat(3,p,iv) = mat(3,p,iv) + mult(iv) * mat(4,i,iv)
        mat(4,p,iv) = mat(4,p,iv) + mult(iv) * mat(5,i,iv)
        rhs(p,iv) = rhs(p,iv) + mult(iv) * rhs(i,iv)

        mult(iv) = -mat(1,r,iv) * fact(iv)
        mat(2,r,iv) = mat(2,r,iv) + mult(iv) * mat(4,i,iv)
        mat(3,r,iv) = mat(3,r,iv) + mult(iv) * mat(5,i,iv)
        rhs(r,iv) = rhs(r,iv) + mult(iv) * rhs(i,iv)

        mult(iv) = -mat(5,s,iv) * fact(iv)
        mat(1,s,iv) = mat(1,s,iv) + mult(iv) * mat(4,i,iv)
        mat(2,s,iv) = mat(2,s,iv) + mult(iv) * mat(5,i,iv)
        rhs(s,iv) = rhs(s,iv) + mult(iv) * rhs(i,iv)

        mult(iv) = -mat(4,t,iv) * fact(iv)
        mat(5,t,iv) = mat(5,t,iv) + mult(iv) * mat(4,i,iv)
        mat(1,t,iv) = mat(1,t,iv) + mult(iv) * mat(5,i,iv)
        rhs(t,iv) = rhs(t,iv) + mult(iv) * rhs(i,iv)

!.... do the (n-3)rd row

        i = n - 3

        p = i + 1
        r = i + 2
        s = i + 3
        
        fact(iv) = one / mat(3,i,iv)

        mult(iv) = -mat(2,p,iv) * fact(iv)
        mat(3,p,iv) = mat(3,p,iv) + mult(iv) * mat(4,i,iv)
        mat(4,p,iv) = mat(4,p,iv) + mult(iv) * mat(5,i,iv)
        rhs(p,iv) = rhs(p,iv) + mult(iv) * rhs(i,iv)

        mult(iv) = -mat(1,r,iv) * fact(iv)
        mat(2,r,iv) = mat(2,r,iv) + mult(iv) * mat(4,i,iv)
        mat(3,r,iv) = mat(3,r,iv) + mult(iv) * mat(5,i,iv)
        rhs(r,iv) = rhs(r,iv) + mult(iv) * rhs(i,iv)

        mult(iv) = -mat(5,s,iv) * fact(iv)
        mat(1,s,iv) = mat(1,s,iv) + mult(iv) * mat(4,i,iv)
        mat(2,s,iv) = mat(2,s,iv) + mult(iv) * mat(5,i,iv)
        rhs(s,iv) = rhs(s,iv) + mult(iv) * rhs(i,iv)

!.... do the (n-2)nd row

        i = n - 2

        p = i + 1
        r = i + 2
        
        fact(iv) = one / mat(3,i,iv)

        mult(iv) = -mat(2,p,iv) * fact(iv)
        mat(3,p,iv) = mat(3,p,iv) + mult(iv) * mat(4,i,iv)
        mat(4,p,iv) = mat(4,p,iv) + mult(iv) * mat(5,i,iv)
        rhs(p,iv) = rhs(p,iv) + mult(iv) * rhs(i,iv)

        mult(iv) = -mat(1,r,iv) * fact(iv)
        mat(2,r,iv) = mat(2,r,iv) + mult(iv) * mat(4,i,iv)
        mat(3,r,iv) = mat(3,r,iv) + mult(iv) * mat(5,i,iv)
        rhs(r,iv) = rhs(r,iv) + mult(iv) * rhs(i,iv)

!.... do the next-to-the-last row
        
        i = n - 1

        p = i + 1
        
        fact(iv) = one / mat(3,i,iv)
        mult(iv) = -mat(2,p,iv) * fact(iv)
        mat(3,p,iv) = mat(3,p,iv) + mult(iv) * mat(4,i,iv)
        rhs(p,iv) = rhs(p,iv) + mult(iv) * rhs(i,iv)

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

        i = n
        
        rhs(i,iv) = rhs(i,iv) / mat(3,i,iv)

!.... do the next-to-the-last row

        i = n - 1
        p = i + 1
        
        rhs(i,iv) = rhs(i,iv) - rhs(p,iv) * mat(4,i,iv)
        rhs(i,iv) = rhs(i,iv) / mat(3,i,iv)

!.... now do the rest of the rows in reverse order

        do i = n - 2, 3, -1
        
          p = i + 1
          r = i + 2
          
          rhs(i,iv) = rhs(i,iv) - rhs(p,iv) * mat(4,i,iv)
          rhs(i,iv) = rhs(i,iv) - rhs(r,iv) * mat(5,i,iv)
          rhs(i,iv) = rhs(i,iv) / mat(3,i,iv)

        end do
        
!.... Do the last two boundary rows

        i = 2
        p = i + 1
        r = i + 2
        s = i + 3
        
        rhs(i,iv) = rhs(i,iv) - rhs(p,iv) * mat(4,i,iv)
        rhs(i,iv) = rhs(i,iv) - rhs(r,iv) * mat(5,i,iv)
        rhs(i,iv) = rhs(i,iv) - rhs(s,iv) * mat(1,i,iv)
        rhs(i,iv) = rhs(i,iv) / mat(3,i,iv)
        
        i = 1
        p = i + 1
        r = i + 2
        s = i + 3
        t = i + 4
        
        rhs(i,iv) = rhs(i,iv) - rhs(p,iv) * mat(4,i,iv)
        rhs(i,iv) = rhs(i,iv) - rhs(r,iv) * mat(5,i,iv)
        rhs(i,iv) = rhs(i,iv) - rhs(s,iv) * mat(1,i,iv)
        rhs(i,iv) = rhs(i,iv) - rhs(t,iv) * mat(2,i,iv)
        rhs(i,iv) = rhs(i,iv) / mat(3,i,iv)

        end do
!=============================================================================!
        return
        end


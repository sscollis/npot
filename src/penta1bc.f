!=============================================================================!
        subroutine penta1bc( nsys, n, mat, rhs, mult, fact )
!=============================================================================!
!  
!  Solves a pentadiagonal system of equations without pivoting.
!  Vectorized over the first index.
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
        real mat(5,nsys,n), rhs(nsys,n)
        real mult(nsys), fact(nsys)
        
!.... useful constants

        real zero, one
        parameter (zero = 0.0, one = 1.0)

!.... local variables
        
        integer nc, i, j, l, m, p, q, r, s, t, iv
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!

!$omp parallel do private(iv,i,p,r,s,t)
        do iv = 1, nsys

!.... do the first row

        i = 1

        p = i + 1
        r = i + 2
        
        fact(iv) = one / mat(3,iv,i)
        
        mult(iv) = -mat(2,iv,p) * fact(iv)
        mat(3,iv,p) = mat(3,iv,p) + mult(iv) * mat(4,iv,i)
        mat(4,iv,p) = mat(4,iv,p) + mult(iv) * mat(5,iv,i)
        mat(5,iv,p) = mat(5,iv,p) + mult(iv) * mat(1,iv,i)
        mat(1,iv,p) = mat(1,iv,p) + mult(iv) * mat(2,iv,i)
        rhs(iv,p) = rhs(iv,p) + mult(iv) * rhs(iv,i)
        
        mult(iv) = -mat(1,iv,r) * fact(iv)
        mat(2,iv,r) = mat(2,iv,r) + mult(iv) * mat(4,iv,i)
        mat(3,iv,r) = mat(3,iv,r) + mult(iv) * mat(5,iv,i)
        mat(4,iv,r) = mat(4,iv,r) + mult(iv) * mat(1,iv,i)
        mat(5,iv,r) = mat(5,iv,r) + mult(iv) * mat(2,iv,i)
        rhs(iv,r) = rhs(iv,r) + mult(iv) * rhs(iv,i)

!.... do the second row

        i = 2
        
        p = i + 1
        r = i + 2
        
        fact(iv) = one / mat(3,iv,i)

        mult(iv) = -mat(2,iv,p) * fact(iv)
        mat(3,iv,p) = mat(3,iv,p) + mult(iv) * mat(4,iv,i)
        mat(4,iv,p) = mat(4,iv,p) + mult(iv) * mat(5,iv,i)
        mat(5,iv,p) = mat(5,iv,p) + mult(iv) * mat(1,iv,i)
        rhs(iv,p) = rhs(iv,p) + mult(iv) * rhs(iv,i)

        mult(iv) = -mat(1,iv,r) * fact(iv)
        mat(2,iv,r) = mat(2,iv,r) + mult(iv) * mat(4,iv,i)
        mat(3,iv,r) = mat(3,iv,r) + mult(iv) * mat(5,iv,i)
        mat(4,iv,r) = mat(4,iv,r) + mult(iv) * mat(1,iv,i)
        rhs(iv,r) = rhs(iv,r) + mult(iv) * rhs(iv,i)

!.... do the interior rows

        do i = 3, n - 5

          p = i + 1
          r = i + 2
          
          fact(iv) = one / mat(3,iv,i)
          
          mult(iv) = -mat(2,iv,p) * fact(iv)
          mat(3,iv,p) = mat(3,iv,p) + mult(iv) * mat(4,iv,i)
          mat(4,iv,p) = mat(4,iv,p) + mult(iv) * mat(5,iv,i)
          rhs(iv,p) = rhs(iv,p) + mult(iv) * rhs(iv,i)
          
          mult(iv) = -mat(1,iv,r) * fact(iv)
          mat(2,iv,r) = mat(2,iv,r) + mult(iv) * mat(4,iv,i)
          mat(3,iv,r) = mat(3,iv,r) + mult(iv) * mat(5,iv,i)
          rhs(iv,r) = rhs(iv,r) + mult(iv) * rhs(iv,i)
              
        end do          ! end of loop on i

!.... do the (n-4)th row

        i = n - 4

        p = i + 1
        r = i + 2
        s = i + 3
        t = i + 4
        
        fact(iv) = one / mat(3,iv,i)
        
        mult(iv) = -mat(2,iv,p) * fact(iv)
        mat(3,iv,p) = mat(3,iv,p) + mult(iv) * mat(4,iv,i)
        mat(4,iv,p) = mat(4,iv,p) + mult(iv) * mat(5,iv,i)
        rhs(iv,p) = rhs(iv,p) + mult(iv) * rhs(iv,i)
        
        mult(iv) = -mat(1,iv,r) * fact(iv)
        mat(2,iv,r) = mat(2,iv,r) + mult(iv) * mat(4,iv,i)
        mat(3,iv,r) = mat(3,iv,r) + mult(iv) * mat(5,iv,i)
        rhs(iv,r) = rhs(iv,r) + mult(iv) * rhs(iv,i)
        
        mult(iv) = -mat(5,iv,s) * fact(iv)
        mat(1,iv,s) = mat(1,iv,s) + mult(iv) * mat(4,iv,i)
        mat(2,iv,s) = mat(2,iv,s) + mult(iv) * mat(5,iv,i)
        rhs(iv,s) = rhs(iv,s) + mult(iv) * rhs(iv,i)
        
        mult(iv) = -mat(4,iv,t) * fact(iv)
        mat(5,iv,t) = mat(5,iv,t) + mult(iv) * mat(4,iv,i)
        mat(1,iv,t) = mat(1,iv,t) + mult(iv) * mat(5,iv,i)
        rhs(iv,t) = rhs(iv,t) + mult(iv) * rhs(iv,i)

!.... do the (n-3)rd row

        i = n - 3

        p = i + 1
        r = i + 2
        s = i + 3
        
        fact(iv) = one / mat(3,iv,i)

        mult(iv) = -mat(2,iv,p) * fact(iv)
        mat(3,iv,p) = mat(3,iv,p) + mult(iv) * mat(4,iv,i)
        mat(4,iv,p) = mat(4,iv,p) + mult(iv) * mat(5,iv,i)
        rhs(iv,p) = rhs(iv,p) + mult(iv) * rhs(iv,i)

        mult(iv) = -mat(1,iv,r) * fact(iv)
        mat(2,iv,r) = mat(2,iv,r) + mult(iv) * mat(4,iv,i)
        mat(3,iv,r) = mat(3,iv,r) + mult(iv) * mat(5,iv,i)
        rhs(iv,r) = rhs(iv,r) + mult(iv) * rhs(iv,i)
        
        mult(iv) = -mat(5,iv,s) * fact(iv)
        mat(1,iv,s) = mat(1,iv,s) + mult(iv) * mat(4,iv,i)
        mat(2,iv,s) = mat(2,iv,s) + mult(iv) * mat(5,iv,i)
        rhs(iv,s) = rhs(iv,s) + mult(iv) * rhs(iv,i)

!.... do the (n-2)nd row

        i = n - 2

        p = i + 1
        r = i + 2
        
        fact(iv) = one / mat(3,iv,i)

        mult(iv) = -mat(2,iv,p) * fact(iv)
        mat(3,iv,p) = mat(3,iv,p) + mult(iv) * mat(4,iv,i)
        mat(4,iv,p) = mat(4,iv,p) + mult(iv) * mat(5,iv,i)
        rhs(iv,p) = rhs(iv,p) + mult(iv) * rhs(iv,i)

        mult(iv) = -mat(1,iv,r) * fact(iv)
        mat(2,iv,r) = mat(2,iv,r) + mult(iv) * mat(4,iv,i)
        mat(3,iv,r) = mat(3,iv,r) + mult(iv) * mat(5,iv,i)
        rhs(iv,r) = rhs(iv,r) + mult(iv) * rhs(iv,i)
        
!.... do the next-to-the-last row
        
        i = n - 1

        p = i + 1
        
        fact(iv) = one / mat(3,iv,i)
        mult(iv) = -mat(2,iv,p) * fact(iv)
        mat(3,iv,p) = mat(3,iv,p) + mult(iv) * mat(4,iv,i)
        rhs(iv,p) = rhs(iv,p) + mult(iv) * rhs(iv,i)

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

        i = n
        
        rhs(iv,i) = rhs(iv,i) / mat(3,iv,i)

!.... do the next-to-the-last row

        i = n - 1
        p = i + 1
        
        rhs(iv,i) = rhs(iv,i) - rhs(iv,p) * mat(4,iv,i)
        rhs(iv,i) = rhs(iv,i) / mat(3,iv,i)

!.... now do the rest of the rows in reverse order

        do i = n - 2, 3, -1
        
          p = i + 1
          r = i + 2
          
          rhs(iv,i) = rhs(iv,i) - rhs(iv,p) * mat(4,iv,i)
          rhs(iv,i) = rhs(iv,i) - rhs(iv,r) * mat(5,iv,i)
          rhs(iv,i) = rhs(iv,i) / mat(3,iv,i)

        end do
        
!.... Do the last two boundary rows

          i = 2
          p = i + 1
          r = i + 2
          s = i + 3
          
          rhs(iv,i) = rhs(iv,i) - rhs(iv,p) * mat(4,iv,i)
          rhs(iv,i) = rhs(iv,i) - rhs(iv,r) * mat(5,iv,i)
          rhs(iv,i) = rhs(iv,i) - rhs(iv,s) * mat(1,iv,i)
          rhs(iv,i) = rhs(iv,i) / mat(3,iv,i)

          i = 1
          p = i + 1
          r = i + 2
          s = i + 3
          t = i + 4
          
          rhs(iv,i) = rhs(iv,i) - rhs(iv,p) * mat(4,iv,i)
          rhs(iv,i) = rhs(iv,i) - rhs(iv,r) * mat(5,iv,i)
          rhs(iv,i) = rhs(iv,i) - rhs(iv,s) * mat(1,iv,i)
          rhs(iv,i) = rhs(iv,i) - rhs(iv,t) * mat(2,iv,i)

          rhs(iv,i) = rhs(iv,i) / mat(3,iv,i)

        end do

!=============================================================================!
        return
        end


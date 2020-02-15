!=============================================================================!
        subroutine solve1( nsys, n, mat, rhs, mult, fact )
!=============================================================================!
!  
!  Solves a tridiagonal system of equations without pivoting.
!  Vectorized over the first index.
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!
!  inout:
!    mat        : the tridiagonal matrix to be solved
!    rhs        : the right hand side
!=============================================================================!
        implicit none

        integer, intent(in) :: n, nsys
        real, intent(inout) :: mat(3,nsys,n), rhs(nsys,n)
        real, intent(inout) :: mult(nsys), fact(nsys)
        
!.... useful constants

        real, parameter :: zero = 0.0, one = 1.0

!.... local variables
        
        integer :: i, p, iv

        !$omp parallel do private(iv,i,p)
        do iv = 1, nsys

          do i = 1, n - 1
            p  = i + 1
            fact(iv) = one / mat(2,iv,i)
            mult(iv) = -mat(1,iv,p) * fact(iv)
            mat(2,iv,p) = mat(2,iv,p) + mult(iv) * mat(3,iv,i)
            rhs(iv,p) = rhs(iv,p) + mult(iv) * rhs(iv,i)
          end do

!.... do the last row first

          i = n
          rhs(iv,i) = rhs(iv,i) / mat(2,iv,i)
          
!.... now do the rest of the rows in reverse order

          do i = n - 1, 1, -1
            p = i + 1
            rhs(iv,i) = rhs(iv,i) - rhs(iv,p) * mat(3,iv,i)
            rhs(iv,i) = rhs(iv,i) / mat(2,iv,i)
          end do

        end do

        end subroutine solve1

!=============================================================================!
        subroutine solve2( n, nsys, mat, rhs, mult, fact )
!=============================================================================!
!  
!  Solves a tridiagonal system of equations without pivoting. 
!  Vectorized over the second index. 
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!
!  inout:
!    mat        : the tridiagonal matrix to be solved
!    rhs        : the right hand side
!=============================================================================!
        implicit none

        integer, intent(in) :: n, nsys
        real, intent(inout) :: mat(3,n,nsys), rhs(n,nsys)
        real, intent(inout) :: mult(nsys), fact(nsys)
        
!.... useful constants

        real, parameter :: zero = 0.0, one = 1.0

!.... local variables
        
        integer :: i, p, iv

        !$omp parallel do private(iv,i,p)
        do iv = 1, nsys

          do i = 1, n - 1
            p  = i + 1
            fact(iv) = one / mat(2,i,iv)
            mult(iv) = -mat(1,p,iv) * fact(iv)
            mat(2,p,iv) = mat(2,p,iv) + mult(iv) * mat(3,i,iv)
            rhs(p,iv) = rhs(p,iv) + mult(iv) * rhs(i,iv)
          end do

!.... do the last row first

          i = n
          rhs(i,iv) = rhs(i,iv) / mat(2,i,iv)
          
!.... now do the rest of the rows in reverse order

          do i = n - 1, 1, -1
            p = i + 1
            rhs(i,iv) = rhs(i,iv) - rhs(p,iv) * mat(3,i,iv)
            rhs(i,iv) = rhs(i,iv) / mat(2,i,iv)
          end do

        end do

        end subroutine solve2

!=============================================================================!
        subroutine penta1bc_blk( nsys, n, nblk, mat, rhs, mult, fact )
!=============================================================================!
!  
!  Solves a block pentadiagonal system of equations without pivoting.
!  Vectorized over the first index.
!
!  This version includes the Boundary Treatment
!
!  WARNING: Since pivoting is not done, this algorithm could fail
!
!  input:
!    nsys       : number of systems to vectorize over
!    n          : number of rows in the systems
!    nblk       : dimension of the blocks in each row
!
!  inout:
!    mat        : the block pentadiagonal matrix to be solved
!    rhs        : the right hand side
!
!  Revised:  7-17-95
!=============================================================================!
        implicit none

        integer n, nblk, nsys
        real mat(5,nblk,nblk,nsys,n), rhs(nblk,nsys,n)
        real mult(nsys), fact(nsys)
        
!.... useful constants

        real zero, one
        parameter (zero = 0.0, one = 1.0)

!.... local variables
        
        integer nc, i, l, m, p, q, r, s, t, iv
!=============================================================================!
!                     F O R W A R D   E L I M I N A T I O N 
!=============================================================================!

!$omp parallel do private(i,l,m,p,q,r,s,t,iv)
        do iv = 1, nsys

!.... do the first row

        i = 1
        
!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult(iv) * mat(5,nc,q,iv,i)
                 mat(1,l,q,iv,i) = mat(1,l,q,iv,i) + mult(iv) * mat(1,nc,q,iv,i)
                 mat(2,l,q,iv,i) = mat(2,l,q,iv,i) + mult(iv) * mat(2,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
           end do             
           
        end do                  ! end of loop on nc
        
!.... for all columns
        
        do nc = 1, nblk
           
           p = i + 1
           r = i + 2
           
           fact(iv) = one / mat(3,nc,nc,iv,i)
           
           do m = 1, nblk
              mult(iv) = -mat(2,m,nc,iv,p) * fact(iv)
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult(iv) * mat(5,nc,q,iv,i)
                 mat(5,m,q,iv,p) = mat(5,m,q,iv,p) + mult(iv) * mat(1,nc,q,iv,i)
                 mat(1,m,q,iv,p) = mat(1,m,q,iv,p) + mult(iv) * mat(2,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)
              
              mult(iv) = -mat(1,m,nc,iv,r) * fact(iv)
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult(iv) * mat(5,nc,q,iv,i)
                 mat(4,m,q,iv,r) = mat(4,m,q,iv,r) + mult(iv) * mat(1,nc,q,iv,i)
                 mat(5,m,q,iv,r) = mat(5,m,q,iv,r) + mult(iv) * mat(2,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult(iv) * rhs(nc,iv,i)
              
           end do

        end do                  ! end of loop on nc

!.... do the second row

        i = 2
        
!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult(iv) * mat(5,nc,q,iv,i)
                 mat(1,l,q,iv,i) = mat(1,l,q,iv,i) + mult(iv) * mat(1,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
           end do             
           
        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk
           
           p = i + 1
           r = i + 2
           
           fact(iv) = one / mat(3,nc,nc,iv,i)
           
           do m = 1, nblk
              mult(iv) = -mat(2,m,nc,iv,p) * fact(iv)
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult(iv) * mat(5,nc,q,iv,i)
                 mat(5,m,q,iv,p) = mat(5,m,q,iv,p) + mult(iv) * mat(1,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)
              
              mult(iv) = -mat(1,m,nc,iv,r) * fact(iv)
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult(iv) * mat(5,nc,q,iv,i)
                 mat(4,m,q,iv,r) = mat(4,m,q,iv,r) + mult(iv) * mat(1,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult(iv) * rhs(nc,iv,i)
              
           end do
           
        end do                  ! end of loop on nc

!.... do the interior rows

        do i = 3, n - 5

!.... do the first four columns

           do nc = 1, nblk - 1
              
              do m = nc, nblk - 1
                 l = m + 1
                 mult(iv) = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
                 do q = nc + 1, nblk
                    mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
                 end do
                 do q = 1, nblk
                    mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult(iv) * mat(4,nc,q,iv,i)
                    mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult(iv) * mat(5,nc,q,iv,i)
                 end do
                 rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
              end do          
              
           end do               ! end of loop on nc
           
!.... for all columns
          
           do nc = 1, nblk
              
              p = i + 1
              r = i + 2
              
              fact(iv) = one / mat(3,nc,nc,iv,i)
              
              do m = 1, nblk
                 mult(iv) = -mat(2,m,nc,iv,p) * fact(iv)
                 do q = nc + 1, nblk
                    mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
                 end do
                 do q = 1, nblk
                    mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult(iv) * mat(4,nc,q,iv,i)
                    mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult(iv) * mat(5,nc,q,iv,i)
                 end do
                 rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)
                 
                 mult(iv) = -mat(1,m,nc,iv,r) * fact(iv)
                 do q = nc + 1, nblk
                    mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult(iv) * mat(3,nc,q,iv,i)
                 end do
                 do q = 1, nblk
                    mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult(iv) * mat(4,nc,q,iv,i)
                    mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult(iv) * mat(5,nc,q,iv,i)
                 end do
                 rhs(m,iv,r) = rhs(m,iv,r) + mult(iv) * rhs(nc,iv,i)
                 
              end do
              
           end do               ! end of loop on nc
           
        end do                  ! end of loop on i

!.... do the (n-4)th row

        i = n - 4

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
           end do             
           
        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk
           
           p = i + 1
           r = i + 2
           s = i + 3
           t = i + 4
           
           fact(iv) = one / mat(3,nc,nc,iv,i)
           
           do m = 1, nblk
              mult(iv) = -mat(2,m,nc,iv,p) * fact(iv)
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)
              
              mult(iv) = -mat(1,m,nc,iv,r) * fact(iv)
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult(iv) * rhs(nc,iv,i)
              
              mult(iv) = -mat(5,m,nc,iv,s) * fact(iv)
              do q = nc + 1, nblk
                 mat(5,m,q,iv,s) = mat(5,m,q,iv,s) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(1,m,q,iv,s) = mat(1,m,q,iv,s) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(2,m,q,iv,s) = mat(2,m,q,iv,s) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,s) = rhs(m,iv,s) + mult(iv) * rhs(nc,iv,i)
              
              mult(iv) = -mat(4,m,nc,iv,t) * fact(iv)
              do q = nc + 1, nblk
                 mat(4,m,q,iv,t) = mat(4,m,q,iv,t) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(5,m,q,iv,t) = mat(5,m,q,iv,t) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(1,m,q,iv,t) = mat(1,m,q,iv,t) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,t) = rhs(m,iv,t) + mult(iv) * rhs(nc,iv,i)
              
           end do
           
        end do                  ! end of loop on nc

!.... do the (n-3)rd row
        
        i = n - 3

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
           end do             
           
        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk
           
           p = i + 1
           r = i + 2
           s = i + 3
           
           fact(iv) = one / mat(3,nc,nc,iv,i)
           
           do m = 1, nblk
              mult(iv) = -mat(2,m,nc,iv,p) * fact(iv)
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)
              
              mult(iv) = -mat(1,m,nc,iv,r) * fact(iv)
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult(iv) * rhs(nc,iv,i)
              
              mult(iv) = -mat(5,m,nc,iv,s) * fact(iv)
              do q = nc + 1, nblk
                 mat(5,m,q,iv,s) = mat(5,m,q,iv,s) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(1,m,q,iv,s) = mat(1,m,q,iv,s) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(2,m,q,iv,s) = mat(2,m,q,iv,s) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,s) = rhs(m,iv,s) + mult(iv) * rhs(nc,iv,i)
              
           end do
           
        end do                  ! end of loop on nc

!.... do the (n-2)nd row

        i = n - 2

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(5,l,q,iv,i) = mat(5,l,q,iv,i) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
           end do             
           
        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk
           
           p = i + 1
           r = i + 2
           
           fact(iv) = one / mat(3,nc,nc,iv,i)
           
           do m = 1, nblk
              mult(iv) = -mat(2,m,nc,iv,p) * fact(iv)
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(4,m,q,iv,p) = mat(4,m,q,iv,p) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)
              
              mult(iv) = -mat(1,m,nc,iv,r) * fact(iv)
              do q = nc + 1, nblk
                 mat(1,m,q,iv,r) = mat(1,m,q,iv,r) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(2,m,q,iv,r) = mat(2,m,q,iv,r) + mult(iv) * mat(4,nc,q,iv,i)
                 mat(3,m,q,iv,r) = mat(3,m,q,iv,r) + mult(iv) * mat(5,nc,q,iv,i)
              end do
              rhs(m,iv,r) = rhs(m,iv,r) + mult(iv) * rhs(nc,iv,i)
              
           end do
           
        end do                  ! end of loop on nc

!.... do the next-to-the-last row
        
        i = n - 1

!.... do the first four columns

        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(4,l,q,iv,i) = mat(4,l,q,iv,i) + mult(iv) * mat(4,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
           end do             
           
        end do                  ! end of loop on nc

!.... for all columns
          
        do nc = 1, nblk
           
           p = i + 1
           
           fact(iv) = one / mat(3,nc,nc,iv,i)
           
           do m = 1, nblk
              mult(iv) = -mat(2,m,nc,iv,p) * fact(iv)
              do q = nc + 1, nblk
                 mat(2,m,q,iv,p) = mat(2,m,q,iv,p) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              do q = 1, nblk
                 mat(3,m,q,iv,p) = mat(3,m,q,iv,p) + mult(iv) * mat(4,nc,q,iv,i)
              end do
              rhs(m,iv,p) = rhs(m,iv,p) + mult(iv) * rhs(nc,iv,i)
           end do
           
        end do                  ! end of loop on nc

!.... do the last row

        i = n
        
        do nc = 1, nblk - 1
           
           do m = nc, nblk - 1
              l = m + 1
              mult(iv) = -mat(3,l,nc,iv,i) / mat(3,nc,nc,iv,i)
              do q = nc + 1, nblk
                 mat(3,l,q,iv,i) = mat(3,l,q,iv,i) + mult(iv) * mat(3,nc,q,iv,i)
              end do
              rhs(l,iv,i) = rhs(l,iv,i) + mult(iv) * rhs(nc,iv,i)
           end do             
           
        end do                  ! end of loop on nc

!=============================================================================!
!                   B A C K W A R D   S U B S T I T U T I O N 
!=============================================================================!

!.... do the last row first

        i = n
        
        do m = nblk, 1, -1
           do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
           end do
           rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
        end do
        
!.... do the next-to-the-last row

        i = n - 1
        p = i + 1
        
        do m = nblk, 1, -1
           do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
           end do
           do q = 1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(4,m,q,iv,i)
           end do
           rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
        end do

!.... now do the rest of the rows in reverse order

        do i = n - 2, 3, -1
           
           p = i + 1
           r = i + 2
           
           do m = nblk, 1, -1
              do q = m+1, nblk
                 rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
              end do
              do q = 1, nblk
                 rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(4,m,q,iv,i)
                 rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,r) * mat(5,m,q,iv,i)
              end do
              rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
           end do
           
        end do
        
!.... Do the last two boundary rows

        i = 2
        p = i + 1
        r = i + 2
        s = i + 3
        
        do m = nblk, 1, -1
           do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
           end do
           do q = 1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(4,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,r) * mat(5,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,s) * mat(1,m,q,iv,i)
           end do
           rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
        end do

        i = 1
        p = i + 1
        r = i + 2
        s = i + 3
        t = i + 4
        
        do m = nblk, 1, -1
           do q = m+1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,i) * mat(3,m,q,iv,i)
           end do
           do q = 1, nblk
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,p) * mat(4,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,r) * mat(5,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,s) * mat(1,m,q,iv,i)
              rhs(m,iv,i) = rhs(m,iv,i) - rhs(q,iv,t) * mat(2,m,q,iv,i)
           end do
           rhs(m,iv,i) = rhs(m,iv,i) / mat(3,m,m,iv,i)
        end do
        
        end do
        
        return
        end


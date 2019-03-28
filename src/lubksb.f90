subroutine lubksb(a,n,indx,b)
!-------------------------------------------------------------------------
!  a(n,n):: input, real values, the LU decomposition of a matrix.
!       n:: input, integer, the dimenstion of the matrix.
! indx(n):: input,  integer values, vector that records the row 
!           permutation effected by the partial pivoting.
!    b(n):: output, real values, the solution vector X for 
!                   linear equations A*X=B.
!-------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.18.
!-------------------------------------------------------------------------
! Dependence:: No.--------------------------------------------------------
!-------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.39 IN Press et al.
! -------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n, indx(n)
    real   (kind=8), intent(in):: a(n,n)
    real   (kind=8), intent(inout):: b(n)
   ! Local variables.
    integer(kind=4):: i, ii, j, ll
    real   (kind=8):: summ

    ii = 0

    do i=1, n

        ll = indx(i)

        summ = b(ll)

        b(ll) = b(i)

        if (ii .ne. 0) then

            do j=ii, i-1

                summ = summ - a(i,j) * b(j)

            end do

        else if (summ .ne. 0.0) then
            
            ii = i

        end if

        b(i) = summ

    end do

    do i=n, 1, -1

        summ = b(i)

        do j=i+1, n

            summ = summ - a(i,j) * b(j)

        end do

        b(i) = summ / a(i,i)

    end do

    return

end subroutine lubksb

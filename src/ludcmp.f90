subroutine ludcmp(a,n,indx,d,flag)
!-------------------------------------------------------------------------
!This routine is used in combination with lubksb to solve 
!linear equations or invert a matrix.
!-------------------------------------------------------------------------
!  a(n,n):: input, real values, a matrix to be decomposed.
!       n:: input, integer, the dimension of the matrix.
! indx(n):: output, integer values, vector that records the row 
!           permutation effected by the partial pivoting.
!       d:: output, integer, output as 1 or -1 depending on whether 
!           the number of row interchanges was even or odd.
!    flag:: output, integer, error message, 0=success, 1=singular matrix.
!-------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.20.
!-------------------------------------------------------------------------
! Dependence:: No.--------------------------------------------------------
!-------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.38 IN Press et al.
! ------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: n
    integer(kind=4), intent(out):: indx(n), flag
    real   (kind=8), intent(inout):: a(n,n)
    real   (kind=8), intent(out):: d
    ! Local variables.
    integer(kind=4):: i, j, k, imax
    real   (kind=8):: aamax, dum, summ, vv(n)
   
    indx = 0

    flag = 0

    d = 1.0

    do i=1, n

        aamax = 0.0

        do j=1, n

            if (abs(a(i,j)) .gt. aamax)  aamax = abs(a(i,j))

        end do

        if (aamax .eq. 0.0) then 
 
            flag = 1

            return

        end if

        vv(i) = 1.0/aamax

    end do

    do j=1, n

        do i=1, j-1

            summ = a(i,j)

            do k=1, i-1

                summ = summ - a(i,k) * a(k,j)

            end do

            a(i,j) = summ

        end do

        aamax = 0.0

        do i=j, n

            summ = a(i,j)

            do k=1, j-1

                summ = summ - a(i,k) * a(k,j)

            end do

            a(i,j) = summ

            dum = vv(i) * abs(summ)

            if (dum .ge. aamax) then
       
                imax = i
          
                aamax = dum

            end if

        end do

        if (j .ne. imax) then

            do k=1, n

                dum = a(imax,k)

                a(imax,k) = a(j,k)

                a(j,k) = dum
  
            end do

            d = -d

            vv(imax) = vv(j)

        end if

        indx(j) = imax

        if (a(j,j) .eq. 0.0) a(j,j) = tiny(0.0D+00)

        if (j .ne. n) then

            dum = 1.0 / a(j,j)

            do i=j+1, n

                a(i,j) = a(i,j) * dum

            end do

        end if

    end do

    return

end subroutine ludcmp    

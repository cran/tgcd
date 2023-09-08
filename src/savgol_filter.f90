subroutine savgol_filter(nl,nr,ld,m,n1,y,flag)
!-----------------------------------------------------------------------------------
! This routine is used to perform the Savitzky-Golay algorithm.
!-----------------------------------------------------------------------------------
!    nl:: input, integer, the number of leftward data points used.
!    nr:: input, integer, the number of rightward data points used.
!    ld:: input, integer, the order of the derivative desired.
!     m:: input, integer, the order of the smoothing polynomial.
!    n1:: input, integer, the number of data points.
! y(n1):: input/output, real values, the data to be smoothed.
!  flag:: output, integer, error message, 0=success, 1=failure.
!-----------------------------------------------------------------------------------
! Author: Peng Jun, 2023.09.07.
!-----------------------------------------------------------------------------------
! Dependence:: subroutine savgol.
! -----------------------------------------------------------------------------------

    implicit none
    integer, intent(in):: nl, nr, ld, m, n1
    real(kind(1.0d0)), intent(inout):: y(n1)
    integer, intent(out):: flag
    ! Local variables.
    integer:: i, j, xl(nl+nr+1)
    real(kind(1.0d0)):: y0(n1), coef(nl+nr+1)

    xl(1) = 0

    y0 = y

    do i=1, nl

        xl(i+1) = -i

    end do

    do i=1, nr

        xl(1+nl+i) = nr-i+1

    end do

    call savgol(nl,nr,ld,m,coef,flag)

    if (flag/=0) return

    do i=1, n1-nr

        y(i) = 0.0

        do j=1, nl+nr+1

            if (i+xl(j) .gt. 0) then

                y(i) = y(i) + coef(j)*y0(i+xl(j))

            end if

        end do

    end do

    if (ld==0) then

        y(1:nl) = y0(1:nl)

        y(n1-nr+1:n1) = y0(n1-nr+1:n1)

    else 

        y(1:nl) = y(nl+1)

        y(n1-nr+1:n1) = y(n1-nr)
 
    end if

    return

end subroutine savgol_filter

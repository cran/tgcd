subroutine savgol(nl,nr,ld,m,coef,flag)
!-----------------------------------------------------------------------------------
! This routine is used to calculate a set of Savitzky-Golay filter coefficients.
!-----------------------------------------------------------------------------------
!            nl:: input, integer, the number of leftward data points used.
!            nr:: input, integer, the number of rightward data points used.
!            ld:: input, integer, the order of the derivative desired.
!             m:: input, integer, the order of the smoothing polynomial.
! coef(nl+nr+1):: output, real values, calculated coefficents in wrap-around order.
!          flag:: output, integer, error message, 0=success, 1=failure.
!-----------------------------------------------------------------------------------
! Author: Peng Jun, 2019.03.20.
!-----------------------------------------------------------------------------------
! Dependence:: subroutine ludcmp;
!              subroutine lubksb.
!-----------------------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.646 IN Press et al.
! -----------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nl, nr, ld, m
    real   (kind=8), intent(inout):: coef(nl+nr+1)
    integer(kind=4), intent(out):: flag
    ! Local variables.
    integer(kind=4):: imj, ipj, k, kk, mm, indx(m+1)

    real   (kind=8):: d, fac, summ, a(m+1,m+1), b(m+1)

    flag = 0

    if (nl < 0 .or. nr < 0 .or. ld > m .or. nl+nr < m) then

        flag = 1

        return

    end if

    do ipj=0, 2*m

        summ = 0.0

        if (ipj .eq. 0) summ = 1.0

        do k=1, nr

            summ = summ + (float(k))**ipj

        end do

        do k=1, nl

            summ = summ + (float(-k))**ipj

        end do

        mm = min(ipj, 2*m-ipj)

        do imj=-mm, mm, 2

            a(1+(ipj+imj)/2,1+(ipj-imj)/2) = summ

        end do

    end do

    call ludcmp(a,m+1,indx,d,flag)

    if (flag .ne. 0) return

    b = 0.0

    b(ld+1) = 1.0

    call lubksb(a,m+1,indx,b)

    coef = 0.0

    do k=-nl, nr

        summ = b(1)

        fac = 1.0

        do mm=1, m

            fac = fac * k

            summ = summ + b(mm+1) * fac

        end do

        kk = mod(nl+nr+1-k, nl+nr+1) + 1

        coef(kk) = summ

    end do

    return

end subroutine savgol

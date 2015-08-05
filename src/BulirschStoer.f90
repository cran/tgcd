subroutine BulirschStoer(nt, vect, y0, tol, &
                         kmax, ff, ae, hr, vecy)
!--------------------------------------------------------------------------------
! Subroutine BulirschStoer is used for solving an ordinary differential equation 
! of the form dn/dT=-n*ff*exp(ae/k/T)/hr using the Bulirsch-Stoer algorithm.
!--------------------------------------------------------------------------------
!        nt, input:: integer, the number of data points.
!  vect(nt), input:: real values, temperature values.
!        y0, input:: real value, starting value of n.
!       tol, input:: real value, relative tolerance.
!      kmax, input:: integer, maximal number of steps.
!        ff, input:: real value, the frequency factor.
!        ae, input:: real value, the activation energy.
!        hr, input:: real value, the heating rate.
! vecy(nt), output:: real values, y values at various temperatures.
!--------------------------------------------------------------------------------
! Author:: Peng Jun, 2015.06.16.
!--------------------------------------------------------------------------------
! Dependence:: subroutine midpoint.
!--------------------------------------------------------------------------------
! NOTE:: THIS SUBROUTINE IS MODIFIED USING THE ORIGINAL R FUNCTION 
!        bulirsch_stoer IN R PACKAGE "pracma" WRITTEN BY Hans Werner Borchers.
!--------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nt, kmax
    real   (kind=8), intent(in):: vect(nt), y0, &
                                  tol, ff, ae, hr
    real   (kind=8), intent(out):: vecy(nt)

    real   (kind=8):: mtol, value1
    integer(kind=4):: i

    mtol = tol/real(nt-1)

    vecy(1) = y0

    do i=2, nt
        call midpoint(vect(i-1), vect(i), vecy(i-1),& 
                      tol, kmax, ff, ae, hr, value1)
        vecy(i) = value1
    end do
    return
end subroutine BulirschStoer
!
subroutine midpoint(t0, tn, y0, tol, kmax,&
                    ff, ae, hr, value1)
!--------------------------------------------------------
! Subroutine midpoint performs a modified midpoint method
! to computes values of y values from t0 to tn using the
! Richardson extrapolation method in each step.
!---------------------------------------------------------
!      t0, input:: real value, starting value of t.
!      tn, input:: real value, terminating value of t.
!      y0, input:: real value, initial y value at t0.
!     tol, input:: real value, relative tolerance.
!    kmax, input:: integer, maximal number of steps.
!      ff, input:: real value, the frequency factor.
!      ae, input:: real value, the activation energy.
!      hr, input:: real value, the heating rate.
! value1, output:: real value, the resulting estimate.
!----------------------------------------------------------
! Author:: Peng Jun, 2015.06.16.
!----------------------------------------------------------
! Dependence:: subroutine midp.
!----------------------------------------------------------
! NOTE:: THIS SUBROUTINE IS MODIFIED USING THE ORIGINAL 
!        R FUNCTION midpoint IN R PACKAGE "pracma" WRITTEN 
!        BY Hans Werner Borchers.
!-----------------------------------------------------------
    implicit none
    real   (kind=8), intent(in):: t0, tn, y0, &
                                  tol, ff, ae, hr
    integer(kind=4), intent(in):: kmax
    real   (kind=8), intent(out):: value1

    integer(kind=4):: nstep, i, j
    real   (kind=8):: rvec(kmax), value, rold, & 
                      cc, dr, errr

    value1 = 2.220446e-16

    nstep = 2
    call midp(t0, tn, y0, nstep,& 
              ff, ae, hr, value)   
    rvec(1) = value

    rold = rvec(1)

    do i= 2, kmax
        nstep = 2*nstep 
        call midp(t0, tn, y0, nstep,&
                  ff, ae, hr, value) 
        rvec(i) = value

        do j=i-1, 1, -1
            cc = 4.0**(i-j)
            rvec(j) = (cc*rvec(j+1)-rvec(j))/(cc-1.0)
        end do

        dr = rvec(1) - rold
        errr = abs(dr)    
        rold = rvec(1)

        if (rvec(1)>=0.0 .and. &
            rvec(1)<=y0)  value1 = rvec(1)

        if (errr<tol) exit
    end do
    
    return
end subroutine midpoint
!
subroutine midp(t0, tn, y0, nstep,&
                ff, ae, hr, value) 
!-------------------------------------------------------
! midp is a subroutine for performing extrapolating.
!-------------------------------------------------------
! Author:: Peng Jun, 2015.06.16.
!-------------------------------------------------------
! Dependence:: inner function func.
!-------------------------------------------------------
! NOTE:: THIS SUBROUTINE IS MODIFIED USING THE ORIGINAL 
!        R FUNCTION midp IN R PACKAGE "pracma" WRITTEN 
!        BY Hans Werner Borchers.
!-------------------------------------------------------
    implicit none
    real   (kind=8), intent(in):: t0, tn, y0, &
                                  ff, ae, hr
    integer(kind=4), intent(in):: nstep
    real   (kind=8), intent(out):: value
    !
    real   (kind=8):: h, t, y1, y2, yy
    real   (kind=8), parameter:: kbz = 8.617385e-05
    integer(kind=4):: i

    h = (tn-t0)/nstep
    t = t0
    y1 = y0
    y2 = y1 + h*func(t, y1)

    do i=1, nstep-1
        t = t + h
        yy = y1 + 2.0*h*func(t, y2)
        y1 = y2
        y2 = yy
    end do
    value = 0.5*(y1+y2+h*func(t, yy))
    return
    contains
        function func(T, n)
            implicit none
            real   (kind=8):: T, n
            real   (kind=8):: func
            !
            func = -n*ff*exp(-ae/kbz/T)/hr 
            return
        end function func
end subroutine midp

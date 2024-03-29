subroutine calcMaty_pdf2(nd,n2,pars,xd,maty,bg)
!---------------------------------------------------------------------
! Subroutine calcMaty_pdf2() is used for calculating
! the signal matrix of optimized TL glow peaks according
! to first-order kinetic using the logistic asymmetric equation.
!---------------------------------------------------------------------
!                  nd:: input, integer, number of data points.
!                  n2:: input, integer, number of pars (<=39+3).
!            pars(n2):: input, real values, pars.
!              xd(nd):: input, real values, observations X.
! maty(nd,(n2-3)/3+1):: output, real values, calculated signals.
!                  bg:: input, integer, subtract background or not,
!                       0=no subtraction, 1=subtraction.
!---------------------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!---------------------------------------------------------------------
! Dependence:: NO.
!---------------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer, intent(in):: nd, n2, bg
    real(kind(1.0d0)), intent(in):: pars(n2), xd(nd)
    real(kind(1.0d0)), intent(out):: maty(nd,(n2-3)/3+1)    
    ! Local variables.
    real(kind(1.0d0)), parameter:: kbz=8.617385e-5
    real(kind(1.0d0)):: xx(39+3), maxi, engy, maxt, a2v, xa(nd), expv1(nd),&
                        ba, bb, bc
    integer:: i, n0
    !
    n0 = n2 - 3
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    !
    do i=1, n0/3
        maxi = xx(i)
        engy = xx(i+n0/3)
        maxt = xx(i+2*n0/3)
        !
        a2v = sqrt(1.189*maxt**4*kbz**2/(engy**2+4.0*engy*maxt*kbz))
        !
        xa = (xd-maxt)/a2v
        !
        expv1 = exp(-(xa+0.38542))
        !
        maty(:,i) = 5.2973*maxi*(1.0+expv1)**(-2.4702)*expv1
    end do
    !
    if (bg==0) then
        !
        maty(:,n0/3+1) = 0.0
    else if (bg==1) then
        ba = xx(n0+1)
        bb = xx(n0+2)
        bc = xx(n0+3)
        !
        maty(:,n0/3+1) = ba + bb * exp(xd/bc)
    end if
    !
    return
    !
end subroutine calcMaty_pdf2

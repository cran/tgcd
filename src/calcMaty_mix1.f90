subroutine calcMaty_mix1(nd,n2,pars,xd,maty,bg)
!------------------------------------------------------------------
! Subroutine calcMaty_mix1() is used for calculating
! the signal matrix of optimized TL glow peaks
! according to the mix-order function (type 1).
!------------------------------------------------------------------
!                  nd:: input, integer, number of data points.
!                  n2:: input, integer, number of pars (<=52+3).
!            pars(n2):: input, real values, pars.
!              xd(nd):: input, real values, observations X.
! maty(nd,(n2-3)/4+1):: output, real values, calculated signals.
!                  bg:: input, integer, subtract background or not,
!                       0=no subtraction, 1=subtraction.
!-----------------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!-----------------------------------------------------------------
! Dependence:: subroutine calcAm.
!-----------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer, intent(in):: nd, n2, bg
    real(kind(1.0d0)), intent(in):: pars(n2),xd(nd)
    real(kind(1.0d0)), intent(out):: maty(nd,(n2-3)/4+1)                 
    ! Local variables.
    real(kind(1.0d0)), parameter:: kbz=8.617385e-5
    real(kind(1.0d0)):: xx(52+3), maxi, engy, maxt, Am, Rm, fmin,&
                        alpha, xa, xb(nd), expv1, expv2(nd), expv3(nd),&
                        ba, bb, bc
    integer:: i, j, n0
    !
    n0 = n2 - 3
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    do i=1, n0/4
        maxi = xx(i)
        engy = xx(i+n0/4)
        maxt = xx(i+2*n0/4)
        alpha = xx(i+3*n0/4)
        !
        xa = 2.0*kbz*maxt/engy
        xb = 2.0*kbz*xd/engy
        !
        call calcAm(1.0D-03,1.0D+01,alpha,maxt,engy,Am,fmin)
        Rm = (Am+alpha)/(Am-alpha)
        !
        expv1 = exp((1.0-xa)/Rm)
        expv2 = exp(engy/kbz/xd*(xd-maxt)/maxt)
        expv3 = exp(xd**2/maxt**2/Rm*expv2*(1.0-xb))
        ! 
        ! Check expv3 for underflow or overflow.
        do j=1, nd
            if (expv3(j)>=huge(0.0D+00)) expv3(j)=huge(0.0D+00)
            if (expv3(j)<=tiny(0.0D+00)) expv3(j)=tiny(0.0D+00)
        end do
        !
        maty(:,i) = maxi*(expv1-alpha)*((expv1-alpha)/expv1)*&
                    (expv2/(expv3-alpha))*(expv3/(expv3-alpha))
        !
    end do
    !
    if (bg==0) then
        !
        maty(:,n0/4+1) = 0.0
    else if (bg==1) then
        ba = xx(n0+1)
        bb = xx(n0+2)
        bc = xx(n0+3)
        !
        maty(:,n0/4+1) = ba + bb * exp(xd/bc)
    end if
    !
    return
    !
end subroutine calcMaty_mix1

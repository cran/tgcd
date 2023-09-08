subroutine calcMaty_mix2(nd,n2,pars,xd,maty,bg)
!-------------------------------------------------------------------
! Subroutine calcMaty_mix2() is used for calculating
! the signal matrix of optimized TL glow peaks
! according to the mix-order function (type 2).
!-------------------------------------------------------------------
!                  nd:: input, integer, number of data points.
!                  n2:: input, integer, number of pars (<=52+3).
!            pars(n2):: input, real values, pars.
!              xd(nd):: input, real values, observations X.
! maty(nd,(n2-3)/4+1):: output, real values, calculated signals.
!                  bg:: input, integer, subtract background or not,
!                       0=no subtraction, 1=subtraction.
!-------------------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!-------------------------------------------------------------------
! Dependence:: subroutine calcFct.
!-------------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer, intent(in):: nd, n2, bg
    real(kind(1.0d0)), intent(in):: pars(n2),xd(nd)
    real(kind(1.0d0)), intent(out):: maty(nd,(n2-3)/4+1)                 
    ! Local variables.
    real(kind(1.0d0)), parameter:: kbz=8.617385e-5
    real(kind(1.0d0)):: xx(52+3), maxi, engy, maxt, alpha,&
                        xa, xb(nd), Rm, fxa, fxb(nd),&
                        expv1(nd), expv2(nd), expv3(nd),&
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
        xa = engy/kbz/maxt
        xb = engy/kbz/xd
        Rm = (1.0-alpha)*(1.0+0.2922*alpha-0.2783*alpha**2)
        !
        call calcFct(xa, fxa)
        !
        do j=1, nd
            call calcFct(xb(j), fxb(j))
        end do
        !
        expv1 = exp(xa-xb)
        expv2 = exp(Rm*xa*(xd/maxt*expv1*fxb-fxa))
        !
        ! Check expv2 for underflow or overflow.
        do j=1, nd
            if (expv2(j)>=huge(0.0D+00)) expv2(j)=huge(0.0D+00)
            if (expv2(j)<=tiny(0.0D+00)) expv2(j)=tiny(0.0D+00)
        end do
        !
        expv3 = (1.0+Rm)*expv2-(1.0-Rm)
        !
        maty(:,i) = maxi*(4.0*Rm**2)*(expv1/expv3)*(expv2/expv3)
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
end subroutine calcMaty_mix2

subroutine calcMaty_gnr2(nd,n2,pars,xd,maty,bg)
!-----------------------------------------------------------------
! Subroutine calcMaty_gnr2() is used for calculating
! the signal matrix of optimized TL glow peaks 
! according to the general-order equation (type 2).
!-----------------------------------------------------------------
!                  nd:: input, integer, number of data points.
!                  n2:: input, integer, number of pars (<=52+4).
!            pars(n2):: input, real values, pars.
!              xd(nd):: input, real values, observations X.
! maty(nd,(n2-4)/4+1):: output, real values, calculated signals.
!                  bg:: input, integer, subtract background or not,
!                       0=no subtraction, 1=subtraction.
!-----------------------------------------------------------------
! Author:: Peng Jun, 2019.03.24.
!-----------------------------------------------------------------
! Dependence:: subroutine calcFct.
!-----------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4), intent(in):: nd, n2, bg
    real   (kind=8), intent(in):: pars(n2), xd(nd)
    real   (kind=8), intent(out):: maty(nd,(n2-4)/4+1)                
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(52+4), maxi, engy, maxt, &
                      bv, xa, xb(nd), expv(nd), fxa, fxb(nd),&
                      ba, bb, bc, bd
    integer(kind=4):: i, j, n0
    !
    n0 = n2 -4
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    do i=1, n0/4
        maxi = xx(i)
        engy = xx(i+n0/4)
        maxt = xx(i+2*n0/4)
        bv = xx(i+3*n0/4)
        !
        xa = engy/kbz/maxt
        xb = engy/kbz/xd
        expv = exp(xa-xb)
        !
        call calcFct(xa, fxa)
        !
        do j=1, nd
            call calcFct(xb(j), fxb(j))
        end do
        !
        maty(:,i) = maxi*expv*(1.0+(bv-1.0)/bv*xa*&
                   (xd/maxt*expv*fxb-fxa))**(-bv/(bv-1.0))
    end do
    !
    if (bg==0) then
        !
        maty(:,n0/4+1) = 0.0
    else if (bg==1) then
        ba = xx(n0+1)
        bb = xx(n0+2)
        bc = xx(n0+3)
        bd = xx(n0+4)
        !
        maty(:,n0/4+1) = ba-bb/(1.0+exp(bc*(xd-bd)))
    end if
    !
    return
    !
end subroutine calcMaty_gnr2

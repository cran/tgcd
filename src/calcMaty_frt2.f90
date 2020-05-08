subroutine calcMaty_frt2(nd,n2,pars,xd,maty,bg)
!------------------------------------------------------------------
! Subroutine calcMaty_frt2() is used for calculating
! the signal matrix of optimized TL glow peaks according 
! to first-order kinetic (type 2).
!------------------------------------------------------------------
!                  nd:: input, integer, number of data points.
!                  n2:: input, integer, number of pars (<=39+3).
!            pars(n2):: input, real values, pars.
!              xd(nd):: input, real values, observations X.
! maty(nd,(n2-3)/3+1):: output, real values, calculated signals.
!                  bg:: input, integer, subtract background or not,
!                       0=no subtraction, 1=subtraction.
!------------------------------------------------------------------
! Author:: Peng Jun, 2020.05.08.
!------------------------------------------------------------------
! Dependence:: NO.
!------------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4), intent(in):: nd, n2, bg
    real   (kind=8), intent(in):: pars(n2), xd(nd)
    real   (kind=8), intent(out):: maty(nd,(n2-3)/3+1)    
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(39+3), maxi, engy, maxt, xa, xb(nd), xv(nd),&
                      ba, bb, bc
    integer(kind=4):: i, n0
    !
    n0 = n2 - 3
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    do i=1, n0/3
        maxi = xx(i)
        engy = xx(i+n0/3)
        maxt = xx(i+2*n0/3)
        !
        xa = 2.0*kbz*maxt/engy
        xb = 2.0*kbz*xd/engy
        xv = engy/kbz/xd*((xd-maxt)/maxt)
        !
        maty(:,i) = maxi*exp(1.0+xv-((xd/maxt)**2)*exp(xv)*(1.0-xb)-xa)
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
end subroutine calcMaty_frt2

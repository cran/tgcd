subroutine calcMaty_frt1(nd,n2,pars,xd,maty,bg)
!--------------------------------------------------------------------
! Subroutine calcMaty_frt1() is used for calculating 
! the signal  matrix of optimized TL glow peaks according 
! to first-order kinetic (type 1).
!--------------------------------------------------------------------
!                  nd:: input, integer, number of data points.
!                  n2:: input, integer, number of pars (<=39+4).
!            pars(n2):: input, real values, pars.
!              xd(nd):: input, real values, observations X.
! maty(nd,(n2-4)/3+1):: output, real values, calculated signals.
!                  bg:: input, integer, subtract background or not,
!                       0=no subtraction, 1=subtraction.
!--------------------------------------------------------------------
! Author:: Peng Jun, 2019.03.24.
!--------------------------------------------------------------------
! Dependence:: NO.
!--------------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4), intent(in):: nd, n2, bg
    real   (kind=8), intent(in):: pars(n2), xd(nd)
    real   (kind=8), intent(out):: maty(nd,(n2-4)/3+1)    
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5, a0=0.267773734, &
                 a1=8.6347608925, a2=18.059016973, a3=8.5733287401, &
                 b0=3.9584969228, b1=21.0996530827, b2=25.6329561486, &
                 b3=9.5733223454
    real   (kind=8):: xx(39+4), maxi, engy, maxt, xa, fxa, xb(nd), fxb(nd),&
                      ba, bb, bc, bd
    integer(kind=4):: i, n0
    !
    n0 = n2 - 4
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
        xa = engy/kbz/maxt
        xb = engy/kbz/xd
        fxa = 1.0-(a0+a1*xa+a2*xa**2+a3*xa**3+xa**4)/&
                  (b0+b1*xa+b2*xa**2+b3*xa**3+xa**4)
        fxb = 1.0-(a0+a1*xb+a2*xb**2+a3*xb**3+xb**4)/&
                  (b0+b1*xb+b2*xb**2+b3*xb**3+xb**4)
        !
        maty(:,i) = maxi*exp(xa-xb)*exp(xa*(fxa-xd/maxt*fxb*exp(xa-xb)))
    end do
    !
    if (bg==0) then
        !
        maty(:,n0/3+1) = 0.0
    else if (bg==1) then
        ba = xx(n0+1)
        bb = xx(n0+2)
        bc = xx(n0+3)
        bd = xx(n0+4)
        !
        maty(:,n0/3+1) = ba-bb/(1.0+exp(bc*(xd-bd)))
    end if
    !
    return
    !
end subroutine calcMaty_frt1

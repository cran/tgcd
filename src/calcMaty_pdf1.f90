subroutine calcMaty_pdf1(nd,n2,pars,xd,maty,bg)
!--------------------------------------------------------------------
! Subroutine calcMaty_pdf1() is used for calculating the signal 
! matrix of optimized TL glow peaks according to first-order 
! kinetic using the weibull distribution.
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
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(39+4), maxi, engy, maxt, bv, xa(nd),&
                      ba, bb, bc, bd
    integer(kind=4):: i, n0
    !
    n0 = n2 - 4
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    do i=1, n0/3
        maxi = xx(i)
        engy = xx(i+n0/3)
        maxt = xx(i+2*n0/3)
        !
        bv = sqrt((maxt*kbz)**2*242.036*maxt**2/((engy+maxt*kbz)**2-7.0*(maxt*kbz)**2))
        !
        xa = (xd-maxt)/bv
        !
        maty(:,i) = 2.713*maxi*(xa+0.996)**15*exp(-(xa+0.996)**16)
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
end subroutine calcMaty_pdf1

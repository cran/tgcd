subroutine calcMaty_lw(nd,n2,pars,xd,maty,bg)
!-------------------------------------------------------------------------
! Subroutine calcMaty_lw() is used for calculating
! the signal matrix of optimized TL glow peaks
! according to the Lambert’s W function.
!-------------------------------------------------------------------------
!                  nd:: input, integer, number of data points.
!                  n2:: input, integer, number of pars (<=52+3).
!            pars(n2):: input, real values, pars.
!              xd(nd):: input, real values, observations X.
! maty(nd,(n2-3)/4+1):: output, real values, calculated signals.
!                  bg:: input, integer, subtract background or not,
!                       0=no subtraction, 1=subtraction.
!-------------------------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!-------------------------------------------------------------------------
! Dependence:: subroutine calcei; 
!              subroutine wrightOmega.
!-------------------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer, intent(in):: nd, n2, bg
    real(kind(1.0d0)), intent(in):: pars(n2),xd(nd)
    real(kind(1.0d0)), intent(out):: maty(nd,(n2-3)/4+1)                 
    ! Local variables.
    real(kind(1.0d0)), parameter:: kbz=8.617385e-5
    real(kind(1.0d0)):: xx(52+3), maxi, engy, maxt, &
                        rv, eiv, eivi, Feivi, ftev(nd), ftem, &
                        z1v(nd), z1m, wz1v(nd), wz1m, xi, wv,&
                        ba, bb, bc
    integer:: i, j, n0
    !
    n0 = n2 - 3
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    xi = minval(xd)
    !
    do i=1, n0/4
        maxi = xx(i)
        engy = xx(i+n0/4)
        maxt = xx(i+2*n0/4)
        rv = xx(i+3*n0/4)
        !
        call calcei(-engy/kbz/xi, eivi, 1)
        Feivi = xi*exp(-engy/kbz/xi) + engy/kbz*eivi
        !
        ! Calculate part1: vector wz1v.
        do j=1, nd
            call calcei(-engy/kbz/xd(j), eiv, 1)
            ftev(j) = (xd(j)*exp(-engy/kbz/xd(j)) + engy/kbz*eiv) - Feivi
        end do
        !
        z1v = rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/&
              kbz/maxt**2/(1.0-1.05*rv**1.26)*ftev
        !
        do j=1, nd
            call wrightOmega(z1v(j), wv)
            wz1v(j) = wv
        end do 
        !
        !
        ! Calculate part2: scalar wz1m.
        call calcei(-engy/kbz/maxt, eiv, 1)
        ftem = (maxt*exp(-engy/kbz/maxt) + engy/kbz*eiv) - Feivi
        !
        z1m = rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/&
              kbz/maxt**2/(1.0-1.05*rv**1.26)*ftem
        !
        call wrightOmega(z1m, wv)
        wz1m = wv
        !
        ! Calculate residuals.
        maty(:,i) = maxi*(wz1m+wz1m**2)/(wz1v+wz1v**2)*&
                    exp(-engy/kbz*(1.0/xd-1.0/maxt))
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
end subroutine calcMaty_lw

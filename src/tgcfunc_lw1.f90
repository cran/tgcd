subroutine tgcfunc_lw1(nd,n2,pars,fvec,iflag,&
                       xd,yd,lower,upper,bg)
!-------------------------------------------------------
! Subroutine tgcfunc_lw1() is used for calculating
! the residual vector of a TL glow curve 
! according to the Lambert’s W function (both branches).
!-------------------------------------------------------
!        nd:: input, integer, number of data points.
!        n2:: input, integer, number of pars (<=52+3).
!  pars(n2):: input, real values, pars.
!  fvec(nd):: output, real values, residuals.
!     iflag:: input, integer, working variable.
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observations Y.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real vlaues, upper bounds.
!        bg:: input, integer, subtract background or not,
!             0=no subtraction, 1=subtraction.
!--------------------------------------------------------
! Author:: Peng Jun, 2020.06.06.
!--------------------------------------------------------
! Dependence:: subroutine calcei; 
!              subroutine wrightOmega;
!              subroutine lambertW.
!--------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag, bg
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)               
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(52+3), maxi, engy, maxt, &
                      rv, eiv, eivi, Feivi, ftev(nd), ftem, &
                      z1v(nd), z1m, wz1v(nd), wz1m, xi, wv, &
                      ba, bb, bc, tinyv
    integer(kind=4):: i, j, n0, ner
    !
    tinyv = tiny(0.0D+00)
    !
    ! Bound constraints.
    do i=1, n2
        if (pars(i)<lower(i))  then
            pars(i) = lower(i)
        else if (pars(i)>upper(i)) then
            pars(i) = upper(i)
        end if
    end do
    !
    n0 = n2 - 3
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    xi = minval(xd)
    !
    if (bg==0) then
        !
        fvec = 0.0
    else if (bg==1) then
        ba = xx(n0+1)
        bb = xx(n0+2)
        bc = xx(n0+3)
        !
        fvec = ba + bb * exp(xd/bc)
    end if
    !
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
        if (rv<1.0) then
            z1v = rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/&
                  kbz/maxt**2/(1.0-1.05*rv**1.26)*ftev
        else 
            z1v = abs(rv/(1.0-rv)) + log(abs((1.0-rv)/rv)) + engy*exp(engy/kbz/maxt)/&
                  kbz/maxt**2/(2.963-3.24*rv**(-0.74))*ftev
        end if
        !
        do j=1, nd
            if (rv<1.0) then
                call wrightOmega(z1v(j), wv)
                wz1v(j) = wv   
            else 
                ! Approach use log function.
                if (exp(-z1v(j)) < tinyv) then
                    wz1v(j) = -z1v(j) - log(z1v(j))
                else 
                    ! Calculate use lambertW.
                    call lambertW(-exp(-z1v(j)), wv, ner)
                    wz1v(j) = wv 
                    if (ner/=0) then
                        fvec = 1.0D+30
                        return
                    end if
                    !
                end if
                !
            end if 
            !
        end do 
        !
        !
        ! Calculate part2: scalar wz1m.
        call calcei(-engy/kbz/maxt, eiv, 1)
        ftem = (maxt*exp(-engy/kbz/maxt) + engy/kbz*eiv) - Feivi
        !
        if (rv<1.0) then
            z1m = rv/(1.0-rv) - log((1.0-rv)/rv) + engy*exp(engy/kbz/maxt)/&
                  kbz/maxt**2/(1.0-1.05*rv**1.26)*ftem
        else 
            z1m = abs(rv/(1.0-rv)) + log(abs((1.0-rv)/rv)) + engy*exp(engy/kbz/maxt)/&
                  kbz/maxt**2/(2.963-3.24*rv**(-0.74))*ftem  
        end if
        !
        if (rv<1.0) then
            call wrightOmega(z1m, wv)
            wz1m = wv
        else 
            ! Approach use log function.
            if (exp(-z1m) < tinyv) then
                wz1m = -z1m - log(z1m)
            else 
                ! Calculate use lambertW.
                call lambertW(-exp(-z1m), wv, ner)
                wz1m = wv
                if (ner/=0) then
                    fvec = 1.0D+30
                    return
                end if
                ! 
            end if 
            !
        end if 
        !
        ! Calculate residuals.
        fvec = fvec + maxi*(wz1m+wz1m**2)/(wz1v+wz1v**2)*&
               exp(-engy/kbz*(1.0/xd-1.0/maxt))
        !
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_lw1

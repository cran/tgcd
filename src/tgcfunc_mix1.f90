subroutine tgcfunc_mix1(nd,n2,pars,fvec,iflag,&
                        xd,yd,lower,upper,bg)
!---------------------------------------------------------
! Subroutine tgcfunc_mix1() is used for calculating
! the residual vector of a TL glow curve according 
! to the mix-order equation (type 1).
!---------------------------------------------------------
!        nd:: input, integer, number of data points.
!        n2:: input, integer, number of pars (<=52+4).
!  pars(n2):: input, real values, pars.
!  fvec(nd):: output, real values, residuals.
!     iflag:: input, integer, working variable.
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observations Y.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real vlaues, upper bounds.
!        bg:: input, integer, subtract background or not,
!             0=no subtraction, 1=subtraction.
!----------------------------------------------------------
! Author:: Peng Jun, 2019.03.24.
!----------------------------------------------------------
! Dependence:: subroutine calcAm.
!----------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag, bg
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)               
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(52+4), maxi, engy, maxt, Am, Rm, fmin,&
                      alpha, xa, xb(nd), expv1, expv2(nd), expv3(nd),&
                      ba, bb, bc, bd
    integer(kind=4):: i, j, n0
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
    n0 = n2 - 4
    !
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    if (bg==0) then
        !
        fvec = 0.0
    else if (bg==1) then
        ba = xx(n0+1)
        bb = xx(n0+2)
        bc = xx(n0+3)
        bd = xx(n0+4)
        !
        fvec = ba-bb/(1.0+exp(bc*(xd-bd)))
    end if
    !
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
        call calcAm(1.0D-02,1.0D+01,alpha,maxt,engy,Am,fmin)
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
        fvec = fvec + maxi*(expv1-alpha)*((expv1-alpha)/expv1)*&
                      (expv2/(expv3-alpha))*(expv3/(expv3-alpha))
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_mix1

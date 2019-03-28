subroutine tgcfunc_mix2(nd,n2,pars,fvec,iflag,&
                        xd,yd,lower,upper,bg)
!-----------------------------------------------------------
! Subroutine tgcfunc_mix2() is used for calculating
! the residual vector of a TL glow curve according 
! to the mix-order equation (type 2).
!-----------------------------------------------------------
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
!-----------------------------------------------------------
! Author:: Peng Jun, 2019.03.24.
!-----------------------------------------------------------
! Dependence:: subroutine calcFct.
!-----------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag, bg
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)               
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(52+4), maxi, engy, maxt, alpha,&
                      xa, xb(nd), Rm, fxa, fxb(nd),&
                      expv1(nd), expv2(nd), expv3(nd),&
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
    n0 = n2 -4
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
        fvec = fvec + maxi*(4.0*Rm**2)*(expv1/expv3)*(expv2/expv3)
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_mix2

subroutine tgcfunc_mix3(nd,n2,pars,fvec,iflag,&
                        xd,yd,lower,upper,bg)
!---------------------------------------------------------
! Subroutine tgcfunc_mix3() is used for calculating
! the residual vector of a TL glow curve according 
! to the mix-order equation (type 3).
!---------------------------------------------------------
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
!---------------------------------------------------------
! Author:: Peng Jun, 2020.05.08.
!---------------------------------------------------------
! Dependence:: subroutine calcLv.
!---------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag, bg
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)               
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(52+3), maxi, engy, maxt, alpha,&
                      xa(nd), xb(nd), expv1(nd), expv2(nd), lv, fmin,&
                      ba, bb, bc
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
    n0 = n2 - 3
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
        !
        fvec = ba + bb * exp(xd/bc)
    end if
    !
    do i=1, n0/4
        maxi = xx(i)
        engy = xx(i+n0/4)
        maxt = xx(i+2*n0/4)
        alpha = xx(i+3*n0/4)
        !
        xa = 2.0*kbz*xd/engy
        xb = (xd-maxt)/maxt
        !
        expv1 = exp(2.0/xa*xb)
        !
        call calcLv(1.0D+00,2.0D+00,alpha,maxt,engy,lv,fmin)  
        !
        expv2 = exp(xd**2/maxt**2*(2.0/lv-1.0)*expv1*(1.0-xa))
        !
        ! Check expv2 for underflow or overflow.
        do j=1, nd
            if (expv2(j)>=huge(0.0D+00)) expv2(j)=huge(0.0D+00)
            if (expv2(j)<=tiny(0.0D+00)) expv2(j)=tiny(0.0D+00)
        end do
        !
        fvec = fvec + alpha*maxi*(2.0-lv)**2/(lv-1.0)*&
                     (expv1/(expv2-alpha))*(expv2/(expv2-alpha))
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_mix3

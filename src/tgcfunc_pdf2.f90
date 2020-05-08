subroutine tgcfunc_pdf2(nd,n2,pars,fvec,iflag,&
                        xd,yd,lower,upper,bg)
!-----------------------------------------------------------------------
! Subroutine tgcfunc_pdf2() is used for calculating
! the residual vector of a optimized TL glow curve according
! to second-order kinetic using the logistic asymmetric equation.
!-----------------------------------------------------------------------
!        nd:: input, integer, number of data points.
!        n2:: input, integer, number of pars (<=39+3).
!  pars(n2):: input, real values, pars.
!  fvec(nd):: output, real values, residuals.
!     iflag:: input, integer, working variable.
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observations Y.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real vlaues, upper bounds.
!        bg:: input, integer, subtract background or not,
!             0=no subtraction, 1=subtraction.
!-----------------------------------------------------------------------
! Author:: Peng Jun, 2020.05.08.
!-----------------------------------------------------------------------
! Dependence:: NO.
!-----------------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag, bg
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)                 
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(39+3), maxi, engy, maxt, a2v, xa(nd), expv1(nd),&
                      ba, bb, bc                  
    integer(kind=4):: i, n0
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
    do i=1, n0/3
        maxi = xx(i)
        engy = xx(i+n0/3)
        maxt = xx(i+2*n0/3)
        !
        a2v = sqrt(1.189*maxt**4*kbz**2/(engy**2+4.0*engy*maxt*kbz))
        !
        xa = (xd-maxt)/a2v
        !
        expv1 = exp(-(xa+0.38542))
        !
        fvec = fvec + 5.2973*maxi*(1.0+expv1)**(-2.4702)*expv1
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_pdf2

subroutine tgcfunc_gnr1(nd,n2,pars,fvec,iflag,&
                       xd,yd,lower,upper,bg)
!-----------------------------------------------------------
! Subroutine tgcfunc_gnr1() is used for calculating
! the residual vector of a TL glow curve according 
! to the general-order equation (type 1).
!-----------------------------------------------------------
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
!------------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!------------------------------------------------------------
! Dependence:: NO.
!------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer:: nd, n2, iflag, bg
    real(kind(1.0d0)):: pars(n2), lower(n2), upper(n2),&
                        fvec(nd), xd(nd), yd(nd)               
    ! Local variables.
    real(kind(1.0d0)), parameter:: kbz=8.617385e-5
    real(kind(1.0d0)):: xx(52+3), maxi, engy, maxt, &
                        bv, xa, xb(nd), expv(nd),&
                        ba, bb, bc
    integer:: i, n0
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
    !
    do i=1, n0/4
        maxi = xx(i)
        engy = xx(i+n0/4)
        maxt = xx(i+2*n0/4)
        bv = xx(i+3*n0/4)
        !
        xa = 2.0*kbz*maxt/engy
        xb = 2.0*kbz*xd/engy
        expv = exp(engy/kbz/xd*(xd-maxt)/maxt)
        !
        fvec = fvec + maxi*(bv**(bv/(bv-1.0)))*expv*&
               ((bv-1.0)*(1.0-xb)*((xd/maxt)**2)*&
               expv+1.0+(bv-1.0)*xa)**(-bv/(bv-1.0))
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_gnr1

subroutine tgcfunc_frt3(nd,n2,pars,fvec,iflag,&
                        xd,yd,lower,upper,bg)
!---------------------------------------------------------
! Subroutine tgcfunc_frt3() is used for calculating
! the residual vector of a optimized TL glow curve 
! according to second-order kinetic.
!---------------------------------------------------------
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
!---------------------------------------------------------
! Author:: Peng Jun, 2020.05.08.
!---------------------------------------------------------
! Dependence:: NO.
!---------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag, bg
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)                 
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(39+3), maxi, engy, maxt, xa, xb(nd), expv(nd),&
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
        xa = 2.0*kbz*maxt/engy
        xb = 2.0*kbz*xd/engy
        expv = exp(engy/kbz/xd*(xd-maxt)/maxt)
        !
        fvec = fvec + 4.0*maxi*expv/(((xd/maxt)**2)*(1.0-xb)*expv+1.0+xa)**2
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_frt3

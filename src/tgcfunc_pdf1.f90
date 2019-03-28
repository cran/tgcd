subroutine tgcfunc_pdf1(nd,n2,pars,fvec,iflag,&
                        xd,yd,lower,upper,bg)
!--------------------------------------------------------------
! Subroutine tgcfunc_pdf1() is used for calculating
! the residual vector of a optimized TL glow curve 
! according to first-order kinetic using weibull distribution.
!--------------------------------------------------------------
!        nd:: input, integer, number of data points.
!        n2:: input, integer, number of pars (<=39+4).
!  pars(n2):: input, real values, pars.
!  fvec(nd):: output, real values, residuals.
!     iflag:: input, integer, working variable.
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observations Y.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real vlaues, upper bounds.
!        bg:: input, integer, subtract background or not,
!             0=no subtraction, 1=subtraction.
!-------------------------------------------------------------
! Author:: Peng Jun, 2019.03.24.
!-------------------------------------------------------------
! Dependence:: NO.
!-------------------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag, bg
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)                 
    ! Local variables.
    real   (kind=8), parameter:: kbz=8.617385e-5
    real   (kind=8):: xx(39+4), maxi, engy, maxt, bv, xa(nd),&
                      ba, bb, bc, bd          
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
    do i=1, n0/3
        maxi = xx(i)
        engy = xx(i+n0/3)
        maxt = xx(i+2*n0/3)
        !
        bv = sqrt((maxt*kbz)**2*242.036*maxt**2/((engy+maxt*kbz)**2-7.0*(maxt*kbz)**2))
        !
        xa = (xd-maxt)/bv
        !
        fvec = fvec + 2.713*maxi*(xa+0.996)**15*exp(-(xa+0.996)**16)
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_pdf1

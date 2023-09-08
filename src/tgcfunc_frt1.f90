subroutine tgcfunc_frt1(nd,n2,pars,fvec,iflag,&
                        xd,yd,lower,upper,bg)
!----------------------------------------------------------
! Subroutine tgcfunc_frt1() is used for calculating
! the residual vector of a optimized TL glow curve 
! according to first-order kinetic (type 1).
!----------------------------------------------------------
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
!----------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!----------------------------------------------------------
! Dependence:: NO.
!----------------------------------------------------------
    ! Arguments.
    implicit none 
    integer:: nd, n2, iflag, bg
    real(kind(1.0d0)):: pars(n2), lower(n2), upper(n2),&
                        fvec(nd), xd(nd), yd(nd)                 
    ! Local variables.
    real(kind(1.0d0)), parameter:: kbz=8.617385e-5, a0=0.267773734, &
                       a1=8.6347608925, a2=18.059016973, a3=8.5733287401, &
                       b0=3.9584969228, b1=21.0996530827, b2=25.6329561486, &
                       b3=9.5733223454
    real(kind(1.0d0)):: xx(39+3), maxi, engy, maxt, xa, fxa, xb(nd), fxb(nd),&
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
    do i=1, n0/3
        maxi = xx(i)
        engy = xx(i+n0/3)
        maxt = xx(i+2*n0/3)
        !
        xa = engy/kbz/maxt
        xb = engy/kbz/xd
        fxa = 1.0-(a0+a1*xa+a2*xa**2+a3*xa**3+xa**4)/&
                  (b0+b1*xa+b2*xa**2+b3*xa**3+xa**4)
        fxb = 1.0-(a0+a1*xb+a2*xb**2+a3*xb**3+xb**4)/&
                  (b0+b1*xb+b2*xb**2+b3*xb**3+xb**4)
        fvec = fvec + maxi*exp(xa-xb)*&
               exp(xa*(fxa-xd/maxt*fxb*exp(xa-xb)))
    end do
    !
    fvec = sqrt(abs(fvec-yd))
    !
    return
    !
end subroutine tgcfunc_frt1

subroutine tgcd(xd,yd,nd,pars,stdp,n2,fmin,&
                message,lower,upper,nstart,mdt)
!---------------------------------------------------
! Subroutine tgcd() is used for thermoluminescence 
! glow curve deconvolution using the 
! Levenberg–Marquardt algorithm.
!---------------------------------------------------
!    xd(nd):: input, real values, observation X.
!    yd(nd):: input, real vlaues, observations Y.
!        nd:: input, integer, number of points.
!  pars(n2):: input/output, paraneters.
!  stdp(n2):: output, real values, errors of pars.
!        n2:: input, integer, number of pars (<=30).
!      fmin:: output, real value, minimized objective.
!   message:: output, integer, error message:
!                     0=success, 1=fail.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real values, upper bounds.
!    nstart:: input, integer, number of trials.
!      mdt:: input, real value, minimum distance
!             between each temperature.
!---------------------------------------------------
! Author:: Peng Jun, 2015.06.01.
!---------------------------------------------------
! Dependence:: subroutine lmtl; subroutine hpSort.
!---------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nd, n2, nstart
    real   (kind=8), intent(in):: xd(nd), yd(nd), mdt
    real   (kind=8), intent(in):: lower(n2), upper(n2)
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: stdp(n2), fmin
    integer(kind=4), intent(out):: message
    ! Local variables.
    real   (kind=8):: ran(n2), ranpars(n2), ranstdp(n2), &
                      ranfmin, minfmin, unifa(n2), unifb(n2),&
                      orderTemp(n2/3), mindist
    integer(kind=4):: i, info, indx(n2/3)
    !
    call lmtl(xd,yd,nd,pars,stdp,n2,&
              fmin,info,lower,upper)
    !
    orderTemp = pars(2*n2/3+1:n2)
    call hpSort(orderTemp, n2/3, indx) 
    mindist = minval(orderTemp(2:n2/3)-&
                     orderTemp(1:n2/3-1))
    !
    if (info==0 .and. mindist>=mdt) then
        message = 0
    else 
        message = 1
    end if
    !
    if (nstart==1)  return
    !
    if (nstart>1) then 
        if (message/=0)  then
            minfmin = 1.0D+20
        else 
            minfmin = fmin          
        end if
        !
        ! INTENS.
        unifa(1:n2/3) = pars(1:n2/3)*0.7
        unifb(1:n2/3) = pars(1:n2/3)*1.0
        ! ENERGY.
        unifa(n2/3+1:2*n2/3) = pars(n2/3+1:2*n2/3)*0.9
        unifb(n2/3+1:2*n2/3) = pars(n2/3+1:2*n2/3)*1.1
        ! TEMPER.
        unifa(2*n2/3+1:n2) = pars(2*n2/3+1:n2)*0.9
        unifb(2*n2/3+1:n2) = pars(2*n2/3+1:n2)*1.1
        !
        ! Try-and-error.
        do i=1,  nstart
            call random_number(ran)
            ranpars = unifa + ran*(unifb-unifa)
            !
            call lmtl(xd,yd,nd,ranpars,ranstdp,n2,&
                      ranfmin,info,lower,upper)
            !
            orderTemp = ranpars(2*n2/3+1:n2)
            call hpSort(orderTemp, n2/3, indx) 
            mindist = minval(orderTemp(2:n2/3)-orderTemp(1:n2/3-1))
            !
            if (info==0 .and. mindist>=mdt)  then
                if (ranfmin<minfmin)  then
                    pars = ranpars
                    stdp = ranstdp
                    fmin = ranfmin
                    minfmin = ranfmin
                    message = 0
                end if
            end if
        end do
    end if
    !
    return
end subroutine tgcd
!
subroutine lmtl(xd,yd,nd,pars,stdp,n2,&
                fmin,message,lower,upper)
!----------------------------------------------------
! Subroutine lmtl() is used for fitting a TL glow 
! curve using the Levenberg-Marquardt algorithm.
!----------------------------------------------------
!    xd(nd):: input, real values, observation X.
!    yd(nd):: input, real vlaues, observations Y.
!        nd:: input, integer, number of points.
!  pars(n2):: input/output, paraneters.
!  stdp(n2):: output, real values, errors of pars.
!        n2:: input, integer, number of pars (<=30).
!      fmin:: output, real value, minimized objective.
!   message:: output, integer, error message:
!                     0=success, 1=fail.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real values, upper bounds.
!------------------------------------------------------
! Author:: Peng Jun, 2015.05.26.
!------------------------------------------------------
! Dependence:: subroutine lmdif1; subroutine tgcfunc;
!              subroutine numHess, subroutine inverse.
!------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nd, n2
    real   (kind=8), intent(in):: xd(nd), yd(nd)
    real   (kind=8), intent(in):: lower(n2), upper(n2)
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: stdp(n2), fmin
    integer(kind=4), intent(out):: message
    !
    integer(kind=4):: info, errorflag, singular, i
    real   (kind=8), parameter:: tol=1.490116D-08
    real   (kind=8):: fvec(nd), hess(n2, n2), diag(n2)
    external:: tgcfunc
    !
    stdp = -99.0
    fmin = -99.0
    !
    call lmdif1(tgcfunc,nd,n2,pars,fvec,&
                tol,info,xd,yd,lower,upper)
    !
    if (info==1 .or. info==2 .or. info==3) then
        message = 0
    else 
        message = 1
        return
    end if
    !
    fmin = sum(fvec**2)
    !
    call numHess(xd,yd,nd,pars,&
                 n2,hess,errorflag)
    if (errorflag/=0) then
        message = 1
        return
    end if
    !
    call inverse(hess, n2, singular)
    if (singular/=0) then
        message = 1
        return
    end if
    !
    do i=1, n2
        diag(i) = hess(i,i) * fmin/real(nd-n2)
    end do
    !if (any(diag<0.0)) then
        !message = 1
        !return
    !end if
    !
    stdp = sqrt(diag)
    !
    return
end subroutine lmtl
!
subroutine tgcfunc(nd,n2,pars,fvec,iflag,&
                   xd,yd,lower,upper)
!---------------------------------------------------
! Subroutine tgcfunc() is used for calculating
! the residual vector of a optimized TL glow curve.
!---------------------------------------------------
!        nd:: input, integer, number of data points.
!        n2:: input, integer, number of pars (<=30).
!  pars(n2):: input, real values, pars.
!  fvec(nd):: output, real values, residuals.
!     iflag:: input, integer, working variable.
!    xd(nd):: input, real values, observations X.
!    yd(nd):: input, real values, observations Y.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real vlaues, upper bounds.
!----------------------------------------------------
! Author:: Peng Jun, 2015.05.26.
!----------------------------------------------------
! Dependence:: NO.
!----------------------------------------------------
    ! Arguments.
    implicit none 
    integer(kind=4):: nd, n2, iflag
    real   (kind=8):: pars(n2), lower(n2), upper(n2),&
                      fvec(nd), xd(nd), yd(nd)                 
    ! Local variables.
    real   (kind=8):: xx(39), maxi, engy, maxt, &
                      xa, fxa, xb(nd), fxb(nd)
    real   (kind=8), parameter:: kbz=8.617385e-5, a0=0.267773734, &
                 a1=8.6347608925, a2=18.059016973, a3=8.5733287401, &
                 b0=3.9584969228, b1=21.0996530827, b2=25.6329561486, &
                 b3=9.5733223454
    integer(kind=4):: i
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
    xx = 0.0
    xx(1:n2) = pars(1:n2)
    !
    fvec = 0.0
    do i=1, n2/3
        maxi = xx(i)
        engy = xx(i+n2/3)
        maxt = xx(i+2*n2/3)
        xa = engy/kbz/maxt
        xb = engy/kbz/xd
        fxa = 1.0-(a0+a1*xa+a2*xa**2+a3*xa**3+xa**4)/&
                  (b0+b1*xa+b2*xa**2+b3*xa**3+xa**4)
        fxb = 1.0-(a0+a1*xb+a2*xb**2+a3*xb**3+xb**4)/&
                  (b0+b1*xb+b2*xb**2+b3*xb**3+xb**4)
        fvec = fvec + maxi*exp(xa-xb)*&
        exp(xa*(fxa-xd/maxt*fxb*exp(xa-xb)))
    end do
    fvec = sqrt(abs(fvec-yd))
    return
    !
end subroutine tgcfunc

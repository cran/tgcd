subroutine numHess(xd,yd,nd,pars,&
                   np,hess,iflag)
!-------------------------------------------------------------
! Subroutine numHess() is used for approximating the 
! Hessian matrix of first-order Thermoluminescence glow model
! using a numerical differention method called "Richardson".
!-------------------------------------------------------------
!      xd(nd):: input, real values, observations x.
!      yd(nd):: input, real values, observations y.
!          nd:: input, integer, number of data points.                   
!    pars(np):: input, real values, parameters.
!          np:: input, integer, number of pars (<=30).
! hess(np,np):: output, real values, the Hessian matrix.
!       iflag:: output, integer, 0=success; 1=fail.
!--------------------------------------------------------
! Author:: Peng Jun, 2015.05.26.
!--------------------------------------------------------
! Dependence:: inner function func.
!--------------------------------------------------------
! Reference:: Paul Gilbert and Ravi Varadhan (2012). 
!             numDeriv: Accurate Numerical
!             Derivatives. R package version 2012.9-1. 
!             Availiable at:
!             http://CRAN.R-project.org/package=numDeriv
!--------------------------------------------------------
! NOTE: THIS SUBROUTINE IS BASED ON FUNCTION 
!       "genD" IN R PACKAGE numDeriv.
!--------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: nd, np
    real   (kind=8), intent(in):: xd(nd), yd(nd),& 
                                  pars(np)
    real   (kind=8), intent(out):: hess(np,np)
    integer(kind=4), intent(out):: iflag
    ! Local variables.
    real   (kind=8), parameter:: eps=1.0D-04,&
                                 d=0.1D+00,&
                                 zTol=1.781029D-05
    integer(kind=4), parameter:: r=4
    real   (kind=8):: h0(np), h(np), Hdiag(np),& 
                      incr1(np), incr2(np),&
                      Dd(np*(np+3)/2), Daprox(r),& 
                      Haprox(r)
    real   (kind=8):: f0, f1, f2, m4
    integer(kind=4):: i, j, k, m, u
    real   (kind=8), parameter:: kbz=8.617385e-5, a0=0.267773734, &
                 a1=8.6347608925, a2=18.059016973, a3=8.5733287401, &
                 b0=3.9584969228, b1=21.0996530827, b2=25.6329561486, &
                 b3=9.5733223454
    real   (kind=8):: xx(39), maxi, engy, maxt, &
                      xa, fxa, xb(nd), fxb(nd)
    !
    iflag = 0
    hess = 0.0D+00
    f0 = func(pars)
    if(f0 .ne. f0 .or. f0+1.0D+00==f0) then
        iflag = 1
        return
    end if
    !
    do i=1, np
        h0(i) = dabs(d*pars(i))
        if (dabs(pars(i))<zTol) h0(i) = h0(i) + eps
    end do
    !
    Dd = 0.0D+00
    Hdiag = 0.0D+00
    Daprox = 0.0D+00
    Haprox = 0.0D+00
    !
    do i=1, np
        h = h0
        incr1 = 0.0D+00 
        do k=1, r
            incr1(i) = h(i)
            f1 = func(pars+incr1) 
            f2 = func(pars-incr1)
            Daprox(k) = (f1-f2) / (2.0D+00*h(i))
            Haprox(k) = (f1-2.0D+00*f0+f2) / (h(i))**2
            h = h / 2.0D+00
        end do
        !
        do m=1, r-1
            m4 = (4.0D+00)**m
            do k=1, r-m
                Daprox(k) = (Daprox(k+1)*m4-Daprox(k))/&
                            (m4-1.0D+00)
                Haprox(k) = (Haprox(k+1)*m4-Haprox(k))/&
                            (m4-1.0D+00)
            end do
        end do
        Dd(i) = Daprox(1)
        Hdiag(i) = Haprox(1)
    end do
    !
    u = np
    do i=1, np
        do j=1, i
            u = u + 1
            if (i==j) then 
                Dd(u) = Hdiag(i)
            else 
                h = h0
                incr1 = 0.0D+00 
                incr2 = 0.0D+00
                do k=1, r
                    incr1(i) = h(i)
                    incr2(j) = h(j)
                    f1 = func(pars+incr1+incr2)
                    f2 = func(pars-incr1-incr2)
                    Daprox(k) = (f1-2.0D+00*f0+f2-&
                                 Hdiag(i)*(h(i))**2-&
                                 Hdiag(j)*(h(j))**2)/&
                                (2.0D+00*h(i)*h(j))
                    h = h / 2.0D+00
                end do
                !
                do m=1, r-1
                    m4 = (4.0D+00)**m
                    do k=1, r-m
                        Daprox(k) = (Daprox(k+1)*m4-Daprox(k))/& 
                                    (m4-1.0D+00)
                    end do 
                end do
                Dd(u) = Daprox(1)
            end if
        end do
    end do
    ! 
    u = np
    do i=1, np
        do j=1, i
            u = u + 1
            hess(i,j) = Dd(u)
        end do
    end do
    hess = hess + transpose(hess)
    do i=1, np
        hess(i,i) = hess(i,i) / 2.0D+00
    end do
    if(any(hess .ne. hess) .or.&
       any(hess+1.0D+00==hess))  iflag=1
    !
    return
    !
    ! Inner function func.
    !---------------------
    contains
        function func(x)
            implicit none
            real   (kind=8):: x(np), func
            real   (kind=8):: vec(nd)
            integer(kind=4):: i
            !
            xx = 0.0D+00
            xx(1:np) = x
            !
            vec = 0.0
            do i=1, np/3
                maxi = xx(i)
                engy = xx(i+np/3)
                maxt = xx(i+2*np/3)
                xa = engy/kbz/maxt
                xb = engy/kbz/xd
                fxa = 1.0-(a0+a1*xa+a2*xa**2+a3*xa**3+xa**4)/&
                          (b0+b1*xa+b2*xa**2+b3*xa**3+xa**4)
                fxb = 1.0-(a0+a1*xb+a2*xb**2+a3*xb**3+xb**4)/&
                          (b0+b1*xb+b2*xb**2+b3*xb**3+xb**4)
                vec = vec + maxi*exp(xa-xb)*&
                exp(xa*(fxa-xd/maxt*fxb*exp(xa-xb)))
            end do
            vec = abs(vec-yd)
            func = sum(vec)
            return
        end function func 
        !----------------      
end subroutine numHess

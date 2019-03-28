subroutine calcLv(ax,bx,alpha,maxt,engy,Lv,fmin)
!-----------------------------------------------------
! Subroutine calcLv() is used for calculating 
! the Am value of a mix-order kinetic model.
!-----------------------------------------------------
!       ax:: input, real value, lower limit.
!       bx:: input, real value, upper limit.
!    alpha:: input, real value, the alpha value.
!     maxt:: input, real value, the Tm value.
!     engy:: input, real value, the E value.
!       Lv:: output, real value, resulting value.
!     fmin:: output, real value, minimized objective.
!-----------------------------------------------------
! Author:: Peng Jun, 2019.03.14.
!-----------------------------------------------------
! Dependence:: inner function fcn.
!-----------------------------------------------------
! NOTE: THIS ROUTINE IS BASED ON THE FREE FORTRAN 77 
!       CODE AT: http://www.netlib.org/fmm/fmin.f .
! Part of the R package, http://www.R-project.org .
!-----------------------------------------------------
    implicit none
    ! Arguments.
    real   (kind=8), intent(in):: ax, bx, alpha, maxt, engy
    real   (kind=8), intent(out):: Lv, fmin
    ! Local variables.
    real   (kind=8):: a, b, c, d, e, p, q, r, u, v, w, x
    real   (kind=8):: t2, fu, fv, fw, fx, xm, eps, tol1, tol3
    real   (kind=8), parameter:: tol=1.490116D-08, kbz=8.617385D-5
    !
    c = (3.0D+00 - dsqrt(5.0D+00)) * 0.5D+00
    eps = EPSILON(0.0D+00) 
    tol1 = eps + 1.0D+00
    eps = dsqrt(eps)

    a = ax
    b= bx
    v = a + c * (b - a)
    w = v
    x = v

    d = 0.0D+00
    e = 0.0D+00
    fx = fcn(x)
    fv = fx
    fw = fx
    tol3 = tol/3.0D+00

    do 
        xm = (a + b) * 0.5D+00
        tol1 = eps * dabs(x) + tol3
        t2 = tol1 * 2.0D+00

        if ( dabs(x - xm) <= t2 - (b - a) * 0.5D+00 ) exit

        p = 0.0D+00
        q = 0.0D+00
        r = 0.0D+00
         
        if (dabs(e) > tol1) then
            r = (x - w) * (fx - fv)
            q = (x - v)*(fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0D+00

            if (q > 0.0D+00) then
                p = -p
            else 
                q = -q
            end if
            r = e
            e = d
        end if
        

        if ((dabs(p) >= dabs(q * 0.5D+00 * r)) .or. &
            (p <= q * (a - x)) .or. (p >= q * (b - x))) then
            if (x < xm) then
                e = b - x
            else 
                e = a - x
            end if
            d = c * e
        else 
            d = p / q
            u = x + d
            
            if ((u - a < t2) .or. (b - u < t2)) then
                d = tol1
                if (x >= xm) d = -d
            end if
        end if

        if (dabs(d) >= tol1) then
            u = x + d
        else if (d > 0.0D+00) then
            u = x + tol1
        else 
            u = x - tol1
        end if

        fu = fcn(u)

        if (fu <= fx) then
            if (u < x) then
                b = x
            else 
                a = x
            end if

            v = w; w = x; x = u
            fv = fw; fw = fx; fx = fu

        else 
            if (u < x) then
                a = u
            else 
                b = u
            end if
                 
            if ((fu <= fw) .or. (w == x)) then
                v = w; fv = fw
                w = u; fw = fu
            else if ((fu <= fv) .or. (v == x) .or. (v == w)) then
                v = u; fv = fu
            end if
        end if

    end do
    
    Lv = x
    fmin = fcn(x)
    
    return   
    !
    contains
    ! Inner function fcn.
    function fcn(x)
        implicit none
        real   (kind=8):: fcn, x
        !
        fcn = (alpha-(x-1.0)*exp((2.0-x)/x*(1.0-2.0*kbz*maxt/engy)))**2
        !
        ! Test for Inf or NaN.
        if ((fcn .ne. fcn) .or. (fcn .gt. huge(0.0D+00))) fcn = huge(0.0D+00)
        !
        return
    end function fcn
end subroutine calcLv

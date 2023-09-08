function bisect ( xx, nb, ner, l )

!*****************************************************************************80
!
!! BISECT approximates the W function using bisection.
!
!  Discussion:
!
!    The parameter TOL, which determines the accuracy of the bisection
!    method, is calculated using NBITS (assuming the final bit is lost
!    due to rounding error).
!
!    N0 is the maximum number of iterations used in the bisection
!    method.
!
!    For XX close to 0 for Wp, the exponential approximation is used.
!    The approximation is exact to O(XX^8) so, depending on the value
!    of NBITS, the range of application of this formula varies. Outside
!    this range, the usual bisection method is used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XX, the argument.
!
!    Input, integer ( kind = 4 ) NB, indicates the branch of the W function.
!    0, the upper branch;
!    nonzero, the lower branch.
!
!    Output, integer ( kind = 4 ) NER, the error flag.
!    0, success;
!    1, the routine did not converge.  Perhaps reduce NBITS and try again.
!
!    Input, integer ( kind = 4 ) L, the offset indicator.
!    1, XX represents the offset of the argument from -exp(-1).
!    not 1, XX is the actual argument.
!
!    Output, real ( kind = 8 ) BISECT, the value of W(X), as determined
!
  implicit none

  integer, parameter:: dbdbx=kind(1.0d0)

  real ( dbdbx ) bisect
  real ( dbdbx ) crude
  real ( dbdbx ) d
  real ( dbdbx ) f
  real ( dbdbx ) fd
  integer i
  integer l
  integer n0
  parameter ( n0 = 500 )
  integer nb
  integer, save :: nbits = 0
  integer ner
  real ( dbdbx ) r
  real ( dbdbx ) test
  real ( dbdbx ) tol
  real ( dbdbx ) u
  real ( dbdbx ) x
  real ( dbdbx ) xx

  bisect = 0.0D+00
  ner = 0

  if ( nbits == 0 ) then
    call nbits_compute ( nbits )
  end if

  if ( l == 1 ) then
    x = xx - exp ( -1.0D+00 )
  else
    x = xx
  end if

  if ( nb == 0 ) then

    test = 1.0D+00 / ( 2.0D+00 ** nbits ) ** ( 1.0D+00 / 7.0D+00 )

    if ( abs ( x ) < test ) then

      bisect = x &
        * exp ( - x &
        * exp ( - x &
        * exp ( - x &
        * exp ( - x &
        * exp ( - x &
        * exp ( - x ))))))

      return

    else

      u = crude ( x, nb ) + 1.0D-03
      tol = abs ( u ) / 2.0D+00 ** nbits
      d = max ( u - 2.0D-03, -1.0D+00 )

      do i = 1, n0

        r = 0.5D+00 * ( u - d )
        bisect = d + r
!
!  Find root using w*exp(w)-x to avoid ln(0) error.
!
        if ( x < exp ( 1.0D+00 ) ) then

          f = bisect * exp ( bisect ) - x
          fd = d * exp ( d ) - x
!
!  Find root using ln(w/x)+w to avoid overflow error.
!
        else

          f = log ( bisect / x ) + bisect
          fd = log ( d / x ) + d

        end if

        if ( f == 0.0D+00 ) then
          return
        end if

        if ( abs ( r ) <= tol ) then
          return
        end if

        if ( 0.0D+00 < fd * f ) then
          d = bisect
        else
          u = bisect
        end if

      end do

    end if

  else

    d = crude ( x, nb ) - 1.0D-03
    u = min ( d + 2.0D-03, -1.0D+00 )
    tol = abs ( u ) / 2.0D+00 ** nbits

    do i = 1, n0

      r = 0.5D+00 * ( u - d )
      bisect = d + r
      f = bisect * exp ( bisect ) - x

      if ( f == 0.0D+00 ) then
        return
      end if

      if ( abs ( r ) <= tol ) then
        return
      end if

      fd = d * exp ( d ) - x

      if ( 0.0D+00 < fd * f ) then
        d = bisect
      else
        u = bisect
      end if

    end do

  end if
!
!  The iteration did not converge.
!
  ner = 1

  return
end
function crude ( xx, nb )

!*****************************************************************************80
!
!! CRUDE returns a crude approximation for the W function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XX, the argument.
!
!    Input, integer ( kind = 4 ) NB, indicates the desired branch.
!    * 0, the upper branch;
!    * nonzero, the lower branch.
!
!    Output, real ( kind = 8 ) CRUDE, the crude approximation to W at XX.
!
  implicit none

  integer, parameter:: dbdbx=kind(1.0d0)

  real ( dbdbx ) an2
  real ( dbdbx ) c13
  real ( dbdbx ) crude
  real ( dbdbx ) em
  real ( dbdbx ) em2
  real ( dbdbx ) em9
  real ( dbdbx ) eta
  integer init
  integer nb
  real ( dbdbx ) reta
  real ( dbdbx ) s2
  real ( dbdbx ) s21
  real ( dbdbx ) s22
  real ( dbdbx ) s23
  real ( dbdbx ) t
  real ( dbdbx ) ts
  real ( dbdbx ) xx
  real ( dbdbx ) zl

  save c13
  save em
  save em2
  save em9
  save init
  save s2
  save s21
  save s22
  save s23

  data init / 0 /

  crude = 0.0D+00
!
!  Various mathematical constants.
!
  if ( init == 0 ) then
    init = 1
    em = - exp ( -1.0D+00 )
    em9 = - exp ( -9.0D+00 )
    c13 = 1.0D+00 / 3.0D+00
    em2 = 2.0D+00 / em
    s2 = sqrt ( 2.0D+00 )
    s21 = 2.0D+00 * s2 - 3.0D+00
    s22 = 4.0D+00 - 3.0D+00 * s2
    s23 = s2 - 2.0D+00
  end if
!
!  Crude Wp.
!
  if ( nb == 0 ) then

    if ( xx <= 20.0D+00 ) then
      reta = s2 * sqrt ( 1.0D+00 - xx / em )
      an2 = 4.612634277343749D+00 * sqrt ( sqrt ( reta + &
        1.09556884765625D+00 ) )
      crude = reta / ( 1.0D+00 + reta / ( 3.0D+00 &
        + ( s21 * an2 + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0D+00
    else
      zl = log ( xx )
      crude = log ( xx / log ( xx &
        / zl ** exp ( -1.124491989777808D+00 / &
        ( 0.4225028202459761D+00 + zl ))))
    end if

  else
!
!  Crude Wm.
!
    if ( xx <= em9 ) then
      zl = log ( -xx )
      t = -1.0D+00 - zl
      ts = sqrt ( t )
      crude = zl - ( 2.0D+00 * ts ) / ( s2 + ( c13 - t &
        / ( 2.7D+02 + ts * 127.0471381349219D+00 ) ) * ts )
    else
      zl = log ( -xx )
      eta = 2.0D+00 - em2 * xx
      crude = log ( xx / log ( - xx / ( ( 1.0D+00 &
        - 0.5043921323068457D+00 * ( zl + 1.0D+00 ) ) &
        * ( sqrt ( eta ) + eta / 3.0D+00 ) + 1.0D+00 ) ) )
     end if

  end if

  return
end
subroutine nbits_compute ( nbits )

!*****************************************************************************80
!
!! NBITS_COMPUTE computes the mantissa length minus one.
!
!  Discussion:
!
!    NBITS is the number of bits (less 1) in the mantissa of the
!    floating point number number representation of your machine.
!    It is used to determine the level of accuracy to which the W
!    function should be calculated.
!
!    Most machines use a 24-bit matissa for single precision and
!    53-56 bits for real ( kind = 8 ). The IEEE standard is 53
!    bits. The Fujitsu VP2200 uses 56 bits. Long word length
!    machines vary, e.g., the Cray X/MP has a 48-bit mantissa for
!    single precision.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) NBITS, the mantissa length, in bits, 
!    minus one.
!
  implicit none

  integer, parameter:: dbdbx=kind(1.0d0)

  real ( dbdbx ) b
  integer i
  integer nbits
  real ( dbdbx ) v

  nbits = 0

  b = 1.0D+00

  do

    b = b / 2.0D+00
    v = b + 1.0D+00

    if ( v == 1.0D+00 ) then
      return
    end if

    nbits = nbits + 1

  end do

  return
end
function wapr ( x, nb, nerror, l )

!*****************************************************************************80
!
!! WAPR approximates the W function.
!
!  Discussion:
!
!    The call will fail if the input value X is out of range.
!    The range requirement for the upper branch is:
!      -exp(-1) <= X.
!    The range requirement for the lower branch is:
!      -exp(-1) < X < 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 June 2014
!
!  Author:
!
!    Original FORTRAN77 version by Andrew Barry, S. J. Barry, 
!    Patricia Culligan-Hensley.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Andrew Barry, S. J. Barry, Patricia Culligan-Hensley,
!    Algorithm 743: WAPR - A Fortran routine for calculating real 
!    values of the W-function,
!    ACM Transactions on Mathematical Software,
!    Volume 21, Number 2, June 1995, pages 172-181.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument.
!
!    Input, integer ( kind = 4 ) NB, indicates the desired branch.
!    * 0, the upper branch;
!    * nonzero, the lower branch.
!
!    Output, integer ( kind = 4 ) NERROR, the error flag.
!    * 0, successful call.
!    * 1, failure, the input X is out of range.
!
!    Input, integer ( kind = 4 ) L, indicates the interpretation of X.
!    * 1, X is actually the offset from -(exp-1), so compute W(X-exp(-1)).
!    * not 1, X is the argument; compute W(X);
!
!    Output, real ( kind = 8 ) WAPR, the approximate value of W(X).
!
  implicit none

  integer, parameter:: dbdbx=kind(1.0d0)

  real ( dbdbx ) an2
  real ( dbdbx ) an3
  real ( dbdbx ) an4
  real ( dbdbx ) an5
  real ( dbdbx ) an6
  real ( dbdbx ) c13
  real ( dbdbx ) c23
  real ( dbdbx ) d12
  real ( dbdbx ) delx
  real ( dbdbx ) em
  real ( dbdbx ) em2
  real ( dbdbx ) em9
  real ( dbdbx ) eta
  integer i
  integer init
  integer l
  integer m
  integer nb
  integer nbits
  integer nerror
  integer niter
  real ( dbdbx ) reta
  real ( dbdbx ) s2
  real ( dbdbx ) s21
  real ( dbdbx ) s22
  real ( dbdbx ) s23
  real ( dbdbx ) t
  real ( dbdbx ) tb
  real ( dbdbx ) tb2
  real ( dbdbx ) temp
  real ( dbdbx ) temp2
  real ( dbdbx ) ts
  real ( dbdbx ) wapr
  real ( dbdbx ) x
  real ( dbdbx ) x0
  real ( dbdbx ) x1
  real ( dbdbx ) xx
  real ( dbdbx ) zl
  real ( dbdbx ) zn

  save an3
  save an4
  save an5
  save an6
  save c13
  save c23
  save d12
  save em
  save em2
  save em9
  save init
  save nbits
  save niter
  save s2
  save s21
  save s22
  save s23
  save tb
  save tb2
  save x0
  save x1

  data init / 0 /
  data niter / 1 /

  wapr = 0.0D+00
  nerror = 0

  if ( init == 0 ) then

    init = 1

    call nbits_compute ( nbits )

    if ( 56 <= nbits ) then
      niter = 2
    end if
!
!  Various mathematical constants.
!
    em = -exp ( -1.0D+00 )
    em9 = -exp ( -9.0D+00 )
    c13 = 1.0D+00 / 3.0D+00
    c23 = 2.0D+00 * c13
    em2 = 2.0D+00 / em
    d12 = -em2
    tb = 0.5D+00 ** nbits
    tb2 = sqrt ( tb )
    x0 = tb ** ( 1.0D+00 / 6.0D+00 ) * 0.5D+00
    x1 = ( 1.0D+00 - 17.0D+00 * tb ** ( 2.0D+00 / 7.0D+00 ) ) * em
    an3 = 8.0D+00 / 3.0D+00
    an4 = 135.0D+00 / 83.0D+00
    an5 = 166.0D+00 / 39.0D+00
    an6 = 3167.0D+00 / 3549.0D+00
    s2 = sqrt ( 2.0D+00 )
    s21 = 2.0D+00 * s2 - 3.0D+00
    s22 = 4.0D+00 - 3.0D+00 * s2
    s23 = s2 - 2.0D+00

  end if

  if ( l == 1 ) then

    delx = x

    if ( delx < 0.0D+00 ) then
      nerror = 1
      return
      !!!write ( *, '(a)' ) ''
      !!!write ( *, '(a)' ) 'WAPR - Fatal error!'
      !!!write ( *, '(a)' ) '  The offset X is negative.'
      !!!write ( *, '(a)' ) '  It must be nonnegative.'
      !!!stop 1
    end if

    xx = x + em

  else

    if ( x < em ) then
      nerror = 1
      return
    else if ( x == em ) then
      wapr = -1.0D+00
      return
    end if

    xx = x
    delx = xx - em

  end if

  if ( nb == 0 ) then
!
!  Calculations for Wp.
!
    if ( abs ( xx ) <= x0 ) then
      wapr = xx / ( 1.0D+00 + xx / ( 1.0D+00 + xx &
        / ( 2.0D+00 + xx / ( 0.6D+00 + 0.34D+00 * xx ))))
      return
    else if ( xx <= x1 ) then
      reta = sqrt ( d12 * delx )
      wapr = reta / ( 1.0D+00 + reta / ( 3.0D+00 + reta / ( reta &
        / ( an4 + reta / ( reta * an6 + an5 ) ) + an3 ) ) ) &
        - 1.0D+00
      return
    else if ( xx <= 20.0D+00 ) then
      reta = s2 * sqrt ( 1.0D+00 - xx / em )
      an2 = 4.612634277343749D+00 * sqrt ( sqrt ( reta + &
        1.09556884765625D+00 ))
      wapr = reta / ( 1.0D+00 + reta / ( 3.0D+00 + ( s21 * an2 &
        + s22 ) * reta / ( s23 * ( an2 + reta )))) - 1.0D+00
    else
      zl = log ( xx )
      wapr = log ( xx / log ( xx &
        / zl ** exp ( -1.124491989777808D+00 / &
        ( 0.4225028202459761D+00 + zl ))))
    end if
!
!  Calculations for Wm.
!
  else

    if ( 0.0D+00 <= xx ) then
      nerror = 1
      return
    else if ( xx <= x1 ) then
      reta = sqrt ( d12 * delx )
      wapr = reta / ( reta / ( 3.0D+00 + reta / ( reta / ( an4 &
        + reta / ( reta * an6 - an5 ) ) - an3 ) ) - 1.0D+00 ) - 1.0D+00
      return
    else if ( xx <= em9 ) then
      zl = log ( -xx )
      t = -1.0D+00 - zl
      ts = sqrt ( t )
      wapr = zl - ( 2.0D+00 * ts ) / ( s2 + ( c13 - t &
        / ( 270.0D+00 + ts * 127.0471381349219D+00 )) * ts )
    else
      zl = log ( -xx )
      eta = 2.0D+00 - em2 * xx
      wapr = log ( xx / log ( -xx / ( ( 1.0D+00 &
        - 0.5043921323068457D+00 * ( zl + 1.0D+00 ) ) &
        * ( sqrt ( eta ) + eta / 3.0D+00 ) + 1.0D+00 )))
    end if

  end if

  do i = 1, niter
    zn = log ( xx / wapr ) - wapr
    temp = 1.0D+00 + wapr
    temp2 = temp + c23 * zn
    temp2 = 2.0D+00 * temp * temp2
    wapr = wapr * ( 1.0D+00 + ( zn / temp ) * ( temp2 - zn ) &
      / ( temp2 - 2.0D+00 * zn ) )
  end do

  return
end

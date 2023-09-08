subroutine comb_next ( n, k, a, done )

!*****************************************************************************80
!
!! COMB_NEXT computes combinations of K things out of N.
!
!  Discussion:
!
!    The combinations are computed one at a time, in lexicographical order.
!
!    10 April 2009: Thanks to "edA-qa mort-ora-y" for supplying a
!    correction to this code!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Charles Mifsud,
!    Algorithm 154:
!    Combination in Lexicographic Order,
!    Communications of the ACM,
!    March 1963.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the total number of things.
!
!    Input, integer ( kind = 4 ) K, the number of things in each combination.
!
!    Input/output, integer ( kind = 4 ) A(K), contains the list of elements in
!    the current combination.
!
!    Input/output, logical DONE.  On first call, set DONE to TRUE,
!    and thereafter, its input value should be the output value from
!    the previous call.  The output value will normally be FALSE,
!    indicating that there are further combinations that can be
!    returned.  When DONE is returned TRUE, the sequence is exhausted.
!
  implicit none

  integer k

  integer a(k)
  logical done
  integer i
  integer j
  integer n

  if ( done ) then

    if ( k <= 0 ) then
      return
    end if

    call i4vec_indicator ( k, a )

    done = .false.

  else

    if ( a(k) < n ) then
      a(k) = a(k) + 1
      return
    end if

    do i = k, 2, -1

      if ( a(i-1) < n-k+i-1 ) then

        a(i-1) = a(i-1) + 1

        do j = i, k
          a(j) = a(i-1) + j - ( i-1 )
        end do

        return

      end if

    end do

    done = .true.

  end if

  return
end

subroutine i4vec_indicator ( n, a )

!*****************************************************************************80
!
!! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of elements of A.
!
!    Output, integer ( kind = 4 ) A(N), the array to be initialized.
!
  implicit none

  integer n

  integer a(n)
  integer i

  do i = 1, n
    a(i) = i
  end do

  return
end

subroutine inverse(mat1,n,singular)
!---------------------------------------------------------
! Subroutine inverse() is used for inverting a square
! matrix using Gaussian-Jordan elimation by full pivoting.
!---------------------------------------------------------
! mat1(n,n):: input/output, real values, a matrix.
!         n:: input, integer, row number of the matrix.
!  singular::  output, integer, 0=non-singular, else 1.
!---------------------------------------------------------
! Author:: Peng Jun, 2014.09.03.
!---------------------------------------------------------
! Dependence:: subroutine gjordan.------------------------
!---------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: n
    real   (kind=8), intent(inout):: mat1(n,n)
    integer(kind=4), intent(out):: singular
    ! Local variables.
    real   (kind=8):: mat2(n,n)
    integer(kind=4):: i
    !
    mat2 = 0.0
    do i=1, n
        mat2(i,i) = 1.0 
    end do
    !
    call gjordan(mat1,mat2,n,n,singular)
    mat1 = mat2
    !
    return
end subroutine inverse
!
subroutine gjordan(mat1,mat2,n,m,singular)
!---------------------------------------------------------
! Subroutine gjordan() is used for performing
! Gaussian-Jordan elimation by full pivoting.
!---------------------------------------------------------
! mat1(n,n):: input/output, real values, a square matrix.
! mat2(n,m):: input/output, real values, a matrix.
!         n:: input, integer, row number of mat1, n>=m.
!         m:: input, integer, column number of mat2, m<=n.
!  singular:: output, integer, 0=non-singular, else 1.
!---------------------------------------------------------
! Author:: Peng Jun, 2014.09.03.
!---------------------------------------------------------
! Dependence:: NO.----------------------------------------
!---------------------------------------------------------
    implicit none
    ! Arguments.
    integer(kind=4), intent(in):: n, m
    real   (kind=8), intent(inout):: mat1(n,n), mat2(n,m)
    integer(kind=4), intent(out):: singular
    ! Local variables.
    real    (kind=8), parameter:: sgTol=1.0D-10
    real    (kind=8):: qq, Dd(n)
    integer(kind=4):: i, j, k, is, js(n)
    !
    singular = 0
    !
    do k=1, n
        qq = 0.0
        do i=k, n
            do j=k, n
                if (abs(mat1(i,j))>=qq) then
                    qq = abs(mat1(i,j))
                    is = i
                    js(k) = j
                end if
            end do
        end do
        !
        if (abs(qq)<sgTol) then
            singular = 1
            return
        end if
        !
        Dd(k:n) = mat1(k,k:n)
        mat1(k,k:n) = mat1(is,k:n)
        mat1(is,k:n) = Dd(k:n) 
        !
        Dd(1:m) = mat2(k,1:m)
        mat2(k,1:m) = mat2(is,1:m)
        mat2(is,1:m) = Dd(1:m)
        !
        Dd(1:n) = mat1(1:n,k)
        mat1(1:n,k) = mat1(1:n,js(k))
        mat1(1:n,js(k)) = Dd(1:n)
        !
        mat1(k,k+1:n) = mat1(k,k+1:n) / mat1(k,k)
        mat2(k,1:m) = mat2(k,1:m) / mat1(k,k)
        !
        do i=1, n
            if (i/=k) then
                mat1(i,k+1:n) = mat1(i,k+1:n) -&
                                mat1(i,k)*mat1(k,k+1:n)
                mat2(i,1:m) = mat2(i,1:m) -&
                              mat1(i,k)*mat2(k,1:m)
            end if
        end do
    end do
    !
    do k=n, 1, -1
        Dd(1:m) = mat2(k,1:m)
        mat2(k,1:m) = mat2(js(k),1:m)
        mat2(js(k),1:m) = Dd(1:m)
    end do
    !
    return
end subroutine gjordan       

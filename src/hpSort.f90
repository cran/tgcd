subroutine hpSort(vec,n,indx) 
!---------------------------------------------------------------------
! Subroutine hpSort() is used for sorting a one-dimension array
! (a vector) in ascending order using the Heapsort algorithm.
!---------------------------------------------------------------------
! vec(n) : input/output, real values, the vector to be sorted.
!      n : input, integer, the dimension of the vector.
!   indx : output, integers, the orders of the elements.
! --------------------------------------------------------------------
! Author: Peng Jun, 2023.09.07.
!---------------------------------------------------------------------
! Dependence:: No.----------------------------------------------------
!---------------------------------------------------------------------
! Reference: Press et al, 1986. Numberic recipes in Fortran 77, 
!            the Art of Scientific Computing, second edition. 
! NOTE: THIS SUBROUTINE IS REMODIFIED FROM PAGE.329 IN Press et al.
! --------------------------------------------------------------------
    implicit none
    integer, intent(in):: n
    real(kind(1.0d0)), intent(inout):: vec(n)
    integer, intent(out):: indx(n)
    ! Local variables
    integer:: i, j, k, ir, rraIndx
    real(kind(1.0d0)):: rra
    !
    if (n<2) then 
        indx = 1
        return
    end if
    !
    do i= 1, n
        indx(i) = i
    end do
    !
    k = n/2 + 1
    ir = n
    100 continue
        if (k>1) then
            k = k - 1
            rra = vec(k)
            rraIndx = indx(k)
        else 
            rra = vec(ir)
            rraIndx = indx(ir)
            vec(ir) = vec(1)
            indx(ir) = indx(1)
            ir = ir - 1
            if (ir==1) then
                vec(1) = rra
                indx(1) = rraIndx
                return
            end if
        end if
        i = k
        j = k + k
    200 if (j<=ir) then
        if (j<ir) then
            if (vec(j)<vec(j+1)) j=j+1
        end if
        if (rra<vec(j)) then
            vec(i) = vec(j)
            indx(i) = indx(j)
            i = j
            j = j + j
        else 
            j = ir + 1
        end if
        go to 200
    end if
    vec(i) = rra 
    indx(i) = rraIndx
    go to 100
    !
    return
end subroutine hpSort

subroutine calcShape(nd,np,xd,maty,matsp,flag)
!------------------------------------------------------------------------
! Subroutine calcShape() is used for calculating the
! shape factors for a series of glow peaks.
!------------------------------------------------------------------------
!          nd:: integer, input, number of data points.
!          np:: integer, input, number of glow peaks.
!      xd(nd):: real values, input, temperature values.
! maty(nd,np):: real values, input, TL values for each component.
! matsp(np,7):: real value, output, calculated shape factors.
!        flag:: integer, output, error flag, 0=success, 1=failure.
!------------------------------------------------------------------------
! Author:: Peng jun, 2019.03.26.
!------------------------------------------------------------------------
! Dependence:: No.
!------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nd, np
    real   (kind=8), intent(in):: xd(nd), maty(nd,np)
    real   (kind=8), intent(out):: matsp(np,7)
    integer(kind=4), intent(out):: flag
    !
    integer(kind=4):: i, j, indx, loc1
    real   (kind=8):: yd(nd), hymax, T1, T2, x0, y0, x1, y1
    !  
    matsp= -99.0
    flag = 0
    !
    do i=1, np
        yd = maty(:,i)
        hymax = maxval(yd)/2.0
        indx = maxloc(yd,dim=1)
        !
        ! Calculate T1 for each component.
        loc1 = -99
        do j=1, indx-1
            if (hymax>=yd(j) .and. hymax<=yd(j+1)) then
                loc1 = j
            end if 
        end do
        !
        if (loc1==-99) then
            flag = 1
            return
        end if
        !
        x0 = yd(loc1)
        y0 = xd(loc1)
        x1 = yd(loc1+1)
        y1 = xd(loc1+1)
        T1 = y0*(hymax-x1)/(x0-x1) + y1*(hymax-x0)/(x1-x0)
        !
        ! Calculate T2 for each component.
        loc1 = -99
        do j=indx, nd-1
            if (hymax<=yd(j) .and. hymax>=yd(j+1)) then
                loc1 = j
            end if 
        end do
        !
        if (loc1==-99) then
            flag = 1
            return
        end if
        !
        x0 = yd(loc1)
        y0 = xd(loc1)
        x1 = yd(loc1+1)
        y1 = xd(loc1+1)
        T2 = y0*(hymax-x1)/(x0-x1) + y1*(hymax-x0)/(x1-x0)
        !
        !
        matsp(i,1) = T1
        matsp(i,2) = T2
        matsp(i,3) = xd(indx)
        matsp(i,4) = xd(indx) - T1
        matsp(i,5) = T2 - xd(indx)
        matsp(i,6) = T2 - T1
        matsp(i,7) = (T2 - xd(indx))/(T2 - T1)
    end do
    !
    return
    !
end subroutine calcShape

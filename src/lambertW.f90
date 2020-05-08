subroutine lambertW(xx,v,ner)
!-------------------------------------------------------------------------------------------
! Subroutine lambertW() is used for calculating the lower branch of the lambert-W function.
!-------------------------------------------------------------------------------------------
!  xx:: input, real value, value from which the lambert-W function will be calculated.
!   v:: output, real value, the calculated lambert-W value.
! ner:: output, integer, error message generated during the calcualtion, 0=success, 1=failed.
!--------------------------------------------------------------------------------------------
! Author:: Peng Jun, 2020.05.07.
!--------------------------------------------------------------------------------------------
! Dependence:: function bisect.
!--------------------------------------------------------------------------------------------
    implicit none
    real   (kind=8), intent(in):: xx
    real   (kind=8), intent(out):: v
    integer(kind=4), intent(out):: ner
    !
    integer(kind=4):: nb, l
    real   (kind=8):: bisect
    !
    nb = -1
    l = 2
    !
    v = bisect(xx, nb, ner, l)
    !
    if (ner/=0) v = 0.0
    !
    return
end subroutine lambertW

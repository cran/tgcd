subroutine calcFct(x, v) 
    implicit none
    real(kind=8), intent(in):: x
    real(kind=8), intent(out):: v
    !
    ! Local variables.
    real(kind=8), parameter:: a0=0.250621, a1=2.334733,&
                              b0=1.681534, b1=3.330657
    !
    v = 1.0 - (a0+a1*x+x**2)/(b0+b1*x+x**2)
    return
end subroutine calcFct

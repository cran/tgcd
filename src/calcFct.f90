subroutine calcFct(x, v) 
    ! 2023.09.07.
    implicit none
    real(kind(1.0d0)), intent(in):: x
    real(kind(1.0d0)), intent(out):: v
    !
    ! Local variables.
    real(kind(1.0d0)), parameter:: a0=0.250621, a1=2.334733,&
                                   b0=1.681534, b1=3.330657
    !
    v = 1.0 - (a0+a1*x+x**2)/(b0+b1*x+x**2)
    return
end subroutine calcFct

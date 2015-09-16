subroutine qeOTOR(nt, vect, y0, Nn, Ah,& 
                  An, ff, ae, hr, vecy, info)
!--------------------------------------------------------------------------------
! Subroutine qeOTOR is used for solving QE approximation equation 
! for the one-trap-one recombination center (OTOR) model.
!--------------------------------------------------------------------------------
!        nt, input:: integer, the number of data points.
!  vect(nt), input:: real values, temperature values.
!        y0, input:: real value, starting value of n.
!        Nn, input:: total concentration of traps in the crystal.
!        Ah, input:: real value, probability of electron recombine with holes.
!        An, input:: real value, probability of electron retrapping.
!        ff, input:: real value, the frequency factor.
!        ae, input:: real value, the activation energy.
!        hr, input:: real value, the heating rate.
! vecy(nt), output:: real values, y values at various temperatures.
!     info, output:: integer, error message, info=2 means a successful work.
!--------------------------------------------------------------------------------
! Author:: Peng Jun, 2015.09.13.
!--------------------------------------------------------------------------------
! Dependence:: subroutine dlsoda.
!--------------------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER nt, info
    DOUBLE PRECISION vect(nt), vecy(nt), y0,& 
                     Nn, Ah, An, ff, ae, hr
    !
    INTEGER NEQ(1), ITOL, ITASK, ISTATE, IOPT, & 
            JT, LRW, LIW, IWORK(25), i
    DOUBLE PRECISION RWORK(22+5*16), Y(1), T, TOUT, & 
                     RTOL(1), ATOL(1), kbz, bv
    EXTERNAL FUN, JAC
    !
    kbz = 8.617385e-05
    !  
    NEQ = 1
    ITOL = 1
    ITASK = 1
    ISTATE = 1
    IOPT = 0
    LRW = 22 + 5*16
    LIW = 25
    JT = 2
    RWORK = 0.0
    RWORK(6) = 1.0e+20
    IWORK = 0
    IWORK(1) = 1
    IWORK(2) = 1
    IWORK(6) = 10000
    IWORK(8) = 1
    IWORK(9) = 1
    !
    T = vect(1)
    Y = y0
    RTOL = 1.0e-10  
    ATOL = 1.0e-10 
    !
    bv = 0.0
    !
    DO i=1, nt
        TOUT = vect(i)
        CALL DLSODA (FUN, NEQ, Y, T, TOUT, ITOL, RTOL,& 
                     ATOL, ITASK, ISTATE, IOPT, RWORK,&
                     LRW, IWORK, LIW, JAC, JT, ff, ae,& 
                     Ah, An, Nn, hr, bv)
        info = ISTATE
        IF (ISTATE .LT. 0) RETURN   
        vecy(i) = Y(1)
    END DO
    !
    RETURN
END SUBROUTINE qeOTOR
!
SUBROUTINE FUN(NEQ, T, Y, YDOT, & 
               ff, ae, Ah, An, Nn, hr, bv)
    IMPLICIT NONE
    INTEGER NEQ
    DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), & 
                     kbz, ff, ae, Ah, An, Nn, hr, bv
    kbz = 8.617385e-05 
    YDOT(1) = -ff*(Y(1))**2*exp(-ae/kbz/T)*&
               Ah/(Ah*Y(1)+(Nn-Y(1))*An)/hr + bv*bv        
    RETURN
END SUBROUTINE FUN
!
SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
    IMPLICIT NONE
    INTEGER NEQ, ML, MU, NROWPD
    DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ)
    RETURN
END SUBROUTINE JAC

subroutine tgcd_drive(xd,yd,nd,pars,n2,fmin,lower,upper,nstart,&
                      mdt,mwt,mr,alw,kkf,ggt,tp,bg,tlsig3,tlsig4,suminfo,message)
!----------------------------------------------------------------------------------
! Subroutine tgcd_drive() is used for thermoluminescence 
! glow curve deconvolution according to various kinetic models
! using the Levenberg–Marquardt algorithm.
!----------------------------------------------------------------------------------
!                 xd(nd):: input, real values, observation X.
!                 yd(nd):: input, real vlaues, observations Y.
!                     nd:: input, integer, number of points.
!               pars(n2):: input/output, paraneters.
!                     n2:: input, integer, number of pars (<=39, 52).
!                   fmin:: output, real value, minimized objective.
!              lower(n2):: input, real values, lower bounds.
!              upper(n2):: input, real values, upper bounds.
!                 nstart:: input, integer, number of trials.
!                    mdt:: input, real value, allowed minimum 
!                          distance between Tm of glow peaks.
!                    mwt:: input, real value, maximum total 
!                          half-width for glow peaks.
!                     mr:: input, real value, allowed minimum
!                          resolution of glow peaks.
!                 alw(3):: input, integer values, whether the thresholds
!                          of mdt, mwt, mr will the applied, set the value
!                          to 1 to applied these thresholds.
!                    kkf:: input, real value (lie between 0 and 1).
!                    ggt:: input, integer (1, 2, or 3), 
!                          type of random initialization.
!                     tp:: input, integer, type of kinetic model,
!                          1=first-order (type 1),
!                          2=first-order (type 2),
!                          3=second-order,
!                          4=general-order (type 1),
!                          5=general-order (type 2),
!                          6=general-order (type 3),
!                          7=LW function,
!                          8=mix-order (type 1),
!                          9=mix-order (type 2),
!                          10=mix-order (type 3),
!                          11=weibull function,
!                          12=logistic asymmetric function;
!                          13=LW function (both branches).
!                     bg:: input, integer, subtract background or not,
!                          0=no subtraction, 1=subtraction.
! tlsig3(nd, (n2-3)/3+1):: output, real values, optimized TL signal values.
! tlsig3(nd, (n2-3)/4+1):: output, real values, optimized TL signal values.
!             suminfo(5):: output, integer values, a summary of error information.
!                message:: output, integer, error message: 0=success, 1=failure.
!----------------------------------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!----------------------------------------------------------------------------------
! Dependence:: subroutine tgcd_frt; 
!              subroutine tgcd_nonfrt.
!----------------------------------------------------------------------------------
    implicit none
    integer, intent(in):: nd, n2, nstart, alw(3), ggt, tp, bg
    real(kind(1.0d0)), intent(in):: xd(nd), yd(nd), mdt, mwt, mr, kkf
    real(kind(1.0d0)), intent(in):: lower(n2), upper(n2)
    real(kind(1.0d0)), intent(inout):: pars(n2)
    real(kind(1.0d0)), intent(out):: fmin, tlsig3(nd,(n2-3)/3+1), tlsig4(nd,(n2-3)/4+1)
    integer, intent(out):: suminfo(5), message
    !
    if (tp==1 .or. tp==2 .or. tp==3 .or. tp==11 .or. tp==12) then
        !
        call tgcd_frt(xd,yd,nd,pars,n2,fmin,lower,upper,nstart,&
                      mdt,mwt,mr,alw,kkf,ggt,tp,bg,tlsig3,suminfo,message)
        tlsig4 = -99.0
        !
    else if (tp==4 .or. tp==5 .or. tp==6 .or. tp==7 .or. tp==8 .or. tp==9 .or. tp==10 .or. tp==13) then
        !
        call tgcd_nonfrt(xd,yd,nd,pars,n2,fmin,lower,upper,nstart,&
                         mdt,mwt,mr,alw,kkf,ggt,tp,bg,tlsig4,suminfo,message)
        tlsig3 = -99.0
        !
    end if
    !
    return
end subroutine tgcd_drive

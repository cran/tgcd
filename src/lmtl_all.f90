subroutine lmtl_all(xd,yd,nd,pars,n2,fmin,&
                    message,lower,upper,tp,bg)
!----------------------------------------------------
! Subroutine lmtl_all() is used for fitting a TL glow 
! curve using the Levenberg-Marquardt algorithm using
! various kinetic models.
!----------------------------------------------------
!    xd(nd):: input, real values, observation X.
!    yd(nd):: input, real vlaues, observations Y.
!        nd:: input, integer, number of points.
!  pars(n2):: input/output, paraneters.
!        n2:: input, integer, number of pars (<=39, 52).
!      fmin:: output, real value, minimized objective.
!   message:: output, integer, error message,
!             0=success, 1=fail.
! lower(n2):: input, real values, lower bounds.
! upper(n2):: input, real values, upper bounds.
!        tp:: input, integer, type of kinetic model,
!             1=first-order (type 1),
!             2=first-order (type 2),
!             3=second-order,
!             4=general-order (type 1),
!             5=general-order (type 2),
!             6=general-order (type 3),
!             7=LW function,
!             8=mix-order (type 1),
!             9=mix-order (type 2),
!             10=mix-order (type 3),
!             11=weibull function,
!             12=logistic asymmetric equation,
!             13=LW function (both branches).
!        bg:: input, integer, subtract background or not,
!             0=no subtraction, 1=subtraction.
!------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!------------------------------------------------------
! Dependence:: subroutine lmdif1; 
!              subroutine tgcfunc_frt1;
!              subroutine tgcfunc_frt2;
!              subroutine tgcfunc_frt3;
!              subroutine tgcfunc_gnr1;
!              subroutine tgcfunc_gnr2;
!              subroutine tgcfunc_gnr3;
!              subroutine tgcfunc_lw;
!              subroutine tgcfunc_mix1;
!              subroutine tgcfunc_mix2;
!              subroutine tgcfunc_mix3;
!              subroutine tgcfunc_pdf1;
!              subroutine tgcfunc_pdf2;
!              subroutine tgcfunc_lw1. 
!------------------------------------------------------
    implicit none
    integer, intent(in):: nd, n2, tp, bg
    real(kind(1.0d0)), intent(in):: xd(nd), yd(nd)
    real(kind(1.0d0)), intent(in):: lower(n2), upper(n2)
    real(kind(1.0d0)), intent(inout):: pars(n2)
    real(kind(1.0d0)), intent(out):: fmin
    integer, intent(out):: message
    !
    integer:: info
    real(kind(1.0d0)), parameter:: tol=1.0e-07
    real(kind(1.0d0)):: fvec(nd)
    !
    external:: tgcfunc_frt1
    external:: tgcfunc_frt2
    external:: tgcfunc_frt3
    external:: tgcfunc_gnr1
    external:: tgcfunc_gnr2
    external:: tgcfunc_gnr3
    external:: tgcfunc_lw
    external:: tgcfunc_mix1
    external:: tgcfunc_mix2
    external:: tgcfunc_mix3
    external:: tgcfunc_pdf1
    external:: tgcfunc_pdf2
    external:: tgcfunc_lw1
    !
    fmin = -99.0
    !
    if (tp==1) then
        !
        call lmdif1(tgcfunc_frt1,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==2) then
        !
        call lmdif1(tgcfunc_frt2,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==3) then
        !
        call lmdif1(tgcfunc_frt3,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==4) then
        !
        call lmdif1(tgcfunc_gnr1,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==5) then
        !
        call lmdif1(tgcfunc_gnr2,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==6) then
        !
        call lmdif1(tgcfunc_gnr3,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==7) then
        !
        call lmdif1(tgcfunc_lw,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==8) then
        !
        call lmdif1(tgcfunc_mix1,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==9) then
        !
        call lmdif1(tgcfunc_mix2,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==10) then
        !
        call lmdif1(tgcfunc_mix3,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==11) then
        !
        call lmdif1(tgcfunc_pdf1,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==12) then
        !
        call lmdif1(tgcfunc_pdf2,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    else if (tp==13) then
        !
        call lmdif1(tgcfunc_lw1,nd,n2,pars,fvec,tol,info,xd,yd,lower,upper,bg)
        !
    end if
    !
    if (info==1 .or. info==2 .or. info==3) then
        message = 0
    else 
        message = 1
        return
    end if
    !
    fmin = sum(fvec**2)
    !
    return
end subroutine lmtl_all

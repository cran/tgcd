subroutine tgcd_nonfrt(xd,yd,nd,pars,n2,fmin,lower,upper,nstart,&
                       mdt,mwt,mr,alw,kkf,ggt,tpnonf,bg,tlsig,suminfo,message)
!-----------------------------------------------------------------------------------
! Subroutine tgcd_nonfrt() is used for thermoluminescence  
! glow curve deconvolution according to non firsr-order function 
! using the Levenberg–Marquardt algorithm.
!-----------------------------------------------------------------------------------
!               xd(nd):: input, real values, observation X.
!               yd(nd):: input, real vlaues, observations Y.
!                   nd:: input, integer, number of points.
!             pars(n2):: input/output, paraneters.
!                   n2:: input, integer, number of pars (<=52+3).
!                 fmin:: output, real value, minimized objective.
!            lower(n2):: input, real values, lower bounds.
!            upper(n2):: input, real values, upper bounds.
!               nstart:: input, integer, number of trials.
!                  mdt:: input, real value, allowed minimum 
!                        distance between Tm of glow peaks.
!                  mwt:: input, real value, allowed maximum total 
!                         half-width for glow peaks.
!                   mr:: input, real value, allowed minimum
!                        resolution of glow peaks.
!               alw(3):: input, integer values, whether the thresholds
!                        of mdt, mwt, mr will the applied, set the value
!                        to 1 to applied these thresholds.
!                  kkf:: input, real value (lie between 0 and 1).
!                  ggt:: input, integer (1 or 2), 
!                        type of random initialization.
!               tpnonf:: input, integer, type of non first-order model,
!                        4=general-order (type 1),
!                        5=general-order (type 2),
!                        6=general-order (type 3),
!                        7=the Lambert’s W,
!                        8=mix-order (type 1),
!                        9=mix-order (type 2),
!                        10=mix-order (type 3).
!                        13=the Lambert’s W (both branches),
!                   bg:: input, integer, subtract background or not,
!                        0=no subtraction, 1=subtraction.
! tlsig(nd,(n2-3)/4+1):: output, real values, optimized TL signal values.
!           suminfo(5):: output, integer values, a summary of error information.
!              message:: output, integer, error message, 0=success, 1=failure.
!--------------------------------------------------------------------------------
! Author:: Peng Jun, 2023.08.28.
!--------------------------------------------------------------------------------
! Dependence:: subroutine lmtl_all; 
!              subroutine hpSort;
!              subroutine calcShape;
!              subroutine calcMaty_gnr1;
!              subroutine calcMaty_gnr2;
!              subroutine calcMaty_gnr3;
!              subroutine calcMaty_lw;
!              subroutine calcMaty_mix1;
!              subroutine calcMaty_mix2;
!              subroutine calcMaty_mix3;
!              subroutine calcMaty_lw1.
!--------------------------------------------------------------------------------
    implicit none
    integer(kind=4), intent(in):: nd, n2, nstart, alw(3), ggt, tpnonf, bg
    real   (kind=8), intent(in):: xd(nd), yd(nd), mdt, mwt, mr, kkf
    real   (kind=8), intent(in):: lower(n2), upper(n2)
    real   (kind=8), intent(inout):: pars(n2)
    real   (kind=8), intent(out):: fmin, tlsig(nd,(n2-3)/4+1)
    integer(kind=4), intent(out):: suminfo(5), message
    ! Local variables.
    real   (kind=8):: ran(n2), ranpars(n2), pars0(n2), ranfmin, minfmin, unifa(n2), unifb(n2),&
                      orderTemp((n2-3)/4), mindist, maxwidth, minresol, resolvec((n2-3)/4-1),&
                      maty(nd,(n2-3)/4+1), matsp((n2-3)/4,7)
    integer(kind=4):: i, j, info, indx((n2-3)/4), flag, icy, n0, seed
    real   (kind=8), parameter:: ceof_a(9)=(/1.58, 1.766, 1.953, 2.141, 2.329,&
                                             2.519, 2.709,2.9,3.283/),&
                                 coef_x(9)=(/1.038,1.038,1.035,1.0325,1.0303,&
                                             1.0284,1.0267,1.0252,1.0226/)
    !
    seed = 123456789
    !
    n0 = n2 - 3
    !
    fmin = -99.0
    tlsig = -99.0
    suminfo = 0
    message = 1
    !
    pars0 = pars
    !
    ! Do optimization for the first time.
    ranpars = pars
    call lmtl_all(xd,yd,nd,ranpars,n2,ranfmin,info,lower,upper,tpnonf,bg)
    if (info/=0) then
        suminfo(1) = suminfo(1) + 1
        goto 999
    end if
    !
    !
    ! Calculate minimum distance between peaks if alw(1)=1.
    if (alw(1)==1) then 
        !
        orderTemp = ranpars(2*n0/4+1:3*n0/4)
        call hpSort(orderTemp, n0/4, indx) 
        mindist = minval(orderTemp(2:n0/4)-orderTemp(1:n0/4-1))
        if (mindist<mdt) then 
            suminfo(2) = suminfo(2) + 1
            goto 999
        end if
        !
    end if
    !
    !
    if (tpnonf==4) then
        !
        call calcMaty_gnr1(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpnonf==5) then
        !
        call calcMaty_gnr2(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpnonf==6) then
        !
        call calcMaty_gnr3(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpnonf==7) then
        !
        call calcMaty_lw(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpnonf==8) then
        !
        call calcMaty_mix1(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpnonf==9) then
        !
        call calcMaty_mix2(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpnonf==10) then
        !
        call calcMaty_mix3(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpnonf==13) then
        !
        call calcMaty_lw1(nd,n2,ranpars,xd,maty,bg)
        !
    end if
    !
    !
    ! Calculate shape parameters if alw(2)=1 or alw(3)=1.
    if (alw(2)==1 .or. (alw(3)==1 .and. n0/4>1)) then
        !
        call calcShape(nd,n0/4,xd,maty(:,1:n0/4),matsp,flag)
        if (flag/=0) then 
            suminfo(3) = suminfo(3) + 1
            goto 999
        end if
        !
    end if
    !
    !
    ! Calculate maximum total half-width of peaks if alw(2)=1.
    if (alw(2)==1) then
        !
        maxwidth = maxval(matsp(:,6))
        if (maxwidth>mwt) then
            suminfo(4) = suminfo(4) + 1
            goto 999
        end if
        !
    end if
    !
    !
    ! Calculate minimum resolution between glow peaks if alw(3)=1 and there are more than one glow peak.
    if (alw(3)==1 .and. n0/4>1) then
        !
        orderTemp = matsp(:,3)
        call hpSort(orderTemp, n0/4, indx) 
        matsp = matsp(indx,:)
        !
        do j=1, n0/4-1
           resolvec(j) = (matsp(j+1,3)-matsp(j,3)) / (matsp(j,5)+matsp(j+1,4))
        end do
        !
        minresol = minval(resolvec)
        if (minresol<mr) then
            suminfo(5) = suminfo(5) + 1
            goto 999
        end if
        !
    end if
    !
    !!!if (mindist>=mdt .and. flag==0 .and. maxwidth<=mwt) then
    !!!end if
    message = 0
    pars = ranpars
    fmin = ranfmin
    tlsig = maty
    !
    999 continue
    !!!if (nstart==1)  return
    !
    if (message/=0)  then
        minfmin = 1.0e+20
    else 
        minfmin = fmin          
    end if
    !
    ! INTENS.
    unifa(1:n0/4) = pars0(1:n0/4)*(1.0-kkf)
    unifb(1:n0/4) = pars0(1:n0/4)*(1.0+kkf)
    !
    ! ENERGY.
    unifa(n0/4+1:2*n0/4) = pars0(n0/4+1:2*n0/4)*(1.0-kkf)
    unifb(n0/4+1:2*n0/4) = pars0(n0/4+1:2*n0/4)*(1.0+kkf)
    !
    ! TEMPER.
    unifa(2*n0/4+1:3*n0/4) = pars0(2*n0/4+1:3*n0/4)*(1.0-kkf)
    unifb(2*n0/4+1:3*n0/4) = pars0(2*n0/4+1:3*n0/4)*(1.0+kkf)
    !
    ! bValue, rValue, or aValue.
    unifa(3*n0/4+1:n0) = pars0(3*n0/4+1:n0)*(1.0-kkf)
    unifb(3*n0/4+1:n0) = pars0(3*n0/4+1:n0)*(1.0+kkf)
    !
    ! Background values.
    unifa(n0+1:n0+3) = pars0(n0+1:n0+3)*(1.0-kkf)
    unifb(n0+1:n0+3) = pars0(n0+1:n0+3)*(1.0+kkf)
    !
    ! Reset bValue, rValue, or aValue if ggt==2.
    if (ggt==2) then
        !
        if (tpnonf==4 .or. tpnonf==5 .or. tpnonf==6) then 
            ! bValue.
            unifa(3*n0/4+1:n0) = 1.01
            unifb(3*n0/4+1:n0) = 1.99
        else if (tpnonf==7 .or. tpnonf==8 .or. tpnonf==9 .or. tpnonf==10 .or. tpnonf==13) then
            ! rValue or aValue.
            unifa(3*n0/4+1:n0) = 0.01
            unifb(3*n0/4+1:n0) = 0.99
            !
            if (tpnonf==13)  unifb(3*n0/4+1:n0) = 1.5
        end if        
        !
    end if
    !
    !
    ! Try-and-error.
    TryError: do i=1,  nstart
        call r8vec_uniform_01(n2, seed, ran)
        ranpars = unifa + ran*(unifb-unifa)
        !
        ! Reset ENERGY if ggt==2.
        if (ggt==2) then
            icy = mod(i,9)
            if (icy==0) icy = 9
            !
            ranpars(n0/4+1:2*n0/4) = ceof_a(icy)*ranpars(2*n0/4+1:3*n0/4)**(coef_x(icy))*1.0e-03
        end if
        !
        call lmtl_all(xd,yd,nd,ranpars,n2,ranfmin,info,lower,upper,tpnonf,bg)
        if (info/=0) then
            suminfo(1) = suminfo(1) + 1
            cycle TryError
        end if
        !
        !
        ! Calculate minimum distance between glow peaks if alw(1)=1.
        if (alw(1)==1) then
            !
            orderTemp = ranpars(2*n0/4+1:3*n0/4)
            call hpSort(orderTemp, n0/4, indx) 
            mindist = minval(orderTemp(2:n0/4)-orderTemp(1:n0/4-1))
            if (mindist<mdt) then
                suminfo(2) = suminfo(2) + 1
                cycle TryError
            end if
            !
        end if
        !
        !
        if (tpnonf==4) then
            !
            call calcMaty_gnr1(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpnonf==5) then
            !
            call calcMaty_gnr2(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpnonf==6) then
            !
            call calcMaty_gnr3(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpnonf==7) then
            !
            call calcMaty_lw(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpnonf==8) then
            !
            call calcMaty_mix1(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpnonf==9) then
            !
            call calcMaty_mix2(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpnonf==10) then
            !
            call calcMaty_mix3(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpnonf==13) then
            !
            call calcMaty_lw1(nd,n2,ranpars,xd,maty,bg)
            !
        end if
        !
        !
        ! Calculate shape parameters if alw(2)=1 or alw(3)=1.
        if (alw(2)==1 .or. (alw(3)==1 .and. n0/4>1)) then
            !
            call calcShape(nd,n0/4,xd,maty(:,1:n0/4),matsp,flag)
            if (flag/=0) then
                suminfo(3) = suminfo(3) + 1
                cycle TryError
            end if
            !
        end if
        !
        !
        ! Calculate maximum total half-width of glow peaks if alw(2)=1.
        if (alw(2)==1) then
            !
            maxwidth = maxval(matsp(:,6))
            if (maxwidth>mwt) then 
                suminfo(4) = suminfo(4) + 1
                cycle TryError
            end if
            !
        end if
        !
        !
        ! Calculate minimum resolution between glow peaks if alw(3)=1 and there are more than one glow peak.
        if (alw(3)==1 .and. n0/4>1) then
            !
            orderTemp = matsp(:,3)
            call hpSort(orderTemp, n0/4, indx) 
            matsp = matsp(indx,:)
            !
            do j=1, n0/4-1
                resolvec(j) = (matsp(j+1,3)-matsp(j,3)) / (matsp(j,5)+matsp(j+1,4))
            end do
            !
            minresol = minval(resolvec)
            if (minresol<mr) then
                suminfo(5) = suminfo(5) + 1
                cycle TryError
            end if
            !
        end if
        !
        !!!if (ranfmin<minfmin .and. mindist>=mdt .and. &
        !!!flag==0 .and. maxwidth<=mwt)  then
        !!!end if
        if (ranfmin<minfmin) then
            message = 0
            pars = ranpars
            fmin = ranfmin
            tlsig = maty
            minfmin = ranfmin
        end if
        !
    end do TryError
    !
    return
end subroutine tgcd_nonfrt

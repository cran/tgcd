subroutine tgcd_frt(xd,yd,nd,pars,n2,fmin,lower,upper,nstart,&
                    mdt,mwt,mr,alw,kkf,ggt,tpf,bg,tlsig,suminfo,message)
!----------------------------------------------------------------------------------
! Subroutine tgcd_frt() is used for thermoluminescence glow curve
! deconvolution according to first-order (or second order) kinetics 
! using the Levenberg–Marquardt algorithm.
!----------------------------------------------------------------------------------
!               xd(nd):: input, real values, observation X.
!               yd(nd):: input, real vlaues, observations Y.
!                   nd:: input, integer, number of points.
!             pars(n2):: input/output, paraneters.
!                   n2:: input, integer, number of pars (<=39+3).
!                 fmin:: output, real value, minimized objective.
!            lower(n2):: input, real values, lower bounds.
!            upper(n2):: input, real values, upper bounds.
!               nstart:: input, integer, number of trials.
!                  mdt:: input, real value, allowed minimum 
!                        distance between Tm of glow peaks.
!                  mwt:: input, real value, allowed maximum total 
!                        half-width for glow peaks.
!                   mr:: input, real value, allowed minimum
!                        resolution of glow peaks.
!               alw(3):: input, integer values, whether the thresholds
!                        of mdt, mwt, mr will the applied, set the value
!                        to 1 to applied these thresholds.
!                  kkf:: input, real value (lie between 0 and 1).
!                  ggt:: input, integer (1 or 2), 
!                        type of random initialization.
!                  tpf:: input, integer, type of kinetic model,
!                        1=first-order (type 1),
!                        2=first-order (type 2),
!                        3=second-order,
!                        11=weibull function,
!                        12=logistic asymmetric equation.
!                   bg:: input, integer, subtract background or not,
!                        0=no subtraction, 1=subtraction.
! tlsig(nd,(n2-3)/3+1):: output, real values, optimized TL signal values.
!           suminfo(5):: output, integer values, a summary of error information.
!              message:: output, integer, error message: 0=success, 1=failure.                  
!---------------------------------------------------------------------------------
! Author:: Peng Jun, 2023.09.07.
!---------------------------------------------------------------------------------
! Dependence:: subroutine lmtl_all; 
!              subroutine hpSort;
!              subroutine r8vec_uniform_01;
!              subrontine calcShape;
!              subroutine calcMaty_frt1;
!              subroutine calcMaty_frt2;
!              subroutine calcMaty_frt3;
!              subroutine calcMaty_pdf1;
!              subroutine calcMaty_pdf2.
!---------------------------------------------------------------------------------
    implicit none
    integer, intent(in):: nd, n2, nstart, alw(3), ggt, tpf, bg
    real(kind(1.0d0)), intent(in):: xd(nd), yd(nd), mdt, mwt, mr, kkf
    real(kind(1.0d0)), intent(in):: lower(n2), upper(n2)
    real(kind(1.0d0)), intent(inout):: pars(n2)
    real(kind(1.0d0)), intent(out):: fmin, tlsig(nd,(n2-3)/3+1)
    integer, intent(out):: suminfo(5), message
    ! Local variables.
    real(kind(1.0d0)):: ran(n2), ranpars(n2), pars0(n2), ranfmin, minfmin, unifa(n2), unifb(n2),&
                        orderTemp((n2-3)/3), mindist, maxwidth, minresol, resolvec((n2-3)/3-1),&
                        maty(nd,(n2-3)/3+1), matsp((n2-3)/3,7)
    integer:: i, j, info, indx((n2-3)/3), flag, icy, n0, seed
    real(kind(1.0d0)), parameter:: ceof_a(9)=(/1.58, 1.766, 1.953, 2.141, 2.329,&
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
    call lmtl_all(xd,yd,nd,ranpars,n2,ranfmin,info,lower,upper,tpf,bg)
    if (info/=0) then 
        suminfo(1) = suminfo(1) + 1
        goto 999
    end if
    !
    !
    ! Calculate minimum distancce between glow peaks if alw(1)=1.
    if (alw(1)==1) then 
        !
        orderTemp = pars(2*n0/3+1:n0)
        call hpSort(orderTemp, n0/3, indx) 
        mindist = minval(orderTemp(2:n0/3)-orderTemp(1:n0/3-1))
        if (mindist<mdt) then 
            suminfo(2) = suminfo(2) + 1
            goto 999
        end if
        !
    end if
    !
    !
    if (tpf==1) then
        !
        call calcMaty_frt1(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpf==2) then
        !
        call calcMaty_frt2(nd,n2,ranpars,xd,maty,bg)
            !
    else if (tpf==3) then
        !
        call calcMaty_frt3(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpf==11) then
        !
        call calcMaty_pdf1(nd,n2,ranpars,xd,maty,bg)
        !
    else if (tpf==12) then
        !
        call calcMaty_pdf2(nd,n2,ranpars,xd,maty,bg)
        !
    end if
    !
    !
    ! Calculate shape parameters if alw(2)=1 or alw(3)=1.
    if (alw(2)==1 .or. (alw(3)==1 .and. n0/3>1)) then
        !
        call calcShape(nd,n0/3,xd,maty(:,1:n0/3),matsp,flag)
        if (flag/=0) then 
            suminfo(3) = suminfo(3) + 1
            goto 999
        end if
        !
    end if 
    !
    !
    ! Calculate maximum tota half-width of peaks if alw(2)=1.
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
    if (alw(3)==1 .and. n0/3>1) then
        !
        orderTemp = matsp(:,3)
        call hpSort(orderTemp, n0/3, indx) 
        matsp = matsp(indx,:)
        !
        do j=1, n0/3-1
            resolvec(j) = (matsp(j+1,3) - matsp(j,3)) / (matsp(j,5)+matsp(j+1,4))
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
    !!!if (mindist>=mdt .and. flag==0 .and. maxwidth<=mwt .and. minresol>=mr) then
    !!!end if
    message = 0
    pars = ranpars
    fmin = ranfmin
    tlsig = maty
    !
    999 continue
    !!!if (nstart==1) return
    !
    if (message/=0) then
        minfmin = 1.0e+20
    else 
        minfmin = fmin
    end if 
    !
    ! INTENS.
    unifa(1:n0/3) = pars0(1:n0/3)*(1.0-kkf)
    unifb(1:n0/3) = pars0(1:n0/3)*(1.0+kkf)
    !
    ! ENERGY.
    unifa(n0/3+1:2*n0/3) = pars0(n0/3+1:2*n0/3)*(1.0-kkf)
    unifb(n0/3+1:2*n0/3) = pars0(n0/3+1:2*n0/3)*(1.0+kkf)
    !
    ! TEMPER.
    unifa(2*n0/3+1:n0) = pars0(2*n0/3+1:n0)*(1.0-kkf)
    unifb(2*n0/3+1:n0) = pars0(2*n0/3+1:n0)*(1.0+kkf)
    !
    ! Background values.
    unifa(n0+1:n0+3) = pars0(n0+1:n0+3)*(1.0-kkf)
    unifb(n0+1:n0+3) = pars0(n0+1:n0+3)*(1.0+kkf)
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
            ranpars(n0/3+1:2*n0/3) = ceof_a(icy)*ranpars(2*n0/3+1:n0)**(coef_x(icy))*1.0e-03
        end if
        !
        call lmtl_all(xd,yd,nd,ranpars,n2,ranfmin,info,lower,upper,tpf,bg)
        if (info/=0) then
            suminfo(1) = suminfo(1) + 1
            cycle TryError
        end if
        !
        ! 
        ! Calculate minimum distance between glow peaks if alw(1)=1.
        if (alw(1)==1) then
            !   
            orderTemp = ranpars(2*n0/3+1:n0)
            call hpSort(orderTemp, n0/3, indx) 
            mindist = minval(orderTemp(2:n0/3)-orderTemp(1:n0/3-1))
            if (mindist<mdt) then
                suminfo(2) = suminfo(2) + 1
                cycle TryError
            end if
            !
        end if
        !
        !
        if (tpf==1) then
            !
            call calcMaty_frt1(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpf==2) then
            !
            call calcMaty_frt2(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpf==3) then
            !
            call calcMaty_frt3(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpf==11) then
            !
            call calcMaty_pdf1(nd,n2,ranpars,xd,maty,bg)
            !
        else if (tpf==12) then
            !
            call calcMaty_pdf2(nd,n2,ranpars,xd,maty,bg)
            !
        end if
        !
        !
        ! Calculate shape parameters if alw(2)=1 or alw(3)=1.
        if (alw(2)==1 .or. (alw(3)==1 .and. n0/3>1)) then
            !
            call calcShape(nd,n0/3,xd,maty(:,1:n0/3),matsp,flag)
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
        if (alw(3)==1 .and. n0/3>1) then
            !
            orderTemp = matsp(:,3)
            call hpSort(orderTemp, n0/3, indx) 
            matsp = matsp(indx,:)
            !
            do j=1, n0/3-1
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
        !!! flag==0 .and. maxwidth<=mwt .and. minresol>=mr)  then
        !!! end if
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
end subroutine tgcd_frt

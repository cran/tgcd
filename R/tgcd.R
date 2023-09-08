###
tgcd <-
function(Sigdata, npeak, model="g1", subBG=FALSE, pickp="d2", 
         pickb="d0", nstart=60, kkf=0.03, mdt=NULL, mwt=NULL, 
         mr=NULL, edit.inis=TRUE, inisPAR=NULL, inisBG=NULL, 
         hr=NULL, hwd=NULL, pod=NULL, plot=TRUE, outfile=NULL)  {
    UseMethod("tgcd")
} #
### 2023.09.07.
tgcd.default <- 
function(Sigdata, npeak, model="g1", subBG=FALSE, pickp="d2", 
         pickb="d0", nstart=60, kkf=0.03, mdt=NULL, mwt=NULL, 
         mr=NULL, edit.inis=TRUE, inisPAR=NULL, inisBG=NULL,  
         hr=NULL, hwd=NULL, pod=NULL, plot=TRUE, outfile=NULL)  {
        ### Stop if not.
        stopifnot(ncol(Sigdata)==2L, 
                  ###all(Sigdata[,1L,drop=TRUE]>0),
                  ###all(Sigdata[,2L,drop=TRUE]>=0),
                  length(npeak)==1L, is.numeric(npeak), 
                  npeak %in% seq(13L), 3L*npeak<nrow(Sigdata), 
                  length(model)==1L, model %in% c("f1","f2","f3","s1","s2","g1","g2","g3","lw","m1","m2","m3","wo"),
                  length(subBG)==1L, is.logical(subBG),
                  length(pickp)==1L, pickp %in% c("d0","d01","d1","d2","d3","d4"),
                  length(pickb)==1L, pickb %in% c("d0","d01"),
                  length(nstart)==1L, is.numeric(nstart), nstart>1L, nstart<=10000L,
                  length(kkf)==1L, is.numeric(kkf), kkf>0, kkf<1.0,
                  is.null(mdt) || is.numeric(mdt),
                  is.null(mwt) || is.numeric(mwt),
                  is.null(mr)  || is.numeric(mr),
                  length(edit.inis)==1L, is.logical(edit.inis),
                  is.null(inisPAR) || is.matrix(inisPAR),
                  is.null(inisBG) || is.numeric(inisBG),
                  is.null(hr) || is.numeric(hr),
                  is.null(hwd) || is.numeric(hwd),
                  is.null(pod) || is.numeric(pod),
                  is.logical(plot), length(plot)==1L,
                  is.null(outfile) || is.character(outfile))
        ###
        ###
        if (!is.null(mdt)) {
            if (length(mdt)!=1L) stop("Error: mdt should be an one-element vector!")
            if (mdt<0) stop("Error: mdt should not be smaller than 0!")
            if(npeak==1L) cat("Note: mdt will not be used because npeak=1!\n")
        } # end if.
        ###
        if (!is.null(mwt)) {
            if (length(mwt)!=1L) stop("Error: mwt should be an one-element vector!")
            if (mwt<=0) stop("Error: mwt should exceed 0!")
            if(npeak==1L) cat("Note: mwt will not be used because npeak=1!\n")
        } # end if.
        ###
        if (!is.null(mr)) {
            if(length(mr)!=1L) stop("Error: mr should be an one-element vector!")
            if (mr<0) stop("Error: mr should not be smaller than 0!")
            if(npeak==1L) cat("Note: mr will not be used because npeak=1!\n")
        } # end if.
        ###
        if (!is.null(inisPAR))  {
            if(dim(inisPAR)[1L]!=npeak) stop("Error: incorrect dimensions of inisPAR!")
            ###
            if(model %in% c("f1","f2","f3","s1","s2") && dim(inisPAR)[2L]!=3L) stop("Error: incorrect dimensions of inisPAR!")
            ### 
            if(model %in% c("g1","g2","g3","lw","m1","m2","m3","wo") && dim(inisPAR)[2L]!=4L) stop("Error: incorrect dimensions of inisPAR!")
            ###
            if (any(inisPAR<=0)) stop("Error: all elements in inisPAR should be larger than 0!")
        } # end if.
        ###
        if(!is.null(inisBG)) {
           if (length(inisBG)!=3L) stop("Error: inisBG should be a 3-element numeric vector!")
           ###
           if (any(inisBG<=0)) stop("Error: all elements in inisBG should be larger than 0!")
        } # end if.
        ###
        if (!is.null(hr)) {
            if(length(hr)!=1L) stop("Error: hr should be a one-element vector!")
            ###
            if(hr<=0) stop("Error: hr should be larger than zero!")
        } # end if.
        ###
        if (!is.null(hwd)) {
            if (length(hwd)!=1L) stop("Error: hwd should be a one-element vector!")
            if (hwd<1L) stop("Error: hwd should not be smaller than 1!")
        } # end if.
        ###
        if (!is.null(pod)) {
            if (length(pod)!=1L) stop("Error: pod should be a one-element vector!")
            if (pod<1L) stop("Error: pod should not be smaller than 1!")
            if (pod>2.0*hwd) stop("Error: pod should not exceed 2*hwd!")
        } # end if.
        ### 
        if (!is.null(outfile)) {
            if (length(outfile)!=1L) stop("Error: outfile should be an one-element vector!")
        } # end if.
        ###
        ###
        ### Temperature and signal values.
        temp <- as.numeric(Sigdata[,1L,drop=TRUE])
        if (min(temp)<273.0) cat("Warning: the minimum temperature is smaller than 273 K!\n")
        ###
        signal <- as.numeric(Sigdata[,2L,drop=TRUE])
        ###
        ###
        if ( is.null(inisPAR) ||  (subBG==TRUE && is.null(inisBG)) || plot==TRUE) {
            ###
            opar <- par("mfrow", "mar")
            on.exit(par(opar))
            ###
        } # end if.
        ###
        ###
        ### If it is NULL the argument inisPAR, starting parameters need be initlized manually.
        if(is.null(inisPAR))  {
            ###
            par(mar=c(4.5,4.5,4,2)+0.1)
            ###
            if (pickp %in% c("d0","d01")) {
                ###
                abzero <- which(signal>(.Machine$double.eps)^0.3)
                ###
                plot(temp[abzero], signal[abzero], type="p", pch=21, cex=1.5, bg="grey50", xlab="Temperature (K)", 
                     log=ifelse(pickp=="d0","","y"), ylab="TL intensity (counts)", main=paste("Click the mouse to select ", 
                     npeak, " peak maxima:", sep=""), cex.lab=1.3)
                ###
                drv <- 0L
                if (is.null(hwd)) hwd <- 3L*(drv+2L)
                ###
                if (is.null(pod)) pod <- 4L
                ###
                dy_signal <- try(savgol(y=signal,drv=drv,hwd=hwd,pod=pod),silent=TRUE)
                if (inherits(dy_signal, what="try-error")==FALSE) points(temp[abzero], dy_signal[abzero], type="l", col="skyblue3", lwd=5)
                ###
            } else if (pickp %in% c("d1","d2","d3","d4")) {
                ###
                if(pickp=="d1") {
                    drv <- 1L
                    YLAB <- "First-order derivative of TL glow curve"
                } else if(pickp=="d2") {
                    drv <- 2L
                    YLAB <- "Second-order derivative of TL glow curve"
                } else if(pickp=="d3") {
                    drv <- 3L
                    YLAB <- "Third-order derivative of TL glow curve"
                } else if(pickp=="d4") {
                    drv <- 4L
                    YLAB <- "Fourth-order derivative of TL glow curve"
                } # end if.
                ###
                if (is.null(hwd)) hwd <- 3L*(drv+2L)
                ###
                if (is.null(pod)) pod <- 4L
                ###
                dy_signal <- try(savgol(y=signal,drv=drv,hwd=hwd,pod=pod),silent=TRUE)
                if (inherits(dy_signal, what="try-error")==TRUE) stop("Error: failed in derivative calculation!")
                ###
                plot(temp, dy_signal, type="l", col="skyblue3", lwd=5.0, xlab="Temperature (K)", ylab=YLAB, 
                     main=paste("Click the mouse to select ", npeak, " peak locations:", sep=""), cex.lab=1.3)
                axis(side=2, col="skyblue3", lwd=3.9)
                ###
                scale_signal <- signal/max(signal)*max(dy_signal)
                points(temp, scale_signal, type="l", col="red", lwd=1.2, lty="dashed")
                ###
                abline(h=0.0, col="purple", lwd=3.0, lty="dashed")
                ###
            } # end if.
            ###
            grid(col="darkviolet", lwd=1.0)
            ###
            sldxy <- try(locator(n=npeak), silent=TRUE)
            if(inherits(sldxy, what="try-error")==TRUE) stop("Error: failed in manual initilization of kinetic parameters!")
            ###
            ###
            sldxy_index <- order(sldxy$x, decreasing=FALSE)
            sldxy$x <- sldxy$x[sldxy_index]
            ###
            #if (pickp %in% c("d0","d01")) {
                #sldxy$y <- sldxy$y[sldxy_index]
            #} else if (pickp %in% c("d1","d2","d3","d4")) {
                #linear_interp <- suppressWarnings(try(approx(x=temp, y=signal, xout=sldxy$x),silent=TRUE))
                #if (class(linear_interp)=="try-error" || !is.finite(linear_interp$y)) {
                    #stop("Error: failed in linear interpolation to obtain intensity (Im) by using temperature (Tm)!")
                #} # end if.
                #sldxy$y <- linear_interp$y
            #} # end if
            ###
            linear_interp <- suppressWarnings(try(approx(x=temp, y=signal, xout=sldxy$x),silent=TRUE))
            if (inherits(linear_interp, what="try-error")==TRUE || any(!is.finite(linear_interp$y))) {
                stop("Error: failed in linear interpolation to obtain intensity (Im) by using temperature (Tm)!")
            } # end if.
            sldxy$y <- linear_interp$y
            ###
            par(mar=c(5,4,4,2)+0.1)
            ###
            ###TmIm_vec <- rbind(sldxy$x, sldxy$y)
            ###rownames(TmIm_vec) <- c("x","y")
            ###colnames(TmIm_vec) <- paste(seq(npeak),"th-Peak",sep="")
            ###cat("\n")
            ###cat("The selected coordinates (Tm, Im) for glow peaks are:\n")
            ###cat("-------------------------------------------------------------------------------------\n")
            ###print(TmIm_vec)
            ###cat("-------------------------------------------------------------------------------------\n\n\n")
            ###
        } # end if.
        ###
        ### If subBG=TRUE and inisBG=NULL, starting parameters need be initlized manually.
        if (subBG==TRUE && is.null(inisBG)) {
            ###
            par(mar=c(4.5,4.5,4,2)+0.1)
            ###
            abzero <- which(signal>(.Machine$double.eps)^0.3)
            ###
            plot(temp[abzero], signal[abzero], type="p", pch=21, cex=1.1, bg="grey70", xlab="Temperature (K)", 
                 cex.lab=1.3, log=ifelse(pickb=="d0","","y"), ylab="TL intensity (counts)", 
                 main="Click the mouse to select 3 points for background initialization")
            ###
            drv <- 0L
            if (is.null(hwd)) hwd <- 3L*(drv+2L)
            ###
            if (is.null(pod)) pod <- 4L
            ###
            dy_signal <- try(savgol(y=signal,drv=drv,hwd=hwd,pod=pod),silent=TRUE)
            if (inherits(dy_signal, what="try-error")==FALSE) points(temp, dy_signal, type="l", col="skyblue3", lwd=3.0)
            ###
            sldxyBG <- try(locator(n=3L), silent=TRUE)
            if(inherits(sldxyBG, what="try-error")==TRUE) stop("Error: failed in automatical initilization of background parameters!")
            ###
            ###---------------------------------------------------------------------------
            bgFUNC <- function(p,x,y,X,Y) {
                p <- abs(p)
                v <- sum((p[1L]+p[2L]*exp(x/p[3L])-y)^2)
                if (!is.finite(v) || any(Y<p[1L]+p[2L]*exp(X/p[3L]))) { 
                    return(.Machine$double.xmax) } else { 
                    return(v) }
            } # end function bgFUNC.
            ###---------------------------------------------------------------------------
            ###
            minVAL <- .Machine$double.xmax
            ntrial <- 600L
            ###
            for (i in seq(ntrial)) {
                ###
                ba <- exp(runif(n=1L, min=log(1.0e-03),max=log(1.0e03)))
                bb <- exp(runif(n=1L, min=log(1.0e-03),max=log(1.0e03)))
                bc <- runif(n=1L, min=10.0,max=max(temp))
                ###
                p <- c(ba,bb,bc)
                ###
                NLM <- try(nlm(f=bgFUNC, p=p, x=sldxyBG$x, y=sldxyBG$y, X=temp, Y=signal), silent=TRUE)
                ###
                if (inherits(NLM, what="try-error")==FALSE && NLM$minimum<minVAL) {
                    minVAL <- NLM$minimum
                    BGpar <- abs(NLM$estimate)
                } # end if.
            } # end for.
            ###
            if (!exists("BGpar")) stop("Error: failed in automatical initilization of background parameters!")
            ###
            ba <- BGpar[1L]
            bb <- BGpar[2L]
            bc <- BGpar[3L]
            ###
            points(temp, ba+bb*exp(temp/bc), type="l", col="red", lwd=2.0, lty="dashed")
            points(temp, signal-(ba+bb*exp(temp/bc)), type="l", col="purple", lwd=2.0)
            ###
            Sys.sleep(3L)
            ###
            par(mar=c(5,4,4,2)+0.1)
        } # end if.
        ###
        ###
        minTEMPER <- min(temp)
        maxTEMPER <- max(temp)
        minINTENS <- min(signal)
        maxINTENS <- max(signal)
        ###
        ### Function used for editing initial parameters (inisPAR and inisBG) manually (interactively).
        setpars <-function(npeak, sldx, sldy)  {
            mat1 <- mat2 <- mat3 <- mat4 <- as.data.frame(matrix(nrow=npeak+1L, ncol=5L))
            ###
            ### Default TL growth peak intensity. 
            mat1[1L,] <-  c("Peak", "INTENS(min)", "INTENS(max)", "INTEN(ini)", "INTENS(fix)")
            mat1[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat1[-1L,2L] <-  round(minINTENS*0.8, 3L)
            mat1[-1L,3L] <-  round(maxINTENS*1.2, 3L) 
            mat1[-1L,4L] <- if (is.null(inisPAR)) {round(sldy, 3L)} else {round(inisPAR[,1L,drop=TRUE], 3L)} # end if.
            mat1[-1L,5L] <- FALSE
            ###
            ### Default activation energy.
            mat2[1L,]<- c("Peak", "ENERGY(min)", "ENERGY(max)", "ENERGY(ini)", "ENERGY(fix)")
            mat2[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat2[-1L,2L] <- 0.3 
            mat2[-1L,3L] <- 4.0
            mat2[-1L,4L] <- if (is.null(inisPAR)) {round(runif(n=npeak,min=1.0, max=2.0), 3L)} else {round(inisPAR[,2L,drop=TRUE], 3L)} # end if.
            mat2[-1L,5L] <- FALSE
            ###
            ### Default temperature at the peak maximum.
            mat3[1L,] <- c("Peak", "TEMPER(min)", "TEMPER(max)", "TEMPER(ini)", "TEMPER(fix)")
            mat3[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
            mat3[-1L,2L] <-  round(minTEMPER, 3L)
            mat3[-1L,3L] <-  round(maxTEMPER, 3L)
            mat3[-1L,4L] <- if (is.null(inisPAR)) {round(sldx, 3L)} else {round(inisPAR[,3L,drop=TRUE], 3L)} # end if.
            mat3[-1L,5L] <- FALSE
            ###
            ### Default bValue, rValue, or aValue for a glow peak.
            if (model %in% c("g1","g2","g3"))  {
                mat4[1L,] <- c("Peak", "bValue(min)", "bValue(max)", "bValue(ini)", "bValue(fix)")
                mat4[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
                mat4[-1L,2L] <- 1.0
                mat4[-1L,3L] <- 2.0
                mat4[-1L,4L] <- if (is.null(inisPAR)) {round(runif(n=npeak,min=1.01, max=1.99), 3L)} else {round(inisPAR[,4L,drop=TRUE], 3L)} # end if.
                mat4[-1L,5L] <- FALSE
            } else if(model %in% c("lw","wo")) {
                mat4[1L,] <- c("Peak", "rValue(min)", "rValue(max)", "rValue(ini)","rValue(fix)")
                mat4[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
                mat4[-1L,2L] <- 1.0e-16
                mat4[-1L,3L] <-  ifelse(model=="wo",0.93,2.0)
                mat4[-1L,4L] <- if (is.null(inisPAR)) { if(model=="wo") round(runif(n=npeak,min=0.01, max=0.5), 3L) else 
                                round(runif(n=npeak,min=0.01, max=1.9), 3L) } else { inisPAR[,4L,drop=TRUE] } # end if.
                    
                mat4[-1L,5L] <- FALSE
            } else if (model %in% c("m1","m2","m3")) {
                mat4[1L,] <- c("Peak", "aValue(min)", "aValue(max)", "aValue(ini)","aValue(fix)")
                mat4[-1L,1L] <- paste(seq(npeak),"th-Peak", sep="")
                mat4[-1L,2L] <- 1.0e-7
                mat4[-1L,3L] <-  0.9999999
                mat4[-1L,4L] <- if (is.null(inisPAR)) {round(runif(n=npeak,min=0.01, max=0.99), 3L)} else {inisPAR[,4L,drop=TRUE]} # end if.
                mat4[-1L,5L] <- FALSE
            } # end if.
            ###
            ### Default parameters for inisBG.
            if (subBG==TRUE) {
                mat5 <- as.data.frame(matrix(nrow=4L, ncol=5L))
                ###
                mat5[1L,] <- c("BG", "BG(min)", "BG(max)", "BG(ini)","BG(fix)")
                mat5[-1L,1L] <- c("A","B","C")
                mat5[-1L,2L] <- c(1.0e-09, 1.0e-09, 1.0e-09)
                mat5[-1L,3L] <- c(1.0e+09, 1.0e+09, 1.0e+09)
                mat5[-1L,4L] <- if (is.null(inisBG)) { BGpar } else { inisBG } # end if.
                mat5[-1L,5L] <- FALSE
            } # end if.
            ###
            if (model %in% c("f1","f2","f3","s1","s2")) {
                ###
                if (subBG==FALSE) mat <- rbind(mat1, rep("    ", 5L), mat2, rep("    ", 5L), mat3)
                if (subBG==TRUE)  mat <- rbind(mat1, rep("    ", 5L), mat2, rep("    ", 5L), mat3, rep("    ", 5L), mat5)
                ###
            } else if (model %in% c("g1","g2","g3","lw","m1","m2","m3","wo")) {
                ###
                if (subBG==FALSE) mat <- rbind(mat1, rep("    ", 5L), mat2, rep("    ", 5L), mat3, rep("    ", 5L), mat4)  
                if (subBG==TRUE)  mat <- rbind(mat1, rep("    ", 5L), mat2, rep("    ", 5L), mat3, rep("    ", 5L), mat4, rep("    ", 5L), mat5) 
                ###
            } # end if.    
            ###
            ###
            if (edit.inis==TRUE) {
                pars <- try(edit(name=mat), silent=TRUE)
                if (inherits(pars, what="try-error")==TRUE) stop("Error: incorrect modification of initial parameters!")
            } else {
                pars <- mat
            } # end if.
            ###
            return(pars)
       } # end function setpars.
       ###
       ###
       ### cat("Set parameter constraints:\n")
       pars <- setpars(npeak=npeak, sldx=sldxy$x, sldy=sldxy$y)
       indx <- seq(from=2L, to=npeak+1L, by=1L)
       ###
       ################################################################################
       ################################################################################
       ### 1. TL growth peak intensity. Check non-finite vlaues.
       intensity1 <- as.numeric(pars[indx,2L,drop=TRUE])
       whichloc <- which(!is.finite(intensity1))
       if (length(whichloc)>=1L) stop("Error: non-finite lower bound of INTENS")
       ###
       intensity2 <- as.numeric(pars[indx,3L,drop=TRUE])
       whichloc <- which(!is.finite(intensity2))
       if (length(whichloc)>=1L) stop("Error: non-finite upper bound of INTENS")
       ###
       intensity3 <- as.numeric(pars[indx,4L,drop=TRUE])
       whichloc <- which(!is.finite(intensity3))
       if (length(whichloc)>=1L) stop("Error: non-finite initial of INTENS")
       ###
       ###
       ### Check bounds.
       whichloc <- which(intensity3<intensity1 | intensity3>intensity2)
       if (length(whichloc)>=1L) stop("Error: unbounded initial of INTENS")
       ### 
       ### Check logical values.
       fix_intensity <- pars[indx,5L,drop=TRUE]
       if (!all(fix_intensity %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of INTENS!")
       } # end if.
       fix_intensity <- as.logical(fix_intensity)
       intensity1[fix_intensity==TRUE] <- intensity3[fix_intensity==TRUE]
       intensity2[fix_intensity==TRUE] <- intensity3[fix_intensity==TRUE]
       ###
       ###
       ### 2. Activation energy. Check non-finite values.
       energy1 <- as.numeric(pars[indx+(npeak+2L),2L,drop=TRUE])
       whichloc <- which(!is.finite(energy1)) 
       if (length(whichloc)>=1L) stop("Error: non-finite lower bound of ENERGY!")
       ###
       energy2 <- as.numeric(pars[indx+(npeak+2L),3L,drop=TRUE])
       whichloc <- which(!is.finite(energy2)) 
       if (length(whichloc)>=1L) stop("Error: non-finite upper bound of ENERGY!")
       ###
       energy3 <- as.numeric(pars[indx+(npeak+2L),4L,drop=TRUE])
       whichloc <- which(!is.finite(energy3)) 
       if (length(whichloc)>=1L) stop("Error: non-finite initial of ENERGY!")
       ###
       ### Check bounds.
       whichloc <- which(energy3<energy1 | energy3>energy2)
       if (length(whichloc)>=1L) stop("Error: unbounded initial of ENERGY")
       ###
       ### Check logical values.
       fix_energy <- pars[indx+(npeak+2L),5L,drop=TRUE]
       if (!all(fix_energy %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of ENERGY!")
       } # end if.
       fix_energy <- as.logical(fix_energy)
       energy1[fix_energy==TRUE] <- energy3[fix_energy==TRUE]
       energy2[fix_energy==TRUE] <- energy3[fix_energy==TRUE]
       ###
       ###
       ### 3. Temperature at the peak maximum. Check non-finite values.
       temperature1 <- as.numeric(pars[indx+2L*(npeak+2L),2L,drop=TRUE])
       whichloc <- which(!is.finite(temperature1))
       if (length(whichloc)>=1L) stop("Error: non-finite lower bound of TEMPER")
       ###
       temperature2 <- as.numeric(pars[indx+2L*(npeak+2L),3L,drop=TRUE])
       whichloc <- which(!is.finite(temperature2))
       if (length(whichloc)>=1L) stop("Error: non-finite upper bound of TEMPER")
       ###
       temperature3 <- as.numeric(pars[indx+2L*(npeak+2L),4L,drop=TRUE])
       whichloc <- which(!is.finite(temperature3))
       if (length(whichloc)>=1L) stop("Error: non-finite initial of TEMPER")
       ###
       ### Check bounds.
       whichloc <- which(temperature3<temperature1 | temperature3>temperature2)
       if (length(whichloc)>=1L) stop("Error: unbounded initial of TEMPER")
       ###
       ### Check logical values.
       fix_temperature <- pars[indx+2L*(npeak+2L),5L,drop=TRUE]
       if (!all(fix_temperature %in% c("TRUE", "FALSE", "T", "F")))  {
           stop("Error: non-logical variable in the 5th column of TEMPER!")
       } # end if.
       fix_temperature <- as.logical(fix_temperature)
       temperature1[fix_temperature==TRUE] <- temperature3[fix_temperature==TRUE]
       temperature2[fix_temperature==TRUE] <- temperature3[fix_temperature==TRUE]
       ###
       ###
       ### 4. bValue, rValue, or aValue for the the glow peak. 
       if (model %in% c("g1","g2","g3","lw","m1","m2","m3","wo")) {
           ###
           label <- if(model %in% c("g1","g2","g3")) {
               "bValue" 
           } else if (model %in% c("lw","wo")) {
               "rValue" 
           } else if (model %in% c("m1","m2","m3")) {
               "aValue" 
           } # end if.
           ###
           ### Check non-finite values.
           bValue1 <- as.numeric(pars[indx+3L*(npeak+2L),2L,drop=TRUE])
           whichloc <- which(!is.finite(bValue1))
           if (length(whichloc)>=1L) stop(paste("Error: non-finite lower bound of ", label, "!", sep=""))
           ###
           bValue2 <- as.numeric(pars[indx+3L*(npeak+2L),3L,drop=TRUE])
           whichloc <- which(!is.finite(bValue2))
           if (length(whichloc)>=1L) stop(paste("Error: non-finite upper bound of ", label, "!", sep=""))
           ###
           bValue3 <- as.numeric(pars[indx+3L*(npeak+2L),4L,drop=TRUE])
           whichloc <- which(!is.finite(bValue3))
           if (length(whichloc)>=1L) stop(paste("Error: non-finite initial of ", label, "!", sep=""))
           ###
           ### Check bounds.
           whichloc <- which(bValue3<bValue1 | bValue3>bValue2)
           if (length(whichloc)>=1L) stop(paste("Error: unbounded initial of ", label, "!", sep=""))
           ###
           ### Check logical values.
           fix_bValue <- pars[indx+3L*(npeak+2L),5L,drop=TRUE]
           if (!all(fix_bValue %in% c("TRUE", "FALSE", "T", "F")))  {
               stop(paste("Error: non-logical variable in the 5th column of ", label, "!", sep=""))
           } # end if.
           fix_bValue <- as.logical(fix_bValue)
           bValue1[fix_bValue==TRUE] <- bValue3[fix_bValue==TRUE]
           bValue2[fix_bValue==TRUE] <- bValue3[fix_bValue==TRUE]
           ###
       } # end if.
       ###
       ### 5. Check parameters A, B, C, and D in the background expression.
       if (subBG==TRUE) {
           ###
           if (model %in% c("f1","f2","f3","s1","s2"))  {
               LLL <- npeak+1L+2L*(npeak+2L)+seq(from=3L,to=5L,by=1L)
           } else if (model %in% c("g1","g2","g3","lw","m1","m2","m3","wo")) {
               LLL <- npeak+1L+3L*(npeak+2L)+seq(from=3L,to=5L,by=1L)
           } # end if.
           ###
           ### Check non-fninte values.
           BGpars1 <- as.numeric(pars[LLL,2L,drop=TRUE])
           whichloc <- which(!is.finite(BGpars1))
           if (length(whichloc)>=1L) stop("Error: non-finite lower bound of BGpars")
           ###
           BGpars2 <- as.numeric(pars[LLL,3L,drop=TRUE])
           whichloc <- which(!is.finite(BGpars2))
           if (length(whichloc)>=1L) stop("Error: non-finite upper bound of BGpars")
           ###
           BGpars3 <- as.numeric(pars[LLL,4L,drop=TRUE])
           whichloc <- which(!is.finite(BGpars3))
           if (length(whichloc)>=1L) stop("Error: non-finite initials of BGpars")
           ###
           ### Check bounds.
           whichloc <- which(BGpars3<BGpars1 | BGpars3>BGpars2)
           if (length(whichloc)>=1L) stop("Error: unbounded initials of BGpars")
           ###
           ### Check logical values.
           fix_BGpars <- pars[LLL,5L,drop=TRUE]
           if (!all(fix_BGpars %in% c("TRUE", "FALSE", "T", "F")))  {
               stop("Error: non-logical variable in the 5th column of BGpars!")
           } # end if.
           fix_BGpars <- as.logical(fix_BGpars)
           BGpars1[fix_BGpars==TRUE] <- BGpars3[fix_BGpars==TRUE]
           BGpars2[fix_BGpars==TRUE] <- BGpars3[fix_BGpars==TRUE]
           ###   
       } else {
           BGpars1 <- BGpars2 <- BGpars3 <- rep(0.0, 3L)
       } # end if.
       ##############################################################################
       ##############################################################################
       ###
       ###
       nd <- length(temp)
       n2 <- ifelse(model %in% c("f1","f2","f3","s1","s2"), 3L*npeak+3L, 4L*npeak+3L)
       fmin <- 0.0
       ###
       if (model %in% c("f1","f2","f3","s1","s2")) {
           lower <- c(intensity1, energy1, temperature1, BGpars1)
           upper <- c(intensity2, energy2, temperature2, BGpars2)
           pars <-  c(intensity3, energy3, temperature3, BGpars3)
       } else if (model %in% c("g1","g2","g3","lw","m1","m2","m3","wo")) {
           lower <- c(intensity1, energy1, temperature1, bValue1, BGpars1)
           upper <- c(intensity2, energy2, temperature2, bValue2, BGpars2)
           pars <-  c(intensity3, energy3, temperature3, bValue3, BGpars3)
       } # end if.
       ###
       nstart <- nstart - 1L
       ###
       alw <- rep(0L, 3L)
       if (!is.null(mdt) && npeak>1L) alw[1L] <- 1L
       if (!is.null(mwt) && npeak>1L) alw[2L] <- 1L
       if (!is.null(mr) && npeak>1L)  alw[3L] <- 1L
       ###
       ggt <- ifelse(!is.null(inisPAR), 1L, 2L)
       ###
       tp <- if (model=="f1") {
           1L
       } else if (model=="f2") {
           2L
       } else if (model=="f3") {
           11L
       } else if (model=="s1") {
           3L
       } else if (model=="s2") {
           12L
       } else if (model=="g1") {
           4L
       } else if (model=="g2") {
           5L
       } else if (model=="g3") {
           6L
       } else if (model=="lw") {
           13L
       } else if (model=="m1") {
           8L
       } else if (model=="m2") {
           9L
       } else if (model=="m3") {
           10L
       } else if (model=="wo") {
           7L
       } # end if.
       ###
       bg <- ifelse(subBG==FALSE, 0L, 1L)
       ###
       tlsig3 <- matrix(0.0, nrow=nd, ncol=(n2-3L)/3L+1L)
       tlsig4 <- matrix(0.0, nrow=nd, ncol=(n2-3L)/4L+1L)
       ###
       suminfo <- rep(0L, 5L)
       message <- 0L
       ###
       cat("\n")
       cat("Thermoluminescence glow curve deconvolution (tgcd) is in progress, please wait, ...\n")
       cat("-------------------------------------------------------------------------------------\n\n")
       ###
       res <- .Fortran("tgcd_drive", as.double(temp), as.double(signal),
                       as.integer(nd), pars=as.double(pars), as.integer(n2), 
                       fmin=as.double(fmin), as.double(lower), as.double(upper), 
                       as.integer(nstart), as.double(mdt), as.double(mwt), as.double(mr), 
                       as.integer(alw), as.double(kkf), as.integer(ggt), as.integer(tp), 
                       as.integer(bg), tlsig3=as.double(tlsig3), tlsig4=as.double(tlsig4), 
                       suminfo=as.integer(suminfo), message=as.integer(message), 
                       PACKAGE="tgcd")
       ###
       cat("\n")
       cat("A summary of information generated from the trial-and-error protocol:\n")
       cat("--------------------------------------------------------------------------------\n")
       ###
       cat("Failed in the Levenberg-Marquardt algorithm:   ",res$suminfo[1L],"\n")
       ###
       if (!is.null(mdt) && npeak>1L) {
           cat("Minimum distance between glow peaks is smaller than mdt=",mdt,":   ",res$suminfo[2L],"\n",sep="")
       } # end if.
       ###
       if ((!is.null(mwt) && npeak>1L) || (!is.null(mr) && npeak>1L)) {
           cat("Failed in calculation of shape parameters of glow peaks:   ",res$suminfo[3L],"\n")
       } # end if.
       ###
       if (!is.null(mwt) && npeak>1L) {
           cat("Maximum total half-width of glow peaks is greater than mwt=",mwt,":   ",res$suminfo[4L],"\n",sep="")
       } # end if.
       ###
       if (!is.null(mr) && npeak>1L) {
           cat("Minimum resolution of glow peaks is smaller than mr=",mr,":   ",res$suminfo[5L],"\n", sep="")
       } # end if.
       ###
       cat("Succeeded in the trial-and-error protocol (nstart=",nstart+1L,"):   ",nstart+1L-sum(res$suminfo),"\n",sep="")
       cat("--------------------------------------------------------------------------------\n\n")
       ###
       if (res$message!=0L) stop("Error: fail in glow curve deconvolution!")
       ###
       pars <- matrix(res$pars[1L:(n2-3L)], ncol=ifelse(model %in% c("f1","f2","f3","s1","s2"),3L,4L))
       index <- order(pars[,3L,drop=TRUE], decreasing=FALSE)
       pars <- pars[index,,drop=FALSE]
       colnames(pars) <- if (model %in% c("f1", "f2","f3","s1","s2")) {
           c("INTENS(Im)", "ENERGY(E)", "TEMPER(Tm)")
       } else if (model %in% c("g1","g2","g3")) {
           c("INTENS(Im)", "ENERGY(E)", "TEMPER(Tm)", "bValue(b)")
       } else if (model %in% c("lw","wo")) {
           c("INTENS(Im)", "ENERGY(E)", "TEMPER(Tm)", "rValue(r)")
       } else if (model %in% c("m1","m2","m3")) {
           c("INTENS(Im)", "ENERGY(E)", "TEMPER(Tm)", "aValue(a)")
       } # end if.
       rownames(pars) <- paste(seq(npeak),"th-Peak",sep="")
       ###
       if (subBG==TRUE) {
           BGpars <- res$pars[(n2-2L):n2]
           names(BGpars) <- c("A","B","C")
       } else {
           BGpars <- NULL
       } # end if.
       ###
       kbz <- 8.617385e-5
       ###
       ###
       if (model %in% c("f1","f2","f3","s1","s2")) {
           CompSig <- matrix(res$tlsig3, ncol=npeak+1L)
       } else if (model %in% c("g1","g2","g3","lw","m1","m2","m3","wo")) {
           CompSig <- matrix(res$tlsig4, ncol=npeak+1L)
       } # end if.
       CompSig[,1L:npeak] <- CompSig[,index,drop=FALSE]
       ###
       rowsumSig <- rowSums(CompSig)
       residuals <- signal - rowsumSig
       SSR <- sum(residuals^2)
       RCS <- SSR/(nd-n2)
       R2 <- (cor(x=signal,y=rowsumSig,method="pearson"))^2
       FOM <- res$fmin/sum(rowsumSig)*100.0
       ###
       ###
       ###################################################################################################
       ### R interfaces for fortran subroutine wrightOmega(), calcei(), and calcAm().
       ###################################################################################################
       wrightOmega <- function(Z) {
           W <- double(1L)
           res <- .Fortran("wrightOmega", as.double(Z), W=as.double(W), PACKAGE="tgcd")
           return(res$W)
       } # end function wrightOmega.
       ###
       calcEi <- function(x) {
           int <- as.integer(1L)
           val <- double(1L)
           res <- .Fortran("calcei", as.double(x), val=as.double(val), as.integer(int), PACKAGE="tgcd")
           return(res$val)
       } # end function calcEi.
       ###
       calcAm <- function(ax, bx, alpha, maxt, engy) {
           Am <- double(1L)
           res <- .Fortran("calcAm", as.double(ax), as.double(bx), as.double(alpha), as.double(maxt),
                           as.double(engy), Am=as.double(Am), as.double(fmin), PACKAGE="tgcd")
           return(res$Am)
       } # end function calcAm.
       ###
       lambertW <- function(xx) {
           v <- double(1L)
           ner <- integer(1L)
           res <- .Fortran("lambertW", as.double(xx), v=as.double(v), ner=as.integer(ner), PACKAGE="tgcd")
           return(res$v)
       } # end function lambertW.
       ###
       ##################################################################################################
       ###
       ### Calculate shape parameters for glow peaks.
       calShape <- function(y, x)  {
           ny <- length(y)
           maxloc <- which.max(y)
           hmaxval <- max(y)/2.0
           Tm <- x[maxloc]
           T1 <- suppressWarnings(try(approx(x=y[1L:maxloc], y=x[1L:maxloc], xout=hmaxval)$y, silent=TRUE))
           T2 <- suppressWarnings(try(approx(x=y[maxloc:ny], y=x[maxloc:ny], xout=hmaxval)$y, silent=TRUE))
           ###
           if (inherits(T1, what="try-error")==FALSE) { d1 <- Tm-T1 } else { T1 <- d1 <- NA } # end if.
           ###
           if (inherits(T2, what="try-error")==FALSE) { d2 <- T2-Tm } else { T2 <- d2 <- NA } # end if.  
           ###          
           thw <- T2-T1
           sf <- d2/thw
           ###         
           return(c("T1"=T1, "T2"=T2, "Tm"=Tm, "d1"=d1, "d2"=d2, "thw"=thw, "sf"=sf))
       } # end function calShape.
       ###
       ###
       sp <- t(apply(CompSig[,-(npeak+1L),drop=FALSE], MARGIN=2L, calShape, temp))
       rownames(sp) <- paste(seq(npeak),"th-Peak",sep="")
       ###
       maybewrongTmidx <- which(abs(pars[,3L]-sp[,3L])>=10.0)
       if (length(maybewrongTmidx)>=1L) {
           cat("Warning: significant inconsistency between Tm in [pars] and Tm in [sp] for glow peak(s)",maybewrongTmidx,"!\n")
       } # end if.
       ###
       if (npeak>1L) {
           resolvec <- name_resolvec <- vector(length=npeak-1L) 
           for (i in seq(npeak-1L)) {
               resolvec[i] <- (sp[i+1L,3L]-sp[i,3L])/(sp[i,5L]+sp[i+1L,4L])
               name_resolvec[i] <- paste("Peak",i,i+1L,sep="")
           } # end for.
           names(resolvec) <- name_resolvec
       } else {
           resolvec <- NULL
       } # end if.
       ###
       ### Calculate frequency factors for glow peaks.
       if (!is.null(hr))  {
           ###
           calff_first <- function(et)  {
               energy <- et[1L]
               temper <- et[2L]
               ff <- hr*energy/kbz/temper^2*exp(energy/kbz/temper)
               return(ff)          
           } # end function calff_first.
           ###
           calff_second <- function(et) {
               energy <- et[1L]
               temper <- et[2L]
               ff <- hr*energy/kbz/temper^2/(1.0+(2.0*kbz*temper/energy))*exp(energy/kbz/temper)
               return(ff)
           } # end function calff_second.
           ###
           calff_general <- function(et) {
               energy <- et[1L]
               temper <- et[2L]
               bv <- et[3L]
               ff <- hr*energy/kbz/temper^2/(1.0+(bv-1.0)*(2.0*kbz*temper/energy))*exp(energy/kbz/temper)   
               return(ff)
           } # end function calff_general.
           ###
           calff_wo <- function(et)  {
               energy <- et[1L]
               temper <- et[2L]
               rv <- et[3L]
               ###
               xi <- min(temp)
               eivi <- calcEi(-energy/kbz/xi)
               Feivi <- xi*exp(-energy/kbz/xi) + energy/kbz*eivi
               eiv <- calcEi(-energy/kbz/temper)
               ftem <- (temper*exp(-energy/kbz/temper) + energy/kbz*eiv) - Feivi
               z1m <- rv/(1.0-rv) - log((1.0-rv)/rv) + energy*exp(energy/kbz/temper)/
                      kbz/temper^2/(1.0-1.05*rv^1.26)*ftem
               wz1m <- wrightOmega(z1m)
               ###
               ff <- hr*energy/kbz/temper^2/exp(-energy/kbz/temper)/(1.0/(1.0-rv)*(1.0+2.0*wz1m)/(1.0+wz1m)^2)
               return(ff)             
           } # end function calff_wo.
           ###
           calff_mix1 <- function(et) {
               energy <- et[1L]
               temper <- et[2L]
               alpha <- et[3L]
               ###
               Am <- calcAm(ax=0.01, bx=10, alpha, temper, energy)
               Rm <- (Am+alpha)/(Am-alpha)
               ###
               ff <- hr*energy/kbz/temper^2*Rm*(alpha/(1.0-alpha))*exp(energy/kbz/temper) 
               return(ff)
           } # end function calff_mix1.
           ###
           calff_mix23 <- function(et) {
               energy <- et[1L]
               temper <- et[2L]
               alpha <- et[3L]
               ###
               Rm <- (1.0-alpha)*(1.0+0.2922*alpha-0.2783*alpha^2)
               ff <- hr*energy/kbz/temper^2*Rm*(alpha/(1.0-alpha))*exp(energy/kbz/temper) 
               return(ff)
           } # end function calff_mix23.
           ###
           calff_lw <- function(et)  {
               energy <- et[1L]
               temper <- et[2L]
               rv <- et[3L]
               ###
               xi <- min(temp)
               eivi <- calcEi(-energy/kbz/xi)
               Feivi <- xi*exp(-energy/kbz/xi) + energy/kbz*eivi
               eiv <- calcEi(-energy/kbz/temper)
               ftem <- (temper*exp(-energy/kbz/temper) + energy/kbz*eiv) - Feivi
               if (rv<1.0) {
                   z1m <- rv/(1.0-rv) - log((1.0-rv)/rv) + energy*exp(energy/kbz/temper)/
                          kbz/temper^2/(1.0-1.05*rv^1.26)*ftem
                   wz1m <- wrightOmega(z1m)
               } else {
                   z1m <- abs(rv/(1.0-rv)) + log(abs((1.0-rv)/rv)) + energy*exp(energy/kbz/temper)/
                          kbz/temper^2/(2.963-3.24*rv^(-0.74))*ftem
                   if (exp(-z1m) < .Machine$double.xmin) {
                       wz1m <- -z1m - log(z1m)
                   } else {
                       wz1m <- lambertW(-exp(-z1m))
                   } # end if.
               } # end if.
               ###
               ff <- hr*energy/kbz/temper^2/exp(-energy/kbz/temper)/(1.0/(1.0-rv)*(1.0+2.0*wz1m)/(1.0+wz1m)^2)
               return(ff)             
           } # end function calff_lw.
           ###
           ###
           if (model %in% c("f1","f2","f3")) {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff_first)
           } else if (model %in% c("s1","s2")) {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff_second)
           } else if (model %in% c("g1","g2","g3"))  {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff_general)
           } else if (model=="wo") {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff_wo)
           } else if (model=="m1") {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff_mix1)
           } else if (model %in% c("m2","m3")) {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff_mix23)
           } else if (model=="lw") {
               ff <- apply(pars[,-1L,drop=FALSE], MARGIN=1L, calff_lw)
           } # end if.
           ### 
       } else {
           ff <- NULL
       } # end if.
       ###
       ###
       ### Plot the results.
       if (plot==TRUE) {
           ###
           layout(cbind(c(rep(1,13), 2, rep(3,6)), c(rep(1,13), 2, rep(3,6))))
           ###
           ### The first plot.
           par(mar=c(0,5.1,3.1,1.1))
           lineCol <- c("deepskyblue", "orangered", "purple", "violetred", "yellowgreen", "lightblue", 
                        "goldenrod", "forestgreen", "blue", "plum", "tan", "violet", "grey40")
           plot(temp, signal, type="p", pch=21, bg="white", cex=1.0, ylab="TL intensity (counts)", 
                las=0, lab=c(7,7,9), xaxt="n", xaxs="r", yaxs="i", cex.lab=2.0, cex.axis=1.5, ylim=c(0.0, max(signal)*1.1))
           box(lwd=2.0)
           XaxisCentral <- median(axTicks(side=1L))
           ###
           for (i in seq(npeak)) {
               points(temp,CompSig[,i,drop=TRUE], type="l", lwd=2.0, col=lineCol[i])
           } # end for.
           ###
           points(temp, rowsumSig, type="l", lwd=2.0, col="black")
           ###
           if (subBG==FALSE) {
               legend(ifelse(temp[which.max(signal)]>XaxisCentral,"topleft","topright"),
                      legend=c("Fitted.Curve", paste(seq(npeak),"th-Peak",sep=""), paste("FOM=",round(FOM,2),"%",sep="")), 
                      col=c("black", lineCol[seq(npeak)], NA), pch=c(21, rep(NA,npeak),NA),
                      lty=c(rep("solid",npeak+1L),NA), yjust=2, ncol=1, cex=2.0, bty="o", 
                      lwd=2.0, pt.bg="white")
           } else {
               points(temp,CompSig[,npeak+1L,drop=TRUE], type="l", lwd=2, col="grey70")
               ###
               legend(ifelse(temp[which.max(signal)]>XaxisCentral,"topleft","topright"),
                      legend=c("Fitted.Curve", paste(seq(npeak),"th-Peak",sep=""),"Background",paste("FOM=",round(FOM,2),"%",sep="")), 
                      col=c("black", lineCol[seq(npeak)], "grey70", NA), pch=c(21, rep(NA,npeak+1L),NA),
                      lty=c(rep("solid",npeak+2L),NA), yjust=2, ncol=1, cex=2.0, bty="o", 
                      lwd=2.0, pt.bg="white")
           } # end if.
           ###
           ### The second plot.
           par(mar=c(0,5.1,0,1.1))
           plot(c(0,0), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
           ###
           ### The third plot.
           par(mar=c(5.1,5.1,0,1.1))
           plot(temp, residuals, type="p", xlab="Temperature (K)", ylab="Residuals",
                las=0, lab=c(7,7,9), xaxs="r", yaxs="i", pch=21, bg="grey", cex=0.8, cex.lab=2, 
                cex.axis=1.5)
           ### abline(h=0)
           box(lwd=2.0)
           ###
       } # end if.
       ###
       ###
       if (subBG==TRUE) {
           CompSig <- cbind(temp, signal, rowsumSig,  CompSig)
           colnames(CompSig) <- c("Temperature", "Obs.Signal", "Fit.Signal", paste("Peak.", seq(npeak), sep = ""), "Background")
       } else {
           CompSig <- cbind(temp, signal, rowsumSig,  CompSig[,-(npeak+1L)])
           colnames(CompSig) <- c("Temperature", "Obs.Signal", "Fit.Signal", paste("Peak.", seq(npeak), sep = "")) 
       } # end if.
       ###
       if (!is.null(outfile)) write.csv(CompSig, file=paste(outfile, ".csv", sep = ""))
       ###
       output <-list("comp.sig"=CompSig, "residuals"=residuals, "pars"=pars, "BGpars"=BGpars, 
                     "ff"=ff, "sp"=sp, "resolution"=resolvec, "SSR"=SSR, "RCS"=RCS, "R2"=R2, "FOM"=FOM)
       ###
       invisible(output) 
       ###                
} # end fucntion tgcd.
###

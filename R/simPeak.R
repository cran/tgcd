###
simPeak <- 
function(temps, tl0, ae, ff, hr=1.0, kmax=20,
         tol=1e-6, outfile=NULL, plot=TRUE) {
    UseMethod("simPeak")
} #
### 2015.06.18.
simPeak.default <- 
function(temps, tl0, ae, ff, hr=1.0, kmax=20,
         tol=1e-6, outfile=NULL, plot=TRUE) {
    ### Stop if not.
    stopifnot(is.numeric(temps), all(temps>0),
              length(tl0)==1L, tl0>0,
              length(ae)==1L, is.numeric(ae), ae>0,
              length(ff)==1L, is.numeric(ff), ff>0,
              length(hr)==1L, is.numeric(hr), hr>0,
              length(kmax)==1L, is.numeric(kmax), kmax>0,
              length(tol)==1L, is.numeric(tol), tol>0,
              is.null(outfile) || is.character(outfile),
              length(plot)==1L, is.logical(plot))
    ###
    if (!is.null(outfile)) {
        if (length(outfile)!=1L) 
            stop("Error: outfile should be an one-element vector!")
    } # end if.
    ###
    nt <- length(temps)
    vecy <- double(nt)
    ###
    res <- .Fortran("BulirschStoer", as.integer(nt), as.double(temps),
                    as.double(tl0), as.double(tol), as.integer(kmax),
                    as.double(ff), as.double(ae), as.double(hr), 
                    vecy=as.double(vecy), PACKAGE="tgcd")
     ###
     tout <- temps[-nt]
     yout <- res$vecy[seq(nt-1L)]-
             res$vecy[2L:nt]
     ###
     if (plot==TRUE)  {
         layout(cbind(c(1L, 1L, 1L, 2L), 
                      c(1L, 1L, 1L, 2L)))
         par(mar = c(0, 5.1, 3.1, 1.1))
         ###
         plot(tout, yout, type="l", lwd=5, col="skyblue3",
              ylab="TL Intensity", las=0, lab=c(7,7,9), 
              xaxt="n", xaxs="r", yaxs="r", cex.lab=2*par("cex"))
         box(lwd=2L)
         ###
         XaxisCentral <- median(axTicks(side=1L))
         ###
         legend(ifelse(tout[which.max(yout)] > XaxisCentral, "topleft", 
                "topright"), legend=c(paste("Initial Concentration: ",
                format(tl0,digits=3,scientific=TRUE)," (1/cm^3)",sep=""),
                paste("Activation Energy: ",round(ae,2L)," (eV)",sep=""),
                paste("Frequency Factor: ",format(ff,digits=3, scientific=TRUE),
                " (1/s)",sep=""), paste("Heating Rate: ", round(hr,2L),
                " (K/s)",sep="")), yjust=2, ncol=1, cex=2*par("cex"), 
                bty="n",  pt.bg="white")
         ###
         par(mar = c(5.1, 5.1, 0, 1.1))
         plot(tout, res$vecy[seq(nt-1L)], type="l", xlab = "Temperature(K)", 
              ylab="Concentration", las=0, lab=c(7,7,9), xaxs="r", yaxs="i",
              ylim=c(-0.1*tl0, 1.1*tl0), col="grey", lwd=5, cex.lab=2*par("cex"))
         box(lwd=2L) 
         par(mar=c(5,4,4,2)+0.1)
         layout(1L)
     } # end if.
     ###
     if (!is.null(outfile))  {
         PeakSig <- cbind(tout, yout)
         colnames(PeakSig) <- c("temps","tl")
         write.csv(PeakSig, file=paste(outfile, ".csv", sep=""))
     } else {
         return(list(temps=tout, tl=yout))
     } # end if.
} # end function simPeak.
###

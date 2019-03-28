###
savgol <- 
function(y,drv,hwd=3*(drv+2),pod=4) {
    UseMethod("savgol")
} # end function savgol.
### 2019.03.26.
savgol.default <- 
function(y,drv,hwd=3*(drv+2),pod=4) {
    stopifnot(is.numeric(y), 
              is.numeric(drv), length(drv)==1L, drv %in% (0L:6L), 
              is.numeric(hwd), length(hwd)==1L, hwd>=1L,
              is.numeric(pod), length(pod)==1L, pod>=1L,
              drv<=pod, 2.0*hwd>=pod)
    ###
    nl <- nr <- hwd
    ###
    ld <- drv
    ###
    m <- pod
    ###
    n1 <- length(y)
    ###
    flag <- 0L
    ###
    res <- .Fortran("savgol_filter", as.integer(nl),as.integer(nr),
                    as.integer(ld),as.integer(m),as.integer(n1), 
                    y=as.double(y),flag=as.integer(flag),PACKAGE="tgcd")
    ###
    if (res$flag==1L) stop("Savitsky-Golay algorithm failed!")
    ###
    return(res$y)
} # end function savgol.default.

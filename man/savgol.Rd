\name{savgol}
\alias{savgol}
\alias{savgol.default}

\title{
Apply a Savitzky-Golay algorithm to smooth thermoluminescence glow curves
}

\description{
Smooth thermoluminescence glow curves with a Savitzky-Golay smoothing filter and calculate derivatives.
}
\usage{
savgol(y, drv, hwd = 3 * (drv + 2), pod = 4)
}


\arguments{
  \item{y}{\link{numeric}(\bold{required}): the data to be filtered}

  \item{drv}{\link{integer}(\bold{required}): the order of the derivative to be calculated}

  \item{hwd}{\link{integer}(with default): half width of the segement used for filter}

  \item{pod}{\link{integer}(with default): order of the polynomial used for filter}
}

\details{
The Savitzky-Golay smoothing algorithm is particularly good at preserving lineshape while removing high frequency squiggles (Press et al., 1986). 
The procedure can be used to calculate derivatives of thermoluminescence data to identify the location of glow peaks.
}

\value{
The filtered signal, which has the same length as \code{y}.
}

\references{
Press WH, Teukolsky SA, Vetterling WT, Flannery BP, 1986. Numberic recipes in Fortran 77, the Art of Scientific Computing, second edition. 
}

\seealso{
\link{tgcd}
}

\examples{
 library(tgcd)
 data(Refglow)

 x <- Refglow$x009[,1]
 y <- Refglow$x009[,2]
 y0 <- savgol(y, drv=0)
 dy <- savgol(y, drv=1)

 plot(x, y, type="p", pch=21, bg="black")
 points(x, y0, type="l", col="blue", lwd=2)

 plot(x, dy, type="l", col="blue", lwd=2)
 abline(h=0, lty="dashed", col="red", lwd=2)

}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Savitzky-Golay smoothing}

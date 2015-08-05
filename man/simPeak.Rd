\name{simPeak}
\alias{simPeak}
\alias{simPeak.default}
\title{
Thermoluminescence glow curve simulation
}
\description{
Simulating a synthetic glow peak according to the Randall-Wilkins model.
}
\usage{
simPeak(temps, tl0, ae, ff, hr = 1, kmax = 20, 
        tol = 1e-06, outfile = NULL, plot = TRUE)
}

\arguments{
  \item{temps}{\link{vector}(\bold{required}): temperature values (K) where the values of the thermoluminescence 
intensity will be computed. It needs to be sorted increasingly. A vector of temperature values may be generated 
using the internal function \link{seq}}

  \item{tl0}{\link{numeric}(\bold{required}): the concentration of trapped electrons (1/cm^3)}

  \item{ae}{\link{numeric}(\bold{required}): the activation energy (eV)}

  \item{ff}{\link{numeric}(\bold{required}): the frequency factor (1/s)}

  \item{hr}{\link{numeric}(with default): the linear heating rate (K/s)}

  \item{kmax}{\link{integer}(with default): allowed maximum number of steps in the Richardson extrapolation}

  \item{tol}{\link{numeric}(with default): relative tolerance in the Richardson extrapolation}

  \item{outfile}{\link{character}(optional): if specified, simulated signal values of a glow peak will 
be written to a file named \code{"outfile"} in \code{CSV} format and saved to the current work directory}

  \item{plot}{\link{logical}(with default): draw a plot according to the simulated glow peak or not}
}

\details{
Function \emph{simPeak} simulates a synthetic glow peak according to the Randall-Wilkins model (Randall and Wilkins, 1945) using the 
Bulirsch-Stoer algorithm (Stoer and Bulirsch, 2002). This function may be used for synthesizing a glow curve consisting of several 
glow peaks.\cr\cr
The Randall-Wilkins model can be described by (Bos et al., 1993):\cr
\eqn{\frac{dn}{dt}=-n S e^{-E/(kT)}}\cr
where \eqn{n} is the concentration of trapped electrons, \eqn{\frac{dn}{dt}} the rate of change of the concentration of trapped electrons,
\eqn{S} the frequency factor, \eqn{E} the activation energy, \eqn{T} the absolute temperature, and \eqn{k} the Boltzmann constant. The differential 
equation can be written in terms of the temperature \eqn{T} instead of the time \eqn{t} (Pagonis et al., 2006) by writing:\cr
\eqn{\frac{dn}{dT}=\frac{dn}{dt}/\frac{dT}{dt}=-\frac{n S e^{-E/(kT)}}{beta}}\cr
where \eqn{beta} denotes the linear heating rate. 
}

\value{
Return a list containing the following elements:
    \item{temps}{a vector of temperature values}
    \item{tl}{values of the thermoluminescence intensity}
}

\references{
Randall JT, Wilkins MHF, 1945. Phosphorescence and electron traps. I. The study of trap distributions. Proceedings of the Royal Society of London. Series A, Mathematical and Physical Sciences, 184(999): 365-389.

Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1993. An intercomparison of glow curve analysis computer programs: I. Synthetic glow curves. Radiation Protection Dosimetry, 47(1-4), 473-477.

Stoer J, Bulirsch R, 2002. Introduction to numerical analysis. Texts in Applied Mathematics. Springer Science & Business.

Pagonis V, Kitis G, Furetta C, 2006. Numerical and practical exercises in thermoluminescence. Springer Science & Business Media. 
}

\seealso{
\link{tgcd}
}

\examples{
    # Synthesizing a glow curve consisting of three glow peaks. 
     temps <- seq(308.15, 558.89, by=1.26)         # Temperature.
             
     res1 <- simPeak(ae=1.94, ff=3.02e21, hr=6.0,  # Peak 1.
                     temps=temps, tl0=531.01)      
     res2 <- simPeak(ae=2.15, ff=1.42e22, hr=6.0,  # Peak 2.
                     temps=temps, tl0=5068.95)
     res3 <- simPeak(ae=2.03, ff=8.62e19, hr=6.0,  # Peak 3.
                     temps=temps, tl0=11880.16)
     noise <- rnorm(n=199, mean=0, sd=0.05)        # Random noise.
    
     res <- cbind(res1$tl, res2$tl, res3$tl, 
                  res1$tl+res2$tl+res3$tl+noise)
     matplot(x=res1$temps, y=res, type="l", lwd=3, lty="solid",
             xlab="Temperature (K)", ylab="TL intensity (counts)")
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow curve simulation}

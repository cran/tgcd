\name{tgcd}
\alias{tgcd}
\alias{tgcd.default}
\title{
Thermoluminescence glow curve deconvolution
}
\description{
Thermoluminescence glow curve deconvolution (TGCD).
}
\usage{
tgcd(Sigdata, npeak, inis = NULL, mdt = 3,
     nstart = 200, model=c("f","g"), 
     elim = NULL, logy = FALSE, hr = NULL, 
     outfile = NULL, plot = TRUE)
}

\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a \code{two-column} matrix, temperature values and thermoluminescence signal values are stored in the first and second column, respectively}

  \item{npeak}{\link{integer}(\bold{required}): number of glow peaks (\code{npeak<=13})}

  \item{inis}{\link{matrix}(optional): a \code{npeak-row} \code{3-column} (or \code{4-column}) matrix used for storing starting parameters \code{Im}, \code{E}, \code{Tm} (or \code{b}). 
If \code{inis=NULL}, the user has to click with a mouse on a plot to locate each glow peak maximum}

  \item{mdt}{\link{numeric}(with default): allowed minimum distance between each optimized temperature at maximum thermoluminescence intensity. A large \code{mdt} prevents 
the appearance of strongly overlapping peaks. A conservative \code{mdt} is adopted by default}

  \item{nstart}{\link{integer}(with default): number of trials, if \code{nstart>1}, a "try-and-error" protocol will be performed, the upper limit for
\code{nstart} is set equal to 50000}

  \item{model}{\link{character}(with default): \code{"f"} means fitting a first-order model, \code{"g"} means fitting a general-order model, default \code{model="f"}}

  \item{elim}{\link{vector}(with default): lower and upper limits for the activation energy, default \code{elim=c(0.5, 5)} }

  \item{logy}{\link{logical}(with default): draw the y-axis of the plot used for locating peak maxima with a logarithmic scale or not}

  \item{hr}{\link{numeric}(optional): the linear heating rate used for calculating the frequency factor}

  \item{outfile}{\link{character}(optional): if specified, fitted signal values for each glow peak will be written to a file named \code{"outfile"} in \code{CSV} format and 
saved to the current work directory}

  \item{plot}{\link{logical}(with default): draw a plot according to the fitting result or not}
}

\details{
Function \emph{tgcd} is used for deconvolving first-order (Bos et al., 1993; Puchalska and Bilski, 2006) and general 
order (Kitis et al., 1998, Pagonis et al., 2006) thermoluminescence glow curves using the Levenberg-Marquardt algorithm 
(plus supports for constraining and fixing parameters).

The first-order model for a glow peak is:

\code{I(T)=Im*exp(xa-xb)*exp(xa*(alpha(xa)-T/Tm*exp(xa-xb)*alpha(xb)))}

where \code{xa=E/k/Tm}, \code{xb=E/k/T}, \code{alpha(x)} is a quotient of fourth order polynomials. 
\code{I} is the glow peak intensity, \code{E} the activation energy in eV, \code{k} the Boltzmann constant in 
eV/K (\code{k=8.617385e-5}), T the temperature in K with constant heating rate in K/s, \code{Tm} the temperature 
at maximum thermoluminescence intensity, and \code{Im} the maximum intensity.

The general-order model for a glow peak is:

\code{I(T)=Im*b^(b/(b-1))*expv*((b-1)*(1-xa)*(T/Tm)^2*expv+Zm)^(-b/(b-1))}

where \code{xa=2*k*T/E}, \code{xb=2*k*Tm/E}, \code{expv=exp(E/(k*T)*(T-Tm)/Tm)}, \code{Zm=1+(b-1)*xb}, kinetic 
parameter \code{b} lies between 1 and 2. \code{I} is the glow peak intensity, \code{E} the activation energy, 
\code{k} the Boltzmann constant, T the temperature with constant heating rate, \code{Tm} the temperature at
maximum thermoluminescence intensity, and \code{Im} the maximum intensity.

The procedure minimizes the objective:

\eqn{fcn=\sum_{i=1}^n |y_i^o-y_i^f|, i=1,...,n}

where \eqn{y_i^o} and \eqn{y_i^f} denote the \emph{i-th} observed and fitted signal value, respectively, 
and \eqn{n} indicates the number of data points.

Starting parameters (\code{inis}) can be specified by the user through argument \code{inis} or by clicking with a mouse on 
the plot of the thermoluminescence glow curve showing peak maxima if \code{inis=NULL}.The Levenberg-Marquardt algorithm 
(More, 1978) (\emph{minpack}: Original Fortran 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. Fortran 90 version 
by John Burkardt freely available at \url{http://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html}) was modified so 
as to supports constraints and fixes of parameters. If argument \code{nstart>1}, a "try-and-error" protocol with starting 
values generated uniformly around the given starting values \code{inis} will be performed repeatedly to search the optimal 
parameters that give a minimum figure of merit (FOM) value.

Parameters can be interactively constrained and fixed by modifying the following elements in a automatically 
generated \bold{Dialog Table} if \code{inis=NULL}:\cr 
\bold{(1)} \emph{INTENS(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{Im} \cr
\bold{(2)} \emph{ENERGY(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{E} \cr
\bold{(3)} \emph{TEMPER(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{Tm} \cr
\bold{(4)} \emph{bValue(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{b} \cr
}

\value{
Return a list containing the following elements:
   \item{pars}{optimized parameters stored in a matrix}
   \item{ff}{resultant frequency factor if \code{hr!=NULL}}
   \item{sp}{parameters used for describing the shape of a glow peak, see function \link{simPeak} for details}
   \item{FOM}{minimized figure of merit}
}

\references{
Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1993. An intercomparison of glow curve analysis computer programs: I. Synthetic glow curves. Radiation Protection Dosimetry, 47(1-4), 473-477.

Kitis G, Gomes-Ros JM, Tuyn JWN, 1998. Thermoluminescence glow curve deconvolution functions for first, second and general orders of kinetics. Journal of Physics D: Applied Physics, 31(19): 2636-2641.

More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory," in Lecture Notes in Mathematics: Numerical Analysis, Springer-Verlag: Berlin. 105-116.

Pagonis V, Kitis G, Furetta C, 2006. Numerical and practical exercises in thermoluminescence. Springer Science & Business Media. 

Puchalska M, Bilski P, 2006. GlowFit-a new tool for thermoluminescence glow-curve deconvolution. Radiation Measurements, 41(6): 659-664. 
Software is freely available at \url{http://www.ifj.edu.pl/dept/no5/nz58/deconvolution.htm}.

\bold{Further reading} 

Chung KS, Choe HS, Lee JI, Kim JL, Chang SY, 2005. A computer program for the deconvolution of thermoluminescence glow curves. Radiation Protection Dosimetry, 115(1-4): 345-349.
Software is freely available at \url{http://physica.gsnu.ac.kr/TLanal}.

Harvey JA, Rodrigues ML, Kearfott JK, 2011. A computerized glow curve analysis (GCA) method for WinREMS thermoluminescent dosimeter data using MATLAB. Applied Radiation and Isotopes, 69(9):1282-1286.
Source codes are freely available at \url{http://www.sciencedirect.com/science/article/pii/S0969804311002685}.

Kiisk V, 2013. Deconvolution and simulation of thermoluminescence glow curves with Mathcad. Radiation Protection Dosimetry, 156(3): 261-267.
Software is freely available at \url{http://www.physic.ut.ee/~kiisk/mcadapps.htm}.
}

\note{
The model to be optimized should not be underdetermined. This means that the number of data points should exceed the number of parameters. 
A lack of background counts in the analyzed data is assumed. To obtain reliable estimate, the presented background may be accounted for 
by subtracting from measured data before analysis.
}

\seealso{
\link{simPeak}; \link{simqOTOR}
}

\examples{
# Load the data.
  data(Refglow)
# Deconvolve Refglow002 with 4 peaks.
  startingPars <- cbind(c(400, 550, 850, 1600), # Im
                        c(1.4, 1.5, 1.6, 2),    # E
                        c(420, 460, 480, 510))  # Tm
  tgcd(Refglow$Refglow002, npeak=4, model="f",
       inis=startingPars, nstart=10)
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow curve deconvolution}

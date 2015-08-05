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
     nstart = 30, elim = NULL, logy = FALSE, 
     hr = NULL, outfile = NULL, plot = TRUE)
}

\arguments{
  \item{Sigdata}{\link{matrix}(\bold{required}): a \code{two-column} matrix, temperature values and thermoluminescence signal values are stored in the first and second column, respectively}

  \item{npeak}{\link{integer}(\bold{required}): number of glow peaks (\code{npeak<=13})}

  \item{inis}{\link{matrix}(optional): a \code{npeak-row} \code{3-column} matrix used for storing starting parameters (\code{Im}, \code{E}, \code{Tm} are stored 
in the first, second, and third column, respectively, see \bold{Examples}). If \code{inis=NULL}, the user has to click with a mouse on a plot to locate each glow peak maximum}

  \item{mdt}{\link{numeric}(with default): allowed minimum distance between each optimized temperature at maximum thermoluminescence intensity. A large \code{mdt} prevents 
the appearance of strongly overlapping peaks. A conservative \code{mdt} is adopted by default}

  \item{nstart}{\link{integer}(with default): number of trials, if \code{nstart>1}, a "try-and-error" protocol will be performed, the upper limit for
\code{nstart} is set equal to 1000}

  \item{elim}{\link{vector}(with default): lower and upper limits for the activation energy, default \code{elim=c(0.5, 5)} }

  \item{logy}{\link{logical}(with default): draw the y-axis of the plot used for locating peak maxima with a logarithmic scale or not}

  \item{hr}{\link{numeric}(optional): the linear heating rate used for calculating the frequency factor}

  \item{outfile}{\link{character}(optional): if specified, fitted signal values for each glow peak will be written to a file named \code{"outfile"} in \code{CSV} format and 
saved to the current work directory}

  \item{plot}{\link{logical}(with default): draw a plot according to the fitting result or not}
}

\details{
Function \emph{tgcd} is used for deconvolving first-order thermoluminescence glow curves approximated by a quotient of 
polynomials (Bos et al., 1993; Puchalska and Bilski, 2006) using the Levenberg-Marquardt algorithm (plus supports for constraining and fixing parameters).\cr\cr
The fitting model for a single glow peak is:\cr
\code{I(T)=Im*exp(xa-xb)*exp(xa*(alpha(xa)-T/Tm*exp(xa-xb)*alpha(xb)))}\cr
where \code{xa=E/k/Tm}, \code{xb=E/k/T}, \code{alpha(x)} is a quotient of fourth order polynomials, 
\code{I} is the glow peak intensity, \code{E} the activation energy in eV, \code{k} the Boltzmann constant in eV/K (\code{k=8.617385e-5}), 
T the temperature in K with constant heating rate in K/s, \code{Tm} the temperature at maximum thermoluminescence intensity, and \code{Im} 
the maximum intensity.\cr\cr
The procedure minimizes the objective:\cr
\eqn{fcn=\sum_{i=1}^n |y_i^o-y_i^f|, i=1,...,n}\cr
where \eqn{y_i^o} and \eqn{y_i^f} denote the \emph{i-th} observed and fitted signal value, respectively, and \eqn{n} indicates 
the number of data points.\cr\cr
Starting parameters (\code{inis}) can be specified by the user through argument \code{inis} or by clicking with a mouse on 
the plot of the thermoluminescence glow curve showing peak maxima if \code{inis=NULL}. A modification of the Levenberg-Marquardt algorithm 
(minpack: Original Fortran 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. Fortran 90 version by John Burkardt, which is freely 
available at \url{http://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.html}) that supports constraints and fixes of parameters is used for 
parameter optimization.\cr\cr
If argument \code{nstart>1}, a "try-and-error" protocol with starting values generated uniformly from the 
space \code{(0.9*inis,1.1*inis)} around the specified starting values \code{inis} will be performed 
repeatedly to search the optimal parameters that results in a minimum figure of merit (FOM) value.\cr\cr
Parameters can be interactively constrained and fixed by modifying the following elements in a automatically 
generated \bold{Dialog Table} if \code{inis=NULL}:\cr 
\bold{(1)} \emph{INTENS(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{Im} \cr
\bold{(2)} \emph{ENERGY(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{E} \cr
\bold{(3)} \emph{TEMPER(min, max, ini, fix)}: lower and upper bounds, starting and fixing values of \code{Tm} \cr
}

\value{
Return a list containing the following elements:
   \item{pars}{optimized parameters stored in a matrix}
   \item{stdpars}{standard errors of parameters estimated by numerical difference-approximation}
   \item{ff}{resultant frequency factor if \code{hr!=NULL}}
   \item{FOM}{minimized figure of merit}
}

\references{
Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1993. An intercomparison of glow curve analysis computer programs: I. Synthetic glow curves. Radiation Protection Dosimetry, 47(1-4), 473-477.

More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory," in Lecture Notes in Mathematics: Numerical Analysis, Springer-Verlag: Berlin. 105-116.

Pagonis V, Kitis G, Furetta C, 2006. Numerical and practical exercises in thermoluminescence. Springer Science & Business Media. 

Puchalska M, Bilski P, 2006. GlowFit-a new tool for thermoluminescence glow-curve deconvolution. Radiation Measurements, 41(6): 659-664. 
Software is freely available at \url{http://www.ifj.edu.pl/dept/no5/nz58/deconvolution.htm}.

\bold{Further reading} 

Abd El-Hafez AI, Yasin MN, Sadek AM, 2011. GCAFIT-A new tool for glow curve analysis in thermoluminescence nanodosimetry. Nuclear Instruments and Methods in Physics Research A, 637(1): 158-163.

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
\link{simPeak}
}

\examples{
# Load the data.
  data(Refglow)

# Deconvolve Refglow002 with 4 peaks.
  startingPars <- cbind(c(400, 550, 850, 1600), # Im
                        c(1.4, 1.5, 1.6, 2),    # E
                        c(420, 460, 480, 510))  # Tm
  tgcd(Refglow$Refglow002, npeak=4, 
       inis=startingPars, nstart=10)

# Deconvolve Refglow009 with 9 peaks, not run.
# startingPars <- cbind(
#      c(9680, 19500, 23420, 60010, 5580, 4250, 6720, 1790, 2260), # Im
#      c(1.2, 1.3, 1.6, 2.0,  2.0, 1.4, 1.6, 1.9, 1.0),            # E
#      c(390, 430, 460, 490, 500, 520, 550, 580, 600)              # Tm 
#  )
#  tgcd(Refglow$Refglow009, npeak=9, 
#       inis=startingPars, nstart=10, elim=c(1,3))
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow curve deconvolution}

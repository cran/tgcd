\name{tgcd}
\alias{tgcd}
\alias{tgcd.default}
\title{
Thermoluminescence glow curve deconvolution (tgcd)
}
\description{
Thermoluminescence glow curve deconvolution according to various first-order, second-order, general-order, and mixed-order empirical expressions.
}
\usage{
tgcd(Sigdata, npeak, model = "g1", subBG = FALSE, pickp = "d2", 
     pickb = "d0", nstart = 60, kkf = 0.03, mdt = NULL, mwt = NULL,  
     mr = NULL, edit.inis = TRUE, inisPAR = NULL, inisBG = NULL, 
     hr = NULL, hwd = NULL, pod = NULL, plot = TRUE, outfile = NULL)
}

\arguments{
\item{Sigdata}{\link{matrix}(\bold{required}): a \code{2}-column matrix, temperature values (in unit K) and thermoluminescence signal values are 
stored in the first and second column, respectively}

\item{npeak}{\link{integer}(\bold{required}): number of glow peaks, the allowed maximum number of glow peaks is set equal to 13}

\item{model}{\link{character}(with default): model used for glow curve deconvolution, \code{"f1"}, \code{"f2"}, and \code{"f3"} for first-order models, \code{"s1"}, 
\code{"s2"} for second-order models, \code{"g1"}, \code{"g2"}, \code{"g3"} for general-order models, \code{"wo"} and \code{"lw"} for the Wright Omega and the Lambert W functions,
\code{"m1"}, \code{"m2"}, and \code{"m3"} for the mixed-order models (see \bold{Details})}

\item{subBG}{\link{logical}(with default): whether the user want to subtract the background during the deconvolution}

\item{pickp}{\link{character}(with default): mode used for initialization of kinetic parameters if \code{inisPAR=NULL}, 
\code{"d0"} and \code{"d01"} prompt the user to click with a mouse on the original and log-scale glow curves respectively to locate each glow peak,
\code{"d1"}, \code{"d2"}, \code{"d3"}, and \code{"d4"} prompt the user to click with a mouse on the first-, second-, third-, 
and fourth-derivative of the glow curve respectively to locate each glow peak}

\item{pickb}{\link{character}(with default): mode used for initialization of background parameters if \code{inisBG=NULL}, 
\code{"d0"} and \code{"d01"} prompt the user to click with a mouse on the original and log-scale glow curves respectively to initilize the background parameters}

\item{nstart}{\link{integer}(with default): number of trials performed in a "trial-and-error" protocol, the upper limit is set equal to 10000}

\item{kkf}{\link{numeric}(with default): factor controlling the range of values from which random starting parameters will be generated during 
the "trial-and-error" protocol, \code{0<kkf<1}. For example, if \code{kkf=0.03} then kinetic parameters will be generated uniformly between 
\code{(1.0-kkf)*inisPAR} and \code{(1.0+kkf)*inisPAR} and background parameters will be generated uniformly between 
\code{(1.0-kkf)*inisBG} and \code{(1.0+kkf)*inisBG}}

\item{mdt}{\link{numeric}(optional): allowed minimum distance between \code{Tm} values of glow peaks. 
A larger \code{mdt} prevents the appearance of strongly overlapping peaks}

\item{mwt}{\link{numeric}(optional): allowed maximum total half-width of deconvoluted glow peaks. 
A smaller \code{mwt} prevents the appearance of glow peaks with large total half-width}

\item{mr}{\link{numeric}(optional): allowed minimum resolution of glow peaks. 
A larger \code{mr} prevents the appearance of strongly overlapping peaks}

\item{edit.inis}{\link{logical}(with default): whether the user want to further modify, constrain, or fix 
the initialized kinetic (and/or background) parameters through an automatically generated \bold{Dialog Table}}

\item{inisPAR}{\link{matrix}(optional): a matrix (\code{3} or \code{4} columns) used for storing initial
kenetic parameters \code{Im}, \code{E}, \code{Tm}, \code{b} (or \code{R}, \code{a}) (see \bold{Examples})}

\item{inisBG}{\link{vector}(optional): a 3-element vector containing initial background parameters 
\code{A}, \code{B}, and \code{C} used for background subtraction (see \bold{Examples})}

\item{hr}{\link{numeric}(optional): linear heating rate used for calculating the frequency factor}

\item{hwd}{\link{integer}(with default): half width (length) of the segement used for Savitzky-Golay smoothing}

\item{pod}{\link{integer}(with default): order of the polynomial used for Savitzky-Golay smoothing}

\item{plot}{\link{logical}(with default): draw a plot according to the fitting result or not}

\item{outfile}{\link{character}(optional): if specified, fitted signal values for each glow peak will be written to 
a file named \code{"outfile"} in \code{CSV} format and saved to the current work directory}
}

\details{
Function \link{tgcd} is used for deconvolving thermoluminescence glow curves according to various kinetic models. 
In the text below, \eqn{I(T)} is the thermoluminescence intensity as function of temperature \eqn{T}, \eqn{E} the activation energyin ev, \eqn{k} the Boltzmann constant in eV/k, 
\eqn{T} the temperature in K with constant heating rate K/s, \eqn{T_m} the temperature at maximum thermoluminescence intensity in K, 
\eqn{I_m} the maximum intensity, \eqn{b} is an extra parameter (the kinetic order) in application of a general-order model, 
\eqn{R} is an extra parameter in application of the Lambert W (Wright Omega) function,
and \eqn{\alpha} is an extra parameter in application of the mixed-order models.\cr

\bold{First-order} glow peaks appear if the recombination probability (\code{A_m}) is greater than that of re-trapping (\code{A_n}) during excitation. 
The three parameters describing a glow peak are: \eqn{I_m}, \eqn{E}, and \eqn{T_m}.
Three empirical expressions describing first-order glow peaks are available in function \link{tgcd}:\cr
\strong{<1>}The first type of first-order empirical expression (\code{model="f1"}) is (Bos et al., 1993a)\cr\cr
\eqn{I(T)=I_mexp(\frac{E}{kT_m}-\frac{E}{kT})exp[\frac{E}{kT_m}\alpha(\frac{E}{kT_m})-(\frac{T}{T_m})\alpha(\frac{E}{kT})exp(\frac{E}{kT_m}-\frac{E}{kT})]}\cr\cr
where \eqn{\alpha(x)} is a quotient of 4th-order polynomials of the form\cr\cr
\eqn{\alpha(x)=\frac{a_0+a_1x+a_2x^2+a_3x^3+x^4}{b_0+b_1x+b_2x^2+b_3x^3+x^4}}\cr\cr
\eqn{a_0=0.267773734}, \eqn{a_1=8.6347608925}, \eqn{a_2=8.059016973}, \eqn{a_3=8.5733287401},\cr
\eqn{b_0=3.9584969228}, \eqn{b_1=21.0996530827}, \eqn{b_2=25.6329561486}, \eqn{b_3=9.5733223454}\cr

\strong{<2>}The second type of first-order empirical expression (\code{model="f2"}) is (Kitis et al., 1998)\cr\cr
\eqn{I(T)=I_mexp[1+\frac{E}{kT}\frac{T-T_m}{T_m}-\frac{T^2}{T_m^2}exp(\frac{E}{kT}\frac{T-T_m}{T_m})(1-\frac{2kT}{E})-\frac{2kT_m}{E}]}\cr

\strong{<3>}The third type of first-order function fits a weibull function (\code{model="f3"}) (Pagonis et al., 2001)\cr\cr
\eqn{I(T)=2.713I_m(\frac{T-T_m}{b}+0.996)^{15}exp[-(\frac{T-T_m}{b}+0.996)^{16}]}\cr\cr
where \eqn{b=\sqrt{\frac{242.036T_m^3k^2}{(E+T_mk)^2-7T_m^2k^2}}}\cr

\bold{Second-order} glow peaks appear if the re-trapping probability is comparable with or greater than that of recombination during excitation. 
The three parameters describing a glow peak are the same as those of first-orders.
Two empirical expressions describing second-order glow peaks are available in function \link{tgcd}:\cr
\strong{<4>}The first type of second-order  empirical expression (\code{model="s1"}) is (Kitis et al., 1998)\cr\cr
\eqn{I(T)=4I_mexp(\frac{E}{kT}\frac{T-T_m}{T_m})[\frac{T^2}{T_m^2}(1-\frac{2kT}{E})exp(\frac{E}{kT}\frac{T-T_m}{T_m})+1+\frac{2kT_m}{E}]^{-2}}\cr

\strong{<5>}The second type of second-order function fit a logistic asymmetric function (\code{model="s2"}) (Pagonis and Kitis, 2001)\cr\cr
\eqn{I(T)=5.2973I_m[1+exp(-(\frac{T-T_m}{a_2}+0.38542))]^{-2.4702}exp[-(\frac{T-T_m}{a_2}+0.38542)]}\cr\cr
where \eqn{a_2=\sqrt{\frac{1.189T_m^4k^2}{E^2+4ET_mk}}}\cr

\bold{General-order} glow peaks are produced in intermediate cases (neither of first-order, nor of second-order). 
The four parameters describing a glow peak are: \eqn{I_m}, \eqn{E}, \eqn{T_m}, and \eqn{b}.
Three empirical expressions describing general-order glow peaks are available in function \link{tgcd}:\cr
\strong{<6>}The first type of general-order empirical expression (\code{model="g1"}) is (Kitis et al., 1998)\cr\cr
\eqn{I(T)=I_mb^{\frac{b}{b-1}}exp(\frac{E}{kT}\frac{T-T_m}{T_m})[(b-1)(1-\frac{2kT}{E})\frac{T^2}{T_m^2}exp(\frac{E}{kT}\frac{T-T_m}{T_m})+1+\frac{2kT_m(b-1)}{E}]^{-\frac{b}{b-1}}}\cr

\strong{<7>}The second type of general-order empirical expression (\code{model="g2"}) is (Gomez-Ros and Kitis, 2002)\cr\cr
\eqn{I(T)=I_mexp(\frac{E}{kT_m}-\frac{E}{kT})[1+\frac{b-1}{b}\frac{E}{kT_m}(\frac{T}{T_m}exp(\frac{E}{kT_m}-\frac{E}{kT})F(\frac{E}{kT})-F(\frac{E}{kT_m}))]^{-\frac{b}{b-1}}}\cr\cr
where \eqn{F(x)} is a rational approximation function of the form\cr\cr
\eqn{F(x)=1-\frac{a_0+a_1x+x^2}{b_0+b_1x+x^2}}\cr\cr
\eqn{a_0=0.250621}, \eqn{a_1=2.334733}, \eqn{b_0=1.681534}, \eqn{b_1=3.330657}\cr

\strong{<8>}The third type of general-order empirical expression (\code{model="g3"}) is (Gomez-Ros and Kitis, 2002)\cr\cr
\eqn{I(T)=I_mexp(\frac{E}{kT_m^2}(T-T_m))[\frac{1}{b}+\frac{b-1}{b}exp(\frac{E}{kT_m^2}(T-T_m))]^{-\frac{b}{b-1}}}\cr

\bold{One trap-one recombination (OTOR) model} based semi-analytical expressions have also been applied to fit glow peaks,
by using either the Lambert W function (Kitis and Vlachos, 2013; Sadek et al., 2015; Kitis et al., 2016) or 
the Wright Omega function (Singh and Gartia, 2013; 2014; 2015).
The four parameters describing a glow peak are: \eqn{I_m}, \eqn{E}, \eqn{T_m}, and \eqn{R=\frac{A_n}{A_m}} 
(where \eqn{A_n} and \eqn{A_m} represent the retrapping and recombination probabilities, respectively).
Two analytical expressions describing the OTOR model are available in function \link{tgcd}:\cr
\strong{<9>}The semi-analytical expression derived using the Wright Omega function (\eqn{R=\frac{A_n}{A_m}<1}) can be described as (\code{model="wo"})\cr\cr
\eqn{I(T)=I_mexp(-\frac{E}{kT}\frac{T_m-T}{T_m})\frac{\omega(Z_m)+[\omega(Z_m)]^2}{\omega(Z)+[\omega(Z)]^2}}\cr\cr
where \eqn{Z_m=\frac{R}{1-R}-log(\frac{1-R}{R})+\frac{Eexp(\frac{E}{kT_m})}{kT_m^2(1-1.05R^{1.26})}F(T_m,E)},\cr\cr
and \eqn{Z=\frac{R}{1-R}-log(\frac{1-R}{R})+\frac{Eexp(\frac{E}{kT_m})}{kT_m^2(1-1.05R^{1.26})}F(T,E)},\cr\cr
\strong{<10.1>}The semi-analytical expression derived using the Lambert W function for \eqn{R=\frac{A_n}{A_m}<1} can be described as (\code{model="lw"})\cr\cr
\eqn{I(T)=I_mexp(-\frac{E}{kT}\frac{T_m-T}{T_m})\frac{W(exp(Z_m))+[W(exp(Z_m))]^2}{W(exp(Z))+[W(exp(Z))]^2}}\cr\cr
where \eqn{Z_m=\frac{R}{1-R}-log(\frac{1-R}{R})+\frac{Eexp(\frac{E}{kT_m})}{kT_m^2(1-1.05R^{1.26})}F(T_m,E)},\cr\cr
and \eqn{Z=\frac{R}{1-R}-log(\frac{1-R}{R})+\frac{Eexp(\frac{E}{kT_m})}{kT_m^2(1-1.05R^{1.26})}F(T,E)},\cr\cr
\strong{<10.2>}The semi-analytical expression derived using the Lambert W function for \eqn{R=\frac{A_n}{A_m}>1} can be described as (\code{model="lw"})\cr\cr
\eqn{I(T)=I_mexp(-\frac{E}{kT}\frac{T_m-T}{T_m})\frac{W(-1,-exp(-Z_m))+[W(-1,-exp(-Z_m))]^2}{W(-1,-exp(-Z))+[W(-1,-exp(-Z))]^2}}\cr\cr
where \eqn{Z_m=|\frac{R}{1-R}|+log(|\frac{1-R}{R}|)+\frac{Eexp(\frac{E}{kT_m})}{kT_m^2(2.963-3.24R^{-0.74})}F(T_m,E)},\cr\cr
and \eqn{Z=|\frac{R}{1-R}|+log(|\frac{1-R}{R}|)+\frac{Eexp(\frac{E}{kT_m})}{kT_m^2(2.963-3.24R^{-0.74})}F(T,E)}\cr\cr
\eqn{F(T_m,E)} and \eqn{F(T,E)} are described as follows\cr\cr
\eqn{F(T_m,E)=T_mexp(-\frac{E}{kT_m})+\frac{E}{k}Ei(-\frac{E}{kT_m})},\cr\cr
\eqn{F(T,E)=Texp(-\frac{E}{kT})+\frac{E}{k}Ei(-\frac{E}{kT})}\cr\cr
where \eqn{\omega(x)} and \eqn{Ei(x)} are the wright Omega function and the exponential integral function for variable \eqn{x}, respectively.
\eqn{W(x)} and \eqn{W(-1,x)} are the principal and the second branches of the Lambert W function, respectively.
The Fortran 90 subroutine used for evaluating the Wright Omega function is transformed from the Matlab code provided by Andrew Horchler (\url{https://github.com/horchler/wrightOmegaq}).
The Fortran 90 subroutine (original Fortran 77 version by William Cody) used for evaluating the Lambert W function written by John Burkardt is available at 
\url{https://people.sc.fsu.edu/~jburkardt/f_src/toms743/toms743.f90}. 
The Fortran 90 subroutine used for evaluating the exponential integral function is written by John Burkardt (original Fortran 77 version by William Cody) 
(\url{https://people.sc.fsu.edu/~jburkardt/f_src/toms715/toms715.f90}).\cr

\bold{Mixed-order} kinetic models introduce an extra parameter \eqn{\alpha=\frac{n_0}{n_0+M}} where \eqn{n_0} is the initial filled concentration of the active traps and \eqn{M}
is the trap concentration of the thermally disconnected deep traps (Sunta et al., 2002). The four parameters describing a glow peak are: \eqn{I_m}, \eqn{E}, \eqn{T_m}, and \eqn{\alpha}.
Three empirical expressions describing mixed-order glow peaks are available in function \link{tgcd}:\cr
\strong{<11>}The first type of mixed-order empirical expression (\code{model="m1"}) is (Kitis and Gomez-Ros, 2000)\cr\cr
\eqn{I(T)=\frac{I_m[exp(\frac{1-\frac{2kT_m}{E}}{R_m})-\alpha]^2exp(\frac{E}{kT}\frac{T-T_m}{T_m})exp[\frac{T^2}{T_m^2R_m}exp(\frac{E}{kT}\frac{T-T_m}{T_m})
(1-\frac{2kT}{E})]}{exp(\frac{1-\frac{2kT_m}{E}}{R_m})[exp[\frac{T^2}{T_m^2R_m}exp(\frac{E}{kT}\frac{T-T_m}{T_m})(1-\frac{2kT}{E})]-\alpha]^2}}\cr\cr\cr
where \eqn{R_m=\frac{A_m+\alpha}{A_m-\alpha}} and \eqn{A_m=exp(\frac{A_m-\alpha}{A_m+\alpha}(1-\frac{2kT_m}{E}))}\cr

\strong{<12>}The second type of mixed-order empirical expression (\code{model="m2"}) is (Gomez-Ros and Kitis, 2002)\cr\cr
\eqn{I(T)=\frac{4I_mR_m^2exp(\frac{E}{kT_m}-\frac{E}{kT})exp[R_m\frac{E}{kT_m}(\frac{T}{T_m}exp(\frac{E}{kT_m}-\frac{E}{kT})F(\frac{E}{kT})-F(\frac{E}{kT_m}))]}
{(1+R_m)[exp[R_m\frac{E}{kT_m}(\frac{T}{T_m}exp(\frac{E}{kT_m}-\frac{E}{kT})F(\frac{E}{kT})-F(\frac{E}{kT_m}))]-(1-R_m)]^2}}\cr\cr
where \eqn{R_m=(1-\alpha)(1+0.2922\alpha-0.2783\alpha^2)}\cr\cr
\eqn{F(x)=1-\frac{a_0+a_1x+x^2}{b_0+b_1x+x^2}}\cr\cr
\eqn{a_0=0.250621}, \eqn{a_1=2.334733}, \eqn{b_0=1.681534}, \eqn{b_1=3.330657}\cr

\strong{<13>}The third type of mixed-order empirical expression (\code{model="m3"}) is (Vejnovic et al., 2008)\cr\cr
\eqn{I(T)=\frac{I_m\alpha(2-l)^2exp[\frac{T^2}{T_m^2}(\frac{2}{l}-1)exp(\frac{E}{kT_m}\frac{T-T_m}{T_m})(1-\frac{2kT}{E})]exp(\frac{E}{kT_m}\frac{T-T_m}{T_m})}{
(l-1)[exp[\frac{T^2}{T_m^2}(\frac{2}{l}-1)exp(\frac{E}{kT_m}\frac{T-T_m}{T_m})(1-\frac{2kT}{E})]-\alpha]^2}}\cr\cr\cr
where \eqn{\alpha=(l-1)exp[\frac{2-l}{l}(1-\frac{2kT_m}{E})]}\cr

The \strong{background} will be subtracted using the following expression if \code{subBG=TRUE} (Horowitz and Yossian, 1995; Kitis et al., 2012)\cr\cr
\eqn{I(T)=A+Bexp(\frac{T}{C})}\cr\cr
where \eqn{A}, \eqn{B}, and \eqn{C} are positive parameters to be optimized.\cr

The procedure minimizes the objective:
\eqn{fcn=\sum_{i=1}^n |y_i^o-y_i^f|, i=1,...,n}
where \eqn{y_i^o} and \eqn{y_i^f} denote the \emph{i-th} observed and fitted signal value, respectively, 
and \eqn{n} indicates the number of data points.\cr

The Levenberg-Marquardt algorithm (More, 1978) (\emph{minpack}: \url{https://netlib.org/minpack/}, original Fortran 77 version by Jorge More, Burton Garbow, Kenneth Hillstrom. 
Fortran 90 version by John Burkardt \url{https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.f90}) was modified so as to supports constraints and fixes of parameters.
A "trial-and-error" protocol with starting values generated uniformly around the given starting parameters \code{inisPAR} and \code{inisBG} is performed 
repeatedly to search the optimal parameters that give a minimum Figure Of Merit (FOM) value (Balian and Eddy, 1977).\cr

Kinetic parameters can be initialized by the user through argument \code{inisPAR} or by clicking with a mouse on 
the plot of the thermoluminescence glow curve showing peak maxima if \code{inisPAR=NULL}.
Background parameters can be initialized by the user through argument \code{inisBG} or by clicking with a mouse on 
the plot of the thermoluminescence glow curve to select 4 data points if \code{inisBG=NULL}.\cr

Parameters can be interactively constrained and fixed by modifying the following elements in a automatically 
generated \bold{Dialog Table} if \code{edit.inis=TRUE}:\cr 
\bold{(1)} \code{INTENS[min,max,ini,fix]}: lower and upper bounds, starting and fixing values of \eqn{I_m} \cr
\bold{(2)} \code{ENERGY[min,max,ini,fix]}: lower and upper bounds, starting and fixing values of \eqn{E} \cr
\bold{(3)} \code{TEMPER[min,max,ini,fix]}: lower and upper bounds, starting and fixing values of \eqn{T_m} \cr
\bold{(4)} \code{bValue[min,max,ini,fix]}: lower and upper bounds, starting and fixing values of \eqn{b} in the general-order model\cr
\bold{(5)} \code{rValue[min,max,ini,fix]}: lower and upper bounds, starting and fixing values of \eqn{R} in the OTOR model\cr
\bold{(6)} \code{aValue[min,max,ini,fix]}: lower and upper bounds, starting and fixing values of \eqn{\alpha} in the mixed-order model\cr
\bold{(7)} \code{BG[min,max,ini,fix]}: lower and upper bounds, starting and fixing values of background parameters \eqn{A}, \eqn{B}, and \eqn{C}\cr
}

\value{
Return an \bold{invisible} list containing the following elements:
   \item{comp.sig}{calculated signal values for each glow peak}
   \item{residuals}{calculated residual values}
   \item{pars}{optimized parameters stored in a matrix}
   \item{BGpars}{optimized background parameters, it returns \code{NULL} if \code{subBG=FALSE}}
   \item{ff}{calculated frequency factors, it returns \code{NULL} if \code{hr=NULL}}
   \item{sp}{parameters used for describing the shape of a glow peak, see function \link{simPeak} for details}
   \item{resolution}{resolutions of optimized glow peaks calculated after Kitis and Pagonis (2019), it returns \code{NULL} if \code{npeak=1}}
   \item{SSR}{Squared Sum of Residuals}
   \item{RCS}{Reduced Chi-Square value}
   \item{R2}{squared "pearson" correlation between observed and fitted signals}
   \item{FOM}{Figure Of Merit value calculated after Balian and Eddy (1977)}
}

\references{
Balian HG, Eddy NW, 1977. Figure-of-merit (FOM), an improved criterion over the normalized chi-squared test for assessing goodness-of-fit of gamma-ray spectral peaks. 
Nuclear and Instruments Methods, 145: 389-395.

Bos AJJ, Piters TM, Gomez-Ros JM, Delgado A, 1993a. An intercomparison of glow curve analysis computer programs. 
IRI-CIEMAT Report, 131-93-005.

Gomez-Ros JM, Kitis G, 2002. Computerised glow curve deconvolution using general and mixed order kinetics. 
Radiation Protection Dosimetry, 101(1-4): 47-52.

Horowitz YS, Yossian D, 1995. Computerized glow curve deconvolution: application to thermoluminescence dosimetry. 
Radiation Protection Dosimetry, 60: 1-114.

Kitis G, Carinou E, Askounis P, 2012. Glow-curve de-convolution analysis of TL glow-curve from constant temperature hot gas TLD readers.
Radiation Measurements, 47: 258-265.

Kitis G, Gomes-Ros JM, 2000. Thermoluminescence glow-curve deconvolution functions for mixed order of kinetics and continuous trap distribution. 
Nuclear Instruments and Methods in Physics Research A, 440: 224-231.

Kitis G, Gomes-Ros JM, Tuyn JWN, 1998. Thermoluminescence glow curve deconvolution functions for first, second and general orders of kinetics. 
Journal of Physics D: Applied Physics, 31(19): 2636-2641.

Kitis G, Pagonis V, 2019. On the resolution of overlapping peaks in complex thermoluminescence glow curves. 
Nuclear Instruments and Methods in Physics Research, A, 913: 78-84.

Kitis G, Polymeris GS, Pagonis V, 2019. Stimulated luminescence emission: From phenomenological models to master analytical equations. 
Applied Radiation and Isotopes, 153: 108797.

Kitis G, Polymeris GS, Sfampa IK, Prokic M, Meric N, Pagonis V, 2016. Prompt isothermal decay of thermoluminescence in MgB4O7:Dy, Na and LiB4O7:Cu, In dosimeters. 
Radiation Measurements, 84: 15-25.

Kitis G, Vlachos ND, 2013. General semi-analytical expressions for TL, OSL and other luminescence stimulation modes derived from the OTOR model using the Lambert W-function. 
Radiation Measurements, 48: 47-54.

More JJ, 1978. "The Levenberg-Marquardt algorithm: implementation and theory," in Lecture Notes in Mathematics: Numerical Analysis, Springer-Verlag: Berlin. 105-116.

Pagonis V, Mian SM, Kitis G, 2001. Fit of first order thermoluminescence glow peaks using the weibull distribution function. 
Radiation Protection Dosimetry, 93(1): 11-17.

Pagonis V, Kitis G, 2001. Fit of second order thermoluminescence glow peaks using the logistic distribution function. 
Radiation Protection Dosimetry, 95(3): 225-229.

Sadek AM, Eissa HM, Basha AM, Carinou E, Askounis P, Kitis G, 2015. The deconvolution of thermoluminescence glow-curves using general expressions 
derived from the one trap-one recombination (OTOR) level model. Applied Radiation and Isotopes, 95: 214-221.

Singh LL, Gartia RK, 2013. Theoretical derivation of a simplified form of the OTOR/GOT differential equation. 
Radiation Measurements, 59: 160-164.

Singh LL, Gartia RK, 2014. Glow-curve deconvolution of thermoluminescence curves in the simplified OTOR equation using the Hybrid Genetic Algorithm. 
Nuclear Instruments and Methods in Physics Research B, 319: 39-43.

Singh LL, Gartia RK, 2015. Derivation of a simplified OSL OTOR equation using Wright Omega function and its application. 
Nuclear Instruments and Methods in Physics Research B, 346: 45-52.

Sunta CM, Ayta WEF, Chubaci JFD, Watanabe S, 2002. General order and mixed order of thermoluminescence glow curves-a comparison. 
Radiation Measurements, 35: 47-57.

Vejnovic Z, Pavlovic MB, Davidovi M, 2008. Thermoluminescence glow curve deconvolution function for the mixed-order kinetics. 
Radiation Measurements, 43: 1325-1330.

\bold{Further reading} 

Afouxenidis D, Polymeris GS, Tsirliganis NC, Kitis G., 2011. Computerised curve deconvolution of TL/OSL curves using a popular spreadsheet program. 
Radiation Protection Dosimetry, 1-8.

Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1993b. An intercomparison of glow curve analysis computer programs: I. Synthetic glow curves. 
Radiation Protection Dosimetry, 47(1-4), 473-477.

Chung KS, Choe HS, Lee JI, Kim JL, Chang SY, 2005. A computer program for the deconvolution of thermoluminescence glow curves. 
Radiation Protection Dosimetry, 115(1-4): 345-349.

Delgado A, Gomez-Ros JM, 2001. Computerised glow curve analysis: a tool for routine thermoluminescence dosimetry. 
Radiation Protection Dosimetry, 96(1-3): 127-132.

El-Hafez AI, Yasin MN, Sadek AM, 2011. GCAFIT-A new tool for glow curve analysis in thermoluminescence nanodosimetry. 
Nuclear Instruments and Methods in Physics Research A, 637: 158-163.

Espinosa G, Castano VM, 1991. Automated derivative analysis of thermoluminescence glow curves. 
Applied Radiation and Isotopes, 42(4): 377-381.

Harvey JA, Rodrigues ML, Kearfott JK, 2011. A computerized glow curve analysis (GCA) method for WinREMS thermoluminescent dosimeter data using MATLAB. 
Applied Radiation and Isotopes, 69(9):1282-1286.

Horowitz YS, Moscovitch M, 2012. Highlights and pitfalls of 20 years of application of computerised glow curve ananlysis to thermoluminescence research and dosimetry. 
Radiation Protection Dosimetry, 1-22.

Horowitz Y, Delgado A, Pradhan AS, Yoder RC, 2002. Topics under debate: the use of computerised glow curve analysis will optimise personal thermoluminescence dosimetry. 
Radiation Protection Dosimetry, 102(3): 269-277.

Karmakar M, Bhattacharyya S, Sarkar A, Mazumdar PS, Singh SD, 2017. Analysis of thermoluminescence glow curves using derivatives of different orders. 
Radiation Protection Dosimetry, 1-10.

Kiisk V, 2013. Deconvolution and simulation of thermoluminescence glow curves with Mathcad. 
Radiation Protection Dosimetry, 156(3): 261-267.

Pagonis V, Kitis G, Furetta C, 2006. Numerical and practical exercises in thermoluminescence. 
Springer Science & Business Media. 

Puchalska M, Bilski P, 2006. GlowFit-a new tool for thermoluminescence glow-curve deconvolution. 
Radiation Measurements, 41(6): 659-664. 

Sadek AM, 2013. Test of the accuracy of the computerized glow curve deconvolution algorithm for the analysis of thermoluminescence glow curves.
Nuclear Instruments and Methods in Physics Research A, 712: 56-61.

Sadek AM, Kitis G., 2018. Impact of non-fulfillment of the super position principle on the analysis of thermoluminescence glow-curve. 
Radiation Measurements, 116: 14-23.

Singh TB, Rey L, Gartia RK, 2011. Applications of PeakFit software in thermoluminescence studies. 
Indian Journal of Pure & Applied Physics, 49: 297-302.

Sunta CM, 2015. Unraveling thermoluminescence. Springer India.
}

\note{
Function \link{tgcd} analyzes only thermoluminescence glow curves recorded with linear heating function (LHF) profile.
The model to be optimized should not be underdetermined. This means that the number of data points (\eqn{n_d}) should exceed the number of parameters (\eqn{n_2}).\cr

If it is not \code{NULL} the argument \code{inisPAR}, the procedure used for initializing kinetic parameters by clicking with a mouse will not be triggered.
Similarly, if it is not \code{NULL} the argment \code{inisBG}, the procedure used for initializing background parameters by clicking with a mouse will not be triggered.\cr

The user is advocated to use \code{mwt}, \code{mdt}, and \code{mr} to specify the allowed maximum total half-width, 
the allowed minimum distance, and the allowed minimum resolution respectively to resolve seriously overlapped glow peaks during the trial-and-error protocol.\cr

Adrie J.J. Bos is appreciated for providing the reference glow curves of the GLOCANIN project to test this routine.\cr

Amr M. Sadek is thanked for providing the Matlab code implementing the Lambert W (Wright Omega) function for reference.\cr

George Kitis is appreciated for giving some useful suggestions to improve the program and providing many experimentally measured thermoluminescence glow curves to check the routine.\cr
}

\seealso{
\link{simPeak}; \link{simqOTOR}; \link{savgol}
}

\examples{
# Load the data.
  data(Refglow)
  data(Kitis)

# (1) Deconvolution of a glow curve using 4 peaks (no background subtraction) with 
#     the Wright Omega function using specified initial kinetic parameters.

  knPars <- 
  cbind(c(400, 550, 850, 1600), # Im
        c(1.4, 1.5, 1.6, 2),    # E
        c(420, 460, 480, 510),  # Tm
        c(0.1, 0.1, 0.1, 0.1))  # R

  dd1 <- tgcd(Refglow$x002, npeak=4, model="wo",
              inisPAR=knPars, nstart=10, edit.inis=FALSE)

  head(dd1$comp.sig)
  dd1$pars
  dd1$sp
  dd1$FOM

# (2) Deconvolution of a glow curve using 5 peaks (with background subtraction) with 
#     a mixed-order model using user-supplied intial kinetic and background parameters. 
  knPars <- 
  cbind(c(46829.06, 187942.43, 121876.22, 110390.55, 67978.33), #Im
        c(1.17, 1.14, 1.57, 0.77, 1.31),                        # E
        c(369.86, 400.69, 428.51,  482.41, 537.28),             # Tm
        c(0.75,  0.81, 0.92, 0.001, 0.29))                      # a

  bgPars <- c(1, 10, 100)  # A, B, C.

  dd2 <- tgcd(Kitis$x009, npeak=5, model="m1", subBG=TRUE, 
              inisPAR=knPars, inisBG=bgPars, nstart=10, edit.inis=FALSE)

  dd2$residual
  dd2$SSR
  dd2$R2
  dd2$BGpars
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow curve deconvolution}

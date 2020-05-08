\name{Kitis}
\alias{Kitis}
\docType{data}
\title{
Thermoluminescence glow curves provided by George Kitis
}

\description{
A total of 22 thermoluminescence glow curves measured from various materials provided by George Kitis.
}

\usage{data(Kitis)}

\format{
A list that contains 22 thermoluminescence glow curves.
}

\details{
This object contains 22 thermoluminescence glow curves (\code{x001} to \code{x022}) provided by George Kitis.
x001 (Al2O3:C), x002 (CaF2:Dy), x003 (LBO), x004 (Background), x005 (MgO), x006 (BeO), x007 (CaF2:Tm), x008 (Salt), 
x009 (CaF2:Dy), x010 to x016 (quartz irradiated with dose of 1, 2, 4, 8, 16, 32, and 64 Gy), x017 and x018 (BeO), 
x019 and x020 (Salt), x021 and x022 (MgO). 
}

\examples{
 # Load package "tgcd".
  require(tgcd)

  data(Kitis)
  names(Kitis)
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow curve simulation}
\keyword{Glow curve deconvolution}

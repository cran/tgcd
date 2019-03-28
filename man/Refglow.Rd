\name{Refglow}
\alias{Refglow}
\docType{data}
\title{
Reference glow curves
}

\description{
Synthetic and measured thermoluminescence glow curves from the GLOCANIN project.
}

\usage{data(Refglow)}

\format{
A list that contains 10 thermoluminescence glow curves.
}

\details{
This object contains 10 thermoluminescence glow curves (\code{x001} to \code{x010}) from the GLOCANIN project (Bos et al., 1993, 1994).
}

\references{
Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1993. An intercomparison of glow curve analysis computer programs: I. Synthetic glow curves. Radiation Protection Dosimetry, 47(1-4), 473-477.

Bos AJJ, Piters TM, Gomez Ros JM, Delgado A, 1994. An intercomparison of glow curve analysis computer programs: II. Measured glow curves. Radiation Protection Dosimetry, 51(4): 257-264.
}

\examples{
 # Load package "tgcd".
  require(tgcd)

  data(Refglow)
  names(Refglow)
}

\keyword{package tgcd}
\keyword{Thermoluminescence}
\keyword{Glow curve simulation}
\keyword{Glow curve deconvolution}

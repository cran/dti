\name{dwiRiceBias-methods}
\docType{methods}
\alias{dwiRiceBias}
\alias{dwiRiceBias-methods}
\alias{dwiRiceBias,ANY-method}
\alias{dwiRiceBias,dtiData-method}
\title{ Correction for Rician Bias }
\description{
 Correction for Rician Bias assuming known variance parameter
}
\usage{
  \S4method{dwiRiceBias}{dtiData}(object, sigma=NULL, ncoils=1)
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{sigma}{ Scale parameter that relates the distribution of the signal to a
               \eqn{\chi_{2L}} distribution}
  \item{ncoils}{ number of effective coils in parallel imaging, the related \eqn{\chi}
  distribution has \code{2*ncoils} degrees of freedom.}
}
\value{ 
  An object of class \code{"dtiData"}.
}
\section{Methods}{
\describe{

\item{object = "ANY"}{ Returns a warning. 
}

\item{object = "dtiData"}{
Returns a \code{dtiData} object with bias-corrected image intensities.
}
}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiData}}, 
  \code{\link{dtiTensor-methods}}, 
  \code{\link{dwiMixtensor-methods}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}},
  \code{\linkS4class{dwiMixtensor}},
} 
\keyword{methods}
\keyword{ manip }

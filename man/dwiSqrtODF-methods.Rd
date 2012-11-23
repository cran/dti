\name{dwiSqrtODF-methods}
\docType{methods}
\alias{dwiSqrtODF}
\alias{dwiSqrtODF-methods}
\alias{dwiSqrtODF,ANY-method}
\alias{dwiSqrtODF,dtiData-method}
\title{ Methods for positive definite EAP and ODF estimation in Package \pkg{dti} }
\description{
Compute a positive definite estimate of the Ensemble Average Propagator (EAP)
and Orientation Density Function (ODF) using the approach of Cjeng et. al (2012).
}
\usage{
  \S4method{dwiSqrtODF}{dtiData}(object,what="sqrtODF",order=4,forder=1,lambda=0,D0=1.4e-3)
}
\arguments{
\item{object}{Object of class \code{"dtiData"}}
\item{what}{Character, currently only "sqrtODF" is possible }
\item{order}{Even integer, Order of spherical harmonics approximation.}
\item{forder}{Integer, Order of radial approximation.}
\item{lambda}{Non-negative, Regularization parameter.}
\item{D0}{Numeric vector, grid of diffusivity parameters, typically about 1e-3.}
}
\section{Methods}{
\describe{

\item{\code{signature(object = "ANY")}}{Returns a warning.}

\item{\code{signature(object = "dtiData")}}{
Compute a positive definite estimate of the Ensemble Average Propagator (EAP)
and Orientation Density Function (ODF) using the approach of Cjeng et. al (2012).
}
}}
\author{
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\references{ J. Cheng, T. Jiang and R. Deriche.
           \emph{Nonnegative Definite EAP and ODF Estimation via a Unified Multi-Shell HARDI Reconstruction}, MICCAI 2012. 
}
\seealso{ 
  \code{\link{dtiData}}, 
  \code{\link{readDWIdata}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dwiQball}}
}
\keyword{methods}
\keyword{models  }

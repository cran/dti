\name{getsdofsb-methods}
\docType{methods}
\alias{getsdofsb}
\alias{getsdofsb-methods}
\alias{getsdofsb,ANY-method}
\alias{getsdofsb,dtiData-method}
\title{ Estimate the noise standard deviation }
\description{
 Estimate the noise standard deviation. Uses an assumption that the standard deviation
 is a linear function of the expected mean for image intensities. \code{qA0} and \code{qA1} define quantiles of observed image intensities that define the range of values where this assumption is made.  
}
\usage{
  \S4method{getsdofsb}{dtiData}(object,qA0=.1,qA1=.98,nsb=NULL,level=NULL)
}
\arguments{
  \item{object}{ Object of class \code{"dtiData"}}
  \item{qA0}{ level for lower quantile of image intensities}
  \item{qA1}{ level for upper quantile of image intensities}
  \item{nsb}{ number of diffusion weighted image to use}
  \item{level}{ level for mask}
}
\value{ 
  An object of class \code{"dtiData"} with results in slot \code{sdcoef} in components
  \code{5}: intercept parameter, \code{6}: slope parameter for linear model,
  \code{7}: lower bound (depending on \code{qA0}) and \code{8}: upper bound (depending on \code{qA1}).
}
\section{Methods}{
\describe{

\item{\code{signature(object) = "ANY"}}{
Returns a warning.
}

\item{\code{signature(object) = "dtiData"}}{
Returns a \code{dtiData} object with estimated standard deviation parameters in 
slot \code{sdcoef}.
}
}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{dtiData}}, 
  \code{\link{dwi.smooth-methods}}, 
  \code{\linkS4class{dtiData}}, 
} 
\keyword{methods}
\keyword{ manip }

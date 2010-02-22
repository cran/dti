\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,ANY-method}
\alias{show,dti-method}
\alias{show,dtiData-method}
\alias{show,dtiTensor-method}
\alias{show,dwiMixtensor-method}
\alias{show,dtiIndices-method}
\alias{show,dwiFiber-method}
\title{ Methods for Function `show' in Package `dti' }
\description{
  The function provides information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}, \code{"dtiData"},
 \code{"dtiTensor"}, \code{"dwiMixtensor"}, \code{"dtiIndices"} or \code{"dwiQball"}.
}
\usage{
  \S4method{show}{dti}(object)
}
\arguments{
  \item{object}{ Object of class \code{dtiIndices}, \code{dtiData} or \code{dtiTensor}}
}
\section{Methods}{
\describe{

  \item{x = "ANY"}{ Generic function. }

  \item{x = "dti"}{ The function provides information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}. }
}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiIndices}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}
  \code{\linkS4class{dwiMixtensor}}
  \code{\linkS4class{dwiQball}}
  \code{\linkS4class{dwiFiber}}
}
\keyword{ methods }
\keyword{ utilities }

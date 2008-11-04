\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,ANY-method}
\alias{show,dti-method}
\alias{show,dtiData-method}
\alias{show,dtiTensor-method}
\alias{show,dtiIndices-method}
\title{ Methods for Function `show' in Package `dti' }
\description{
  The function provides information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}, \code{"dtiData"},
 \code{"dtiTensor"} and \code{"dtiIndices"}.
}
\usage{
  \S4method{show}{dti}(object)
}
\arguments{
  \item{object}{ Object of class \code{dtiIndices}, \code{dtiData} or \code{dtiTensor}}
}
\section{Methods}{
\describe{

  \item{x = "ANY"}{ Generic function: see \code{\link[base]{show}}. }

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
}
\keyword{ methods }
\keyword{ utilities }

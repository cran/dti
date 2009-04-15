\name{print-methods}
\docType{methods}
\alias{print-methods}
\alias{print,ANY-method}
\alias{print,dwi-method}
\alias{print,dtiData-method}
\alias{print,dtiTensor-method}
\alias{print,dtiIndices-method}
\alias{print,dwiQball-method}
\title{ Methods for Function `print' in Package `dti' }
\description{
  The function provides information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}, \code{"dtiData"},
 \code{"dtiTensor"} and \code{"dtiIndices"}.
}
\usage{
  \S4method{print}{dwi}(x)
}
\arguments{
  \item{x}{ Object of class \code{dtiIndices}, \code{dtiData}, \code{dtiTensor} or \code{dwiQball}}
}
\section{Methods}{
\describe{

  \item{x = "ANY"}{ Generic function: see \code{\link[base]{print}}. }

  \item{x = "dwi"}{ The function provides information on data dimensions, data source and existing slot-names for objects of class \code{"dwi"}. }
}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiIndices}}, 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}
  \code{\linkS4class{dwiQball}}
}
\keyword{ methods }
\keyword{ utilities }

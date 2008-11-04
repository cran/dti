\name{print-methods}
\docType{methods}
\alias{print-methods}
\alias{print,ANY-method}
\alias{print,dti-method}
\alias{print,dtiData-method}
\alias{print,dtiTensor-method}
\alias{print,dtiIndices-method}
\title{ Methods for Function `print' in Package `dti' }
\description{
  The function provides information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}, \code{"dtiData"},
 \code{"dtiTensor"} and \code{"dtiIndices"}.
}
\usage{
  \S4method{print}{dti}(x)
}
\arguments{
  \item{x}{ Object of class \code{dtiIndices}, \code{dtiData} or \code{dtiTensor}}
}
\section{Methods}{
\describe{

  \item{x = "ANY"}{ Generic function: see \code{\link[base]{print}}. }

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

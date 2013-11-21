\name{print-methods}
\docType{methods}
\alias{print-methods}
\alias{print,ANY-method}
\alias{print,dwi-method}
\alias{print,dtiData-method}
\alias{print,dtiTensor-method}
\alias{print,dkiIndices-method}
\alias{print,dkiTensor-method}
\alias{print,dwiMixtensor-method}
\alias{print,dtiIndices-method}
\alias{print,dwiQball-method}
\alias{print,dwiFiber-method}
\title{ Methods for Function `print' in Package `dti' }
\description{
  The function provides information on data dimensions, data source and existing slot-names for objects of class \code{"dti"}, \code{"dtiData"},
 \code{"dtiTensor"}, \code{"dtiIndices"}, \code{"dkiIndices"}, \code{"dkiTensor"}, \code{"dwiMixtensor"},  \code{"dwiQball"} and \code{"dwiFiber"}.
}
\usage{
  \S4method{print}{dwi}(x)
}
\arguments{
  \item{x}{ Object of class \code{"dtiIndices"}, \code{"dtiData"}, \code{"dtiTensor"}, \code{"dkiIndices"}, \code{"dkiTensor"}, \code{"dwiFiber"}, \code{"dwiMixtensor"} or \code{"dwiQball"}}
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
  \code{\linkS4class{dwiMixtensor}}
  \code{\linkS4class{dwiQball}}
  \code{\linkS4class{dwiFiber}}
}
\keyword{ methods }
\keyword{ utilities }

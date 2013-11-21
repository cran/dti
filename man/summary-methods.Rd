\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,dwi-method}
\alias{summary,dtiData-method}
\alias{summary,dtiTensor-method}
\alias{summary,dkiTensor-method}
\alias{summary,dwiMixtensor-method}
\alias{summary,dtiIndices-method}
\alias{summary,dkiIndices-method}
\alias{summary,dwiQball-method}
\alias{summary,dwiFiber-method}
\title{ Methods for Function `summary' in Package `dti' }
\description{
  The method provides summary information for objects of class \code{"dti"}.
}
\usage{
  \S4method{summary}{dwi}(object, ...)
}
\arguments{
  \item{object}{ Object of class \code{"dti"}, \code{"dtiData"}, \code{"dtiTensor"}, \code{"dwiMixtensor"}, \code{"dtiIndices"}, \code{"dkiIndices"}, \code{"dkiTensor"},  \code{"dwiQball"} or \code{"dwiFiber"}. }
  \item{\dots}{ Additional arguments in \code{\dots} are passed to function \code{quantile}, e.g. argument \code{probs} may be specified here.}
}
\section{Methods}{
\describe{

  \item{object = "ANY"}{ Generic function: see \code{\link[base]{summary}}. }

  \item{object = "dwi"}{  The function provides summary information for objects of class \code{"dwi"}, \code{"dtiData"}, \code{"dtiTensor"}, \code{"dwiMixtensor"}, \code{"dtiIndices"}, \code{"dkiIndices"}, \code{"dkiTensor"}, \code{"dwiQball"} 
 and , \code{"dwiFiber"} }
}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiIndices}}, ]
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}
  \code{\linkS4class{dwiMixtensor}}
  \code{\linkS4class{dwiQball}}
  \code{\linkS4class{dwiFiber}}
}
\keyword{ methods }
\keyword{ utilities }

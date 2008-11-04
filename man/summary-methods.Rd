\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,dti-method}
\alias{summary,dtiData-method}
\alias{summary,dtiTensor-method}
\alias{summary,dtiIndices-method}
\title{ Methods for Function `summary' in Package `dti' }
\description{
  The method provides summary information for objects of class \code{"dti"}.
}
\usage{
  \S4method{summary}{dti}(object, ...)
}
\arguments{
  \item{object}{ Object of class \code{"dti"}, \code{"dtiData"}, \code{"dtiTensor"}, or \code{"dtiIndices"}. }
  \item{\dots}{ Additional arguments in \code{\dots} are passed to function \code{quantile}, e.g. argument \code{probs} may be specified here.}
}
\section{Methods}{
\describe{

  \item{object = "ANY"}{ Generic function: see \code{\link[base]{summary}}. }

  \item{object = "dti"}{  The function provides summary information for objects of class \code{"dti"}, \code{"dtiData"}, \code{"dtiTensor"} and \code{"dtiIndices"}. }
}}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiIndices}}, ]
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}
}
\keyword{ methods }
\keyword{ utilities }

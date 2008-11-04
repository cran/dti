\name{dtiIndices-methods}
\docType{methods}
\alias{dtiIndices}
\alias{dtiIndices-methods}
\alias{dtiIndices,ANY-method}
\alias{dtiIndices,dtiTensor-method}
\title{ Methods for Function `dtiIndices' in Package `dti' }
\description{
  The method creates estimates of the fractional anisotropy (FA) and relative anisotropy (RA) indices, the main directions of anisotropy and several statistics used for visualization.
}
\usage{
  \S4method{dtiIndices}{dtiTensor}(object, which)
}
\arguments{
  \item{object}{ Object of class \code{"dtiTensor"} }
  \item{which}{ Indices should be created, currently not implemented. }
}
\value{
  An object of class \code{"dtiIndices"}.
}
\section{Methods}{
\describe{

\item{obj = "ANY"}{ Returns a warning. }

\item{obj = "dtiTensor"}{ Estimate tensor indices like trace, fractional and geodesic anisotropy, main diffusion direction and shape parameters. }

}
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\link{medinria}}, 
  \code{\link{dtiTensor-methods}}, 
  \code{\linkS4class{dtiTensor}}, 
  \code{\linkS4class{dtiIndices}}
}
\examples{
  \dontrun{demo(dti_art)}
}
\keyword{ manip }
\keyword{ models }

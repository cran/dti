\name{dti-class}
\docType{class}
\alias{dti-class}
\alias{plot,dti-method}
\alias{show,dti-method}

\title{Class "dti" }
\description{ Basic description of Diffusion Weighted Image (DWI) data}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("dti", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} usually empty }
    \item{\code{btb}:}{Object of class \code{"matrix"} matrix of dimension \code{c(6,ngrad)} obtained from gradient directions.}
    \item{\code{ngrad}:}{Object of class \code{"integer"} number of gradients (including zero gradients) }
    \item{\code{s0ind}:}{Object of class \code{"integer"} index of zero gradients within sequence \code{1:ngrad} }
    \item{\code{replind}:}{Object of class \code{"integer"} index (identifier) of unique 
           gradient directions. Used to charactreize replications in the gradient design by identical indices. length is \code{ngrad}}
    \item{\code{ddim}:}{Object of class \code{"integer"} dimension of original image cubes. Integer vector of length 3. }
    \item{\code{ddim0}:}{Object of class \code{"integer"} dimension of subcube defined by \code{xind}, \code{yind} and \code{zind}. }
    \item{\code{xind}:}{Object of class \code{"integer"} index for subcube definition in x-direction  }
    \item{\code{yind}:}{Object of class \code{"integer"} index for subcube definition in y-direction }
    \item{\code{zind}:}{Object of class \code{"integer"} index for subcube definition in z-direction }
    \item{\code{voxelext}:}{Object of class \code{"numeric"} voxel extensions in x-, y- and z-direction. vector of length 3. }
    \item{\code{source}:}{Object of class \code{"character"} name of the imgfile used to create the data. }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(x = "dti")}: not yet implemented }
    \item{show}{\code{signature(object = "dti")}: not yet implemented }
	 }
}
\references{ 
\url{http://www.wias-berlin.de/projects/matheon_a3/} 
}
\author{ Karsten Tabelow tabelow@wias-berlin.de, J\"org Polzehl polzehl@wias-berlin.de }
\seealso{\code{\link{dtiData}}, \code{\linkS4class{dtiData}}, \code{\linkS4class{dtiTensor}}, \code{\linkS4class{dtiIndices}}
}
\examples{
showClass("dti")
}
\keyword{classes}

\name{extract-methods}
\docType{methods}
\alias{extract}
\alias{extract-methods}
\alias{extract,ANY-method}
\alias{extract,dtiData-method}
\alias{extract,dtiIndices-method}
\alias{extract,dtiTensor-method}
\alias{[-methods}
\alias{[,ANY-method}
\alias{[,dtiData-method}
\alias{[,dtiIndices-method}
\alias{[,dtiTensor-method}
\title{ Methods for Function `extract' and `[' in Package `dti' }
\description{
  The methods extract and/or compute specified statistics from object of class \code{"dtiData"}, \code{"dtiTensor"}, and \code{"dtiIndices"}. This can be restricted to a subset of voxel.
}
\usage{
  \S4method{extract}{dtiData}(x, what="data", xind, yind, zind)
  \S4method{extract}{dtiTensor}(x, what="tensor", xind, yind, zind)
  \S4method{extract}{dtiIndices}(x, what=c("fa","andir"), xind, yind, zind)
  \S4method{[}{dtiData}(x, i, j, k, drop=FALSE)
  \S4method{[}{dtiTensor}(x, i, j, k, drop=FALSE)
  \S4method{[}{dtiIndices}(x, i, j, k, drop=FALSE)
}
\arguments{
  \item{x}{ Object of class \code{dti} }
  \item{i}{ vector of x-coordinates, defaults to whole range. }
  \item{j}{ vector of y-coordinates, defaults to whole range. }
  \item{k}{ vector of z-coordinates, defaults to whole range. }
  \item{xind}{ vector of x-coordinates, defaults to whole range. }
  \item{yind}{ vector of y-coordinates, defaults to whole range. }
  \item{zind}{ vector of z-coordinates, defaults to whole range. }
  \item{what}{ Statistic to extract. See Methods Section for details. }
  \item{drop}{ unused. }
}
\section{Methods}{
\describe{
  The generic extract function \code{"["} does what it is expected to do: it extracts parts of the object specified by \code{i}, \code{j}, and \code{k}.

  \item{x = "ANY"}{ Returns a warning for \code{extract}. Generic funtion for \code{"["} returns an object of same class with data clipped to the indices specified in arguments \code{i}, \code{j} and \code{k}.}

  \item{x = "dtiData"}{ Extraction of squared gradient matrix (btb) or of S0 (S0), Sb (Sb) or all images (Si) restricted to the cube defined by arguments \code{i}, \code{j} and \code{k}.}

  \item{x = "dtiIndices"}{ Returns an array containing the specified statistics, i.e. fractional anisotropy (fa), geodesic anisotropy (ga), mean diffusivity (md), main direction of anisotropy (andir) and/or shape parameters (bary), as specified in argument \code{what}. Information is extracted for voxel within the cube defined by \code{xind}, \code{yind}, and \code{zind}. }

  \item{x = "dtiTensor"}{ Returns an array containing the specified statistics, i.e. fractional anisotropy (fa), geodesic anisotropy (ga), mean diffusivity (md), eigenvalues (evalues), main direction of anisotropy (andir), the tensor (tensor) the estimated S0 image (s0) and/or the mask used to restrict computations (mask), as specified in argument \code{what}. Information is extracted for voxel within the  cube defined by arguments \code{xind}, \code{yind} and \code{zind}. }
}
}
\value{
  For function \code{extract} a list with components carrying the names of the options specified in
  argument \code{what}.  For indexing function, the cutted object.
}
\author{ 
  Karsten Tabelow \email{tabelow@wias-berlin.de}\cr 
  J\"org Polzehl \email{polzehl@wias-berlin.de}
}
\seealso{ 
  \code{\linkS4class{dtiData}}, 
  \code{\linkS4class{dtiTensor}}, 
  \code{\linkS4class{dtiIndices}}  
}
\keyword{ methods }
\keyword{ manip }
